from __future__ import annotations

import json
import math
import os
import socket
import subprocess
import sys
import threading
import time
from collections import deque

# Configuración de salida UTF-8 para evitar errores de codificación en Windows
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
if hasattr(sys.stderr, "reconfigure"):
    try:
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

# --- MÉTODOS NUMÉRICOS ---
from orbital_metrics import romberg_integration
from reconstruction import neville_interpolate

# --------------------------------------------------------------------------- #
# Configuración y Constantes
# --------------------------------------------------------------------------- #
UDP_HOST, UDP_PORT = "127.0.0.1", 5005
PORT_SAT = 5006

BG = "#0B1220"
SURFACE = "#111B2F"
LINE = "#22304D"
TEXT = "#E7EDFB"
MUTED = "#8A97B8"
ICE = "#7CC4FF"
AMBER = "#FFB454"
ALERT = "#FF6B6B"
OK = "#5FD4A0"

# Ángulos orbitales de Kepler (Inclinación tipo ISS)
INCLINATION = math.radians(51.6)
RAAN = math.radians(30.0)
ARG_PERIAPSIS = math.radians(45.0)
EARTH_RADIUS = 6371.0

st.set_page_config(page_title="Gemelo Digital · MÉTODOS NUMÉRICOS", page_icon="🛰️", layout="wide")


# --------------------------------------------------------------------------- #
# Transformación Orbital 2D → 3D
# --------------------------------------------------------------------------- #
def orbital_to_3d(x_orb: float, y_orb: float) -> tuple[float, float, float]:
    """Transforma coordenadas orbitales planas a coordenadas inerciales 3D.
    
    Aplica la matriz de rotación compuesta: R = R_z(-Ω) · R_x(-i) · R_z(-ω)
    """
    w = ARG_PERIAPSIS
    i = INCLINATION
    omega = RAAN

    # Rotación en el plano orbital por argumento del periapsis
    x1 = x_orb * math.cos(w) - y_orb * math.sin(w)
    y1 = x_orb * math.sin(w) + y_orb * math.cos(w)

    # Inclinación orbital
    x2 = x1
    y2 = y1 * math.cos(i)
    z2 = y1 * math.sin(i)

    # Rotación por longitud del nodo ascendente (RAAN)
    x_3d = x2 * math.cos(omega) - y2 * math.sin(omega)
    y_3d = x2 * math.sin(omega) + y2 * math.cos(omega)
    z_3d = z2

    return x_3d, y_3d, z_3d


@st.cache_data
def make_earth_wireframe_cached(
    radius: float = EARTH_RADIUS, n_meridians: int = 16, n_parallels: int = 12, pts: int = 50
) -> list:
    """Genera la geometría fija de la Tierra y la almacena en caché."""
    traces = []

    # Meridianos
    for j in range(n_meridians):
        lon = j * (2.0 * math.pi / n_meridians)
        mx, my, mz = [], [], []
        for k in range(pts + 1):
            lat = -math.pi / 2 + k * (math.pi / pts)
            mx.append(radius * math.cos(lat) * math.cos(lon))
            my.append(radius * math.cos(lat) * math.sin(lon))
            mz.append(radius * math.sin(lat))
        traces.append(go.Scatter3d(
            x=mx, y=my, z=mz, mode="lines",
            line=dict(color="#1a5276", width=1),
            hoverinfo="skip", showlegend=False,
        ))

    # Paralelos
    for j in range(1, n_parallels):
        lat = -math.pi / 2 + j * (math.pi / n_parallels)
        px, py, pz = [], [], []
        r_lat = radius * math.cos(lat)
        z_lat = radius * math.sin(lat)
        for k in range(pts + 1):
            lon = k * (2.0 * math.pi / pts)
            px.append(r_lat * math.cos(lon))
            py.append(r_lat * math.sin(lon))
            pz.append(z_lat)
        is_eq = abs(lat) < 0.15
        traces.append(go.Scatter3d(
            x=px, y=py, z=pz, mode="lines",
            line=dict(color="#2e86c1" if is_eq else "#1a5276", width=2.5 if is_eq else 1),
            hoverinfo="skip", showlegend=False,
        ))

    # Malla volumétrica translúcida
    sphere_x, sphere_y, sphere_z = [], [], []
    n_u, n_v = 18, 18
    for ui in range(n_u + 1):
        u = ui * (2.0 * math.pi / n_u)
        for vi in range(n_v + 1):
            v = vi * (math.pi / n_v)
            sphere_x.append(radius * math.sin(v) * math.cos(u))
            sphere_y.append(radius * math.sin(v) * math.sin(u))
            sphere_z.append(radius * math.cos(v))

    i_list, j_list, k_list = [], [], []
    for ui in range(n_u):
        for vi in range(n_v):
            p0 = ui * (n_v + 1) + vi
            p1 = p0 + 1
            p2 = (ui + 1) * (n_v + 1) + vi
            p3 = p2 + 1
            i_list.extend([p0, p0])
            j_list.extend([p1, p2])
            k_list.extend([p2, p3])

    traces.append(go.Mesh3d(
        x=sphere_x, y=sphere_y, z=sphere_z,
        i=i_list, j=j_list, k=k_list,
        color="#0d3b66", opacity=0.20,
        hoverinfo="skip", showlegend=False, flatshading=True,
    ))

    # Anillo ecuatorial de referencia inercial (Z = 0)
    eq_x = [(radius + 800) * math.cos(t * 2 * math.pi / pts) for t in range(pts + 1)]
    eq_y = [(radius + 800) * math.sin(t * 2 * math.pi / pts) for t in range(pts + 1)]
    eq_z = [0.0] * (pts + 1)
    traces.append(go.Scatter3d(
        x=eq_x, y=eq_y, z=eq_z, mode="lines",
        line=dict(color="rgba(124, 196, 255, 0.3)", width=1, dash="dash"),
        name="Plano Ecuatorial (Z=0)", showlegend=True, hoverinfo="skip"
    ))

    return traces


def make_satellite_3d_mesh(
    pos: tuple[float, float, float],
    prev_pos: tuple[float, float, float] | None = None,
    scale: float = 600.0  # Tamaño visual en km para que sea visible frente a los 6371 km de la Tierra
) -> list:
    """Genera la geometría 3D completa de un satélite (cuerpo cúbico + paneles solares).
    
    Alinea el cuerpo y las alas según los vectores de velocidad y radio orbital.
    """
    x0, y0, z0 = pos

    # 1. Base ortonormal local para orientar el satélite
    # Si tenemos posición previa, calculamos el vector velocidad tangente
    if prev_pos is not None:
        vx = x0 - prev_pos[0]
        vy = y0 - prev_pos[1]
        vz = z0 - prev_pos[2]
        v_norm = math.sqrt(vx**2 + vy**2 + vz**2)
        if v_norm > 1e-6:
            tx, ty, tz = vx / v_norm, vy / v_norm, vz / v_norm
        else:
            tx, ty, tz = 1.0, 0.0, 0.0
    else:
        tx, ty, tz = 1.0, 0.0, 0.0

    # Vector normal radial (hacia afuera de la Tierra)
    r_norm = math.sqrt(x0**2 + y0**2 + z0**2)
    rx, ry, rz = (x0 / r_norm, y0 / r_norm, z0 / r_norm) if r_norm > 0 else (0.0, 0.0, 1.0)

    # Vector transversal (alas/paneles solares hacia los lados): nx = T x R
    nx = ty * rz - tz * ry
    ny = tz * rx - tx * rz
    nz = tx * ry - ty * rx
    n_norm = math.sqrt(nx**2 + ny**2 + nz**2)
    if n_norm > 1e-6:
        nx, ny, nz = nx / n_norm, ny / n_norm, nz / n_norm
    else:
        nx, ny, nz = 0.0, 1.0, 0.0

    # Re-ortogonalizar el vector vertical del chasis: up = N x T
    ux = ny * tz - nz * ty
    uy = nz * tx - nx * tz
    uz = nx * ty - ny * tx

    traces = []

    # 2. CUERPO PRINCIPAL (Cubo dorado / dorado metálico)
    # Definimos los 8 vértices en coordenadas locales y los proyectamos
    s = scale * 0.4
    vertices_body = []
    for dx in (-s, s):
        for dy in (-s, s):
            for dz in (-s, s):
                # Rotación al marco inercial + traslación a la posición orbital
                px = x0 + dx * tx + dy * nx + dz * ux
                py = y0 + dx * ty + dy * ny + dz * uy
                pz = z0 + dx * tz + dy * nz + dz * uz
                vertices_body.append((px, py, pz))

    bx = [v[0] for v in vertices_body]
    by = [v[1] for v in vertices_body]
    bz = [v[2] for v in vertices_body]

    # Triangulación de las 6 caras del cubo (12 triángulos)
    i_body = [0, 0, 4, 4, 0, 0, 2, 2, 1, 1, 3, 3]
    j_body = [2, 3, 6, 7, 4, 5, 6, 7, 5, 7, 7, 6]
    k_body = [1, 2, 5, 6, 1, 4, 3, 6, 3, 5, 2, 7]

    traces.append(go.Mesh3d(
        x=bx, y=by, z=bz,
        i=i_body, j=j_body, k=k_body,
        color="#F5B041",  # Dorado satelital / capa térmica de kapton
        opacity=1.0,
        flatshading=True,
        name="Cuerpo Satélite",
        showlegend=True,
        hoverinfo="skip"
    ))

    # 3. PANELES SOLARES (Alas rectangulares cian/azul reflectante a los lados)
    wing_span = scale * 1.8
    wing_chord = scale * 0.5
    gap = s * 1.1

    def make_solar_panel(side: float):
        # side = +1 para ala derecha, -1 para ala izquierda
        start_span = gap * side
        end_span = (gap + wing_span) * side

        w_pts = [
            # 4 vértices del rectángulo plano
            (x0 + (-wing_chord)*tx + start_span*nx, y0 + (-wing_chord)*ty + start_span*ny, z0 + (-wing_chord)*tz + start_span*nz),
            (x0 + (wing_chord)*tx  + start_span*nx, y0 + (wing_chord)*ty  + start_span*ny, z0 + (wing_chord)*tz  + start_span*nz),
            (x0 + (wing_chord)*tx  + end_span*nx,   y0 + (wing_chord)*ty  + end_span*ny,   z0 + (wing_chord)*tz  + end_span*nz),
            (x0 + (-wing_chord)*tx + end_span*nx,   y0 + (-wing_chord)*ty + end_span*ny,   z0 + (-wing_chord)*tz + end_span*nz),
        ]
        wx = [p[0] for p in w_pts]
        wy = [p[1] for p in w_pts]
        wz = [p[2] for p in w_pts]

        return go.Mesh3d(
            x=wx, y=wy, z=wz,
            i=[0, 0], j=[1, 2], k=[2, 3],
            color="#00D2FF",  # Azul eléctrico celdas fotovoltaicas
            opacity=0.9,
            flatshading=True,
            showlegend=False,
            hoverinfo="skip"
        )

    # Panel solar izquierdo y panel solar derecho
    traces.append(make_solar_panel(side=1.0))
    traces.append(make_solar_panel(side=-1.0))

    return traces



# --------------------------------------------------------------------------- #
# Backend: Hub de Telemetría Thread-Safe
# --------------------------------------------------------------------------- #
class OrbitalTelemetryHub:
    def __init__(self, host: str, port: int, maxlen: int = 2000):
        self.host, self.port = host, port
        self._buf: deque[dict] = deque(maxlen=maxlen)
        self.interpolations: list[dict] = []
        self._lock = threading.Lock()
        self._stop = threading.Event()

        self.expected_id = 0
        self.dropped_packets = 0
        self.bind_error: str | None = None

        threading.Thread(target=self._rx_loop, daemon=True, name="udp-rx-3d").start()

    def _rx_loop(self) -> None:
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        try:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            sock.bind((self.host, self.port))
        except OSError as exc:
            self.bind_error = f"No se pudo abrir {self.host}:{self.port} ({exc})"
            return

        sock.settimeout(0.5)
        last_ping = 0.0

        while not self._stop.is_set():
            now = time.time()
            if now - last_ping > 1.0:
                try:
                    sock.sendto(b"GROUND_READY", (self.host, PORT_SAT))
                except Exception:
                    pass
                last_ping = now

            try:
                data, _ = sock.recvfrom(16384)
            except (socket.timeout, ConnectionResetError):
                continue
            except OSError:
                break

            try:
                pkt = json.loads(data.decode("utf-8"))
                pkt["rx"] = time.time()
            except (ValueError, UnicodeDecodeError):
                continue

            with self._lock:
                p_id = pkt["id"]

                # Detección de reinicio de misión
                if p_id == 0 and self.expected_id > 0:
                    self._buf.clear()
                    self.interpolations.clear()
                    self.dropped_packets = 0
                    self.expected_id = 0

                # Detección de salto en secuencia (Activación de Neville para cada paquete perdido)
                if p_id > self.expected_id:
                    drops = p_id - self.expected_id
                    self.dropped_packets += drops

                    if len(self._buf) >= 3:
                        # Nodos de soporte que acotan la brecha (anteriores + actual recibido)
                        nodes_x = [(p["timestamp"], p["x"]) for p in list(self._buf)[-3:]] + [(pkt["timestamp"], pkt["x"])]
                        nodes_y = [(p["timestamp"], p["y"]) for p in list(self._buf)[-3:]] + [(pkt["timestamp"], pkt["y"])]

                        last_pkt = self._buf[-1]
                        t_prev = last_pkt["timestamp"]
                        id_prev = self.expected_id - 1

                        # Interpolar cada uno de los paquetes perdidos dentro de la brecha
                        for lost_id in range(self.expected_id, p_id):
                            frac = (lost_id - id_prev) / (p_id - id_prev)
                            t_target = t_prev + frac * (pkt["timestamp"] - t_prev)

                            res_x = neville_interpolate(nodes_x, t_target)
                            res_y = neville_interpolate(nodes_y, t_target)

                            self.interpolations.append({
                                "id": lost_id,
                                "timestamp": t_target,
                                "x": res_x.estimated_value,
                                "y": res_y.estimated_value,
                                "error": max(res_x.estimated_error, res_y.estimated_error),
                                "nodes_x": nodes_x,
                                "nodes_y": nodes_y,
                                "tableau_x": res_x.tableau,
                                "tableau_y": res_y.tableau,
                            })

                self._buf.append(pkt)
                self.expected_id = p_id + 1

        sock.close()

    def reset(self) -> None:
        """Reinicia el buffer y métricas para una nueva misión."""
        with self._lock:
            self._buf.clear()
            self.interpolations.clear()
            self.dropped_packets = 0
            self.expected_id = 0

    def snapshot(self) -> tuple[list[dict], list[dict], int, int]:
        with self._lock:
            return list(self._buf), list(self.interpolations), self.dropped_packets, self.expected_id



@st.cache_resource
def get_hub() -> OrbitalTelemetryHub:
    return OrbitalTelemetryHub(UDP_HOST, UDP_PORT)


# Importamos la función de simulación directamente del script del satélite
try:
    from satellite_uplink import run_satellite
except ImportError:
    run_satellite = None


class SatelliteLauncher:
    """Gestor del hilo worker del satélite para ejecución en background."""
    def __init__(self):
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()

    def launch(self, hub: OrbitalTelemetryHub | None = None) -> bool:
        """Lanza run_satellite en un hilo daemon. Retorna True si inició, False si ya corre."""
        with self._lock:
            if self._thread is not None and self._thread.is_alive():
                return False  # Misión en curso

            if run_satellite is None:
                st.error("No se encontró 'satellite_uplink.py' en la carpeta actual.")
                return False

            if hub is not None:
                hub.reset()

            def _worker():
                try:
                    run_satellite()
                except Exception as exc:
                    print(f"Error en worker del satélite: {exc}", file=sys.stderr)

            # Inicia el emisor orbital como hilo secundario
            self._thread = threading.Thread(
                target=_worker,
                daemon=True,
                name="sat-worker-thread"
            )
            self._thread.start()
            return True

    @property
    def is_running(self) -> bool:
        with self._lock:
            return self._thread is not None and self._thread.is_alive()
@st.cache_resource
def get_launcher() -> SatelliteLauncher:
    return SatelliteLauncher()


# --------------------------------------------------------------------------- #
# Estilos CSS
# --------------------------------------------------------------------------- #
CSS = f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Sora:wght@400;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap');
html, body, [class*="st-"], .stApp {{ font-family: 'Sora', system-ui, sans-serif; }}
.stApp {{ background: radial-gradient(1200px 600px at 85% -10%, rgba(124,196,255,.10), transparent 60%), {BG}; color: {TEXT}; }}
.kpi {{ background:{SURFACE}; border:1px solid {LINE}; border-radius:14px; padding:1rem 1.1rem; position:relative; overflow:hidden; }}
.kpi::before {{ content:""; position:absolute; left:0; top:0; bottom:0; width:3px; background:var(--c); }}
.kpi .lbl {{ color:{MUTED}; font-size:.82rem; }}
.kpi .val {{ font-family:'IBM Plex Mono', monospace; font-size:1.80rem; font-weight:500; margin:.15rem 0; letter-spacing:-.02em; }}
.kpi .val small {{ font-size:.9rem; color:{MUTED}; margin-left:.25rem; }}
.kpi .sub {{ font-size:.8rem; color:var(--c); }}
.empty {{ padding:2.5rem 1rem; text-align:center; color:{MUTED}; }}
</style>
"""


def kpi(label: str, value: str, unit: str, state: str, color: str) -> str:
    return (
        f'<div class="kpi" style="--c:{color}">'
        f'<div class="lbl">{label}</div>'
        f'<div class="val">{value}<small>{unit}</small></div>'
        f'<div class="sub">{state}</div>'
        f'</div>'
    )


# --------------------------------------------------------------------------- #
# Pizarra: desarrollo numérico en vivo
# --------------------------------------------------------------------------- #
def _newton_log_line(pkt: dict) -> str:
    kepler = pkt.get("kepler") or {}
    steps = kepler.get("steps") or []
    ident = f"PKT {int(pkt['id']):02d}"
    if not steps:
        return f"{ident} · posición recibida x = {pkt['x']:.2f} km, y = {pkt['y']:.2f} km"
    last = steps[-1]
    return (
        f"{ident} · Kepler–Newton · M = {kepler['M']:.4f} rad · "
        f"E{last['k'] - 1} = {last['E']:.6f} · "
        f"f = {last['f']:.2e} · f' = {last['fp']:.6f} · "
        f"delta = {last['delta']:.2e} · E{last['k']} = {last['E_next']:.6f} · "
        f"x = {pkt['x']:.2f} km · y = {pkt['y']:.2f} km"
    )


def _neville_log_line(item: dict) -> str:
    return (
        f"Hueco PKT {int(item['id']):02d} · Neville · t* = {item['timestamp']:.4f} rad · "
        f"Q(X) = {item['x']:.2f} km · Q(Y) = {item['y']:.2f} km · "
        f"cota = {item['error']:.4f} km"
    )


def render_direct_ticker(packets: list[dict], interps: list[dict]) -> None:
    """Línea de cálculo que cambia con el paquete que acaba de llegar."""
    pkt = packets[-1]
    kepler = pkt.get("kepler") or {}
    steps = kepler.get("steps") or []
    with st.container(border=True):
        st.markdown("##### En directo")
        if steps:
            last = steps[-1]
            k = int(last["k"])
            st.latex(
                rf"E_{{{k}}}={last['E_next']:.8f}"
                rf"=E_{{{k - 1}}}-\left(\frac{{{last['f']:.3e}}}{{{last['fp']:.6f}}}\right)"
                rf"={last['E']:.8f}-({last['delta']:.3e})"
            )
            if "a" in kepler and "E" in kepler:
                a = float(kepler["a"])
                e = float(kepler["e"])
                big_e = float(kepler["E"])
                x_calc = a * (math.cos(big_e) - e)
                y_calc = a * math.sqrt(max(0.0, 1.0 - e * e)) * math.sin(big_e)
                st.latex(
                    rf"x=a(\cos E-e)={x_calc:.2f}\ \mathrm{{km}},\quad "
                    rf"y=a\sqrt{{1-e^2}}\sin E={y_calc:.2f}\ \mathrm{{km}}"
                )
            st.caption(
                f"Paquete {int(pkt['id']):02d} · {len(steps)} iteraciones · "
                f"residuo {float(kepler.get('residual', 0.0)):.2e}"
            )
        else:
            st.markdown(
                f"Paquete {int(pkt['id']):02d} en antena · "
                f"x = {pkt['x']:.2f} km · y = {pkt['y']:.2f} km"
            )

        if interps:
            item = interps[-1]
            st.latex(
                rf"Q(t^*)\ \mathrm{{en}}\ t^*={item['timestamp']:.4f}"
                rf"\rightarrow X={item['x']:.2f}\ \mathrm{{km}},\ "
                rf"Y={item['y']:.2f}\ \mathrm{{km}},\ "
                rf"|Q_{{n,n}}-Q_{{n,n-1}}|={item['error']:.4f}\ \mathrm{{km}}"
            )
            st.caption(f"Último hueco reconstruido: paquete {int(item['id']):02d}.")


def _render_newton_process(pkt: dict) -> None:
    kepler = pkt.get("kepler") or {}
    steps = kepler.get("steps") or []
    st.markdown(f"**Kepler–Newton · paquete {int(pkt['id']):02d}**")
    if not steps or "M" not in kepler:
        st.caption("Este paquete llegó sin la traza de Newton. Lanza de nuevo la misión para ver las iteraciones.")
        return

    m_anom = float(kepler["M"])
    ecc = float(kepler["e"])
    axis = float(kepler["a"])
    big_e = float(kepler["E"])
    st.latex(r"f(E)=E-e\sin(E)-M=0,\qquad f'(E)=1-e\cos(E)")
    st.latex(
        rf"M={m_anom:.6f},\quad e={ecc:.6f},\quad a={axis:.2f}\ \mathrm{{km}}"
    )
    if ecc < 0.8:
        st.caption("Semilla E₀ = M, porque la excentricidad es menor que 0.8.")
    else:
        st.caption("Semilla ajustada: excentricidad alta, E₀ parte de π o de M + e sin M.")

    bullets = []
    for step in steps:
        k = int(step["k"])
        bullets.append(
            f"- Iteración {k}: "
            f"$E_{{{k}}}={step['E_next']:.8f}"
            f"=E_{{{k - 1}}}-\\left(\\dfrac{{{step['f']:.3e}}}{{{step['fp']:.6f}}}\\right)"
            f"={step['E']:.8f}-({step['delta']:.3e})$"
        )
    st.markdown("\n".join(bullets))

    frame = pd.DataFrame(
        [
            {
                "Iteración": int(step["k"]),
                "E actual": f"{step['E']:.8f}",
                "f(E)": f"{step['f']:.3e}",
                "f'(E)": f"{step['fp']:.8f}",
                "Corrección": f"{step['delta']:.3e}",
                "E siguiente": f"{step['E_next']:.8f}",
            }
            for step in steps
        ]
    )
    st.table(frame, hide_index=True, alt="Iteraciones de Newton-Raphson del paquete en antena")

    last = steps[-1]
    stopped = abs(float(last["delta"])) < 1e-10
    if stopped:
        st.caption(f"Parada: |ΔE| = {abs(float(last['delta'])):.2e} < 10⁻¹⁰. Residuo |f(E)| = {float(kepler['residual']):.2e}.")
    else:
        st.caption(f"Última corrección |ΔE| = {abs(float(last['delta'])):.2e}. Residuo |f(E)| = {float(kepler['residual']):.2e}.")

    x_calc = axis * (math.cos(big_e) - ecc)
    y_calc = axis * math.sqrt(max(0.0, 1.0 - ecc * ecc)) * math.sin(big_e)
    st.latex(
        rf"x=a(\cos E-e)={axis:.2f}\left(\cos({big_e:.6f})-{ecc:.6f}\right)={x_calc:.2f}\ \mathrm{{km}}"
    )
    st.latex(
        rf"y=a\sqrt{{1-e^2}}\sin E={axis:.2f}\sqrt{{1-({ecc:.6f})^2}}\sin({big_e:.6f})={y_calc:.2f}\ \mathrm{{km}}"
    )
    if "nu" in kepler:
        st.caption(f"Anomalía verdadera ν = {float(kepler['nu']):.6f} rad.")


def _neville_substitutions(points: list, tableau: tuple, target_t: float) -> str:
    lines = []
    n = len(points)
    for degree in range(1, n):
        lines.append(f"**Grado {degree}**")
        for i in range(degree, n):
            t_i = float(points[i][0])
            t_prev = float(points[i - degree][0])
            q_right = float(tableau[i][degree - 1])
            q_left = float(tableau[i - 1][degree - 1])
            q_ij = float(tableau[i][degree])
            lines.append(
                f"- $Q_{{{i},{degree}}}"
                f"=\\dfrac{{({target_t:.4f}-{t_prev:.4f})\\cdot ({q_right:.2f})"
                f"-({target_t:.4f}-{t_i:.4f})\\cdot ({q_left:.2f})}}"
                f"{{{t_i:.4f}-{t_prev:.4f}}}"
                f"={q_ij:.2f}$ km"
            )
    return "\n".join(lines)


def _neville_frame(points: list, tableau: tuple) -> pd.DataFrame:
    n = len(points)
    rows = []
    for i in range(n):
        row = {"Nodo": i, "t (rad)": f"{float(points[i][0]):.4f}"}
        for degree in range(n):
            header = f"Grado {degree}"
            if degree <= i:
                row[header] = f"{float(tableau[i][degree]):.2f}"
            else:
                row[header] = "—"
        rows.append(row)
    return pd.DataFrame(rows)


def _render_neville_process(interps: list[dict]) -> None:
    st.markdown("**Neville · último paquete perdido**")
    st.latex(
        r"Q_{i,j}=\frac{(t^*-t_{i-j})\,Q_{i,j-1}-(t^*-t_i)\,Q_{i-1,j-1}}{t_i-t_{i-j}}"
    )
    if not interps:
        st.caption("Neville entra al detectar un salto de identificador, con tres paquetes ya recibidos y el que cierra la brecha.")
        return

    item = interps[-1]
    nodes_x = item.get("nodes_x")
    tableau_x = item.get("tableau_x")
    if not nodes_x or not tableau_x:
        st.caption(
            f"Hueco del paquete {int(item['id']):02d} en t* = {item['timestamp']:.4f} rad · "
            f"X = {item['x']:.2f} km · Y = {item['y']:.2f} km. "
            "Lanza de nuevo la misión para ver cada celda de la tabla triangular."
        )
        return

    axis = st.segmented_control(
        "Coordenada",
        ["X", "Y"],
        default="X",
        key="eje_neville_vivo",
    )
    chosen = axis or "X"
    points = item["nodes_x"] if chosen == "X" else item["nodes_y"]
    tableau = item["tableau_x"] if chosen == "X" else item["tableau_y"]
    target_t = float(item["timestamp"])
    n = len(points)

    st.caption(
        f"Paquete {int(item['id']):02d} · t* = {target_t:.4f} rad · "
        f"nodos que sí llegaron, coordenada {chosen}."
    )
    st.markdown("**Grado 0** · los valores son la telemetría recibida.")
    degree_zero = ", ".join(
        f"$Q_{{{i},0}}={float(points[i][1]):.2f}$" for i in range(n)
    )
    st.markdown(degree_zero)
    st.markdown(_neville_substitutions(points, tableau, target_t))
    st.table(
        _neville_frame(points, tableau),
        hide_index=True,
        alt=f"Tabla triangular de Neville para la coordenada {chosen}",
    )
    if n > 1:
        estimate = float(tableau[n - 1][n - 1])
        previous = float(tableau[n - 1][n - 2])
        bound = abs(estimate - previous)
        st.latex(
            rf"\mathrm{{cota}}=|Q_{{{n - 1},{n - 1}}}-Q_{{{n - 1},{n - 2}}}|"
            rf"=|{estimate:.2f}-{previous:.2f}|={bound:.4f}\ \mathrm{{km}}"
        )


def _render_romberg_process(romberg_res, bounds: tuple[float, float] | None) -> None:
    st.markdown("**Romberg · distancia recorrida**")
    st.latex(
        r"S=\int_{t_0}^{t_f}v(t)\,dt,\qquad "
        r"R_{i,j}=\frac{4^j R_{i,j-1}-R_{i-1,j-1}}{4^j-1}"
    )
    if romberg_res is None or not romberg_res.tableau:
        st.caption("Romberg arranca al acumular tres muestras de velocidad a lo largo del arco.")
        return

    tableau = romberg_res.tableau
    if bounds is not None:
        t0, t1 = bounds
        st.latex(
            rf"R_{{0,0}}=\frac{{t_f-t_0}}{{2}}\left(v(t_0)+v(t_f)\right)"
            rf"={tableau[0][0]:.2f}\ \mathrm{{km}}"
            rf",\quad t_0={t0:.4f},\ t_f={t1:.4f}"
        )

    order = len(tableau) - 1
    rows = []
    for i, row_vals in enumerate(tableau):
        row = {"Nivel": i}
        for j in range(order + 1):
            header = "Trapecio" if j == 0 else f"Richardson {j}"
            if j < len(row_vals):
                row[header] = f"{row_vals[j]:.2f}"
            else:
                row[header] = "—"
        rows.append(row)
    st.table(pd.DataFrame(rows), hide_index=True, alt="Tabla de extrapolación de Romberg del arco orbital")

    if order >= 1:
        factor = 4 ** order
        left = float(tableau[order][order - 1])
        up = float(tableau[order - 1][order - 1])
        value = (factor * left - up) / (factor - 1.0)
        st.latex(
            rf"R_{{{order},{order}}}"
            rf"=\frac{{{factor}\cdot {left:.2f}-{up:.2f}}}{{{factor}-1}}"
            rf"={value:.2f}\ \mathrm{{km}}"
        )
    st.caption(
        f"{romberg_res.method} · {romberg_res.function_evaluations} evaluaciones de v(t) · "
        f"cota {romberg_res.estimated_error:.3e} km."
    )


def render_calculation_process(
    packets: list[dict],
    interps: list[dict],
    romberg_res,
    romberg_bounds: tuple[float, float] | None,
) -> None:
    st.subheader("Cálculo en vivo", divider=True, icon=":material/function:")
    st.caption("Cada refresco reescribe el paquete en antena, el último hueco de Neville y el arco de Romberg.")

    st.markdown("**Bitácora de la misión**")
    recent_packets = [_newton_log_line(pkt) for pkt in packets[-4:]]
    recent_packets.reverse()
    st.markdown(f"**{recent_packets[0]}**")
    if len(recent_packets) > 1:
        st.markdown("\n".join(f"- {line}" for line in recent_packets[1:]))
    if interps:
        recent_gaps = [_neville_log_line(item) for item in interps[-4:]]
        recent_gaps.reverse()
        st.markdown("**Huecos reconstruidos**")
        st.markdown("\n".join(f"- {line}" for line in recent_gaps))

    _render_newton_process(packets[-1])
    st.divider()
    _render_neville_process(interps)
    st.divider()
    _render_romberg_process(romberg_res, romberg_bounds)


# --------------------------------------------------------------------------- #
# Fragmento Reactivo en Tiempo Real
# --------------------------------------------------------------------------- #
def telemetry_dashboard(hub: OrbitalTelemetryHub, launcher: SatelliteLauncher) -> None:
    packets, interps, dropped, next_id = hub.snapshot()
    total_rx = len(packets)
    connected = total_rx > 0 and (time.time() - packets[-1]["rx"]) < 2.5

    sat_title = packets[-1].get("sat_name", "ISS (ZARYA)") if total_rx > 0 else "ISS (ZARYA) · CelesTrak"

    # Cabecera de estado de enlace
    c_status = OK if connected else (ALERT if total_rx > 0 else MUTED)
    t_status = "Enlace Activo (UDP :5005)" if connected else ("Enlace en Standby" if total_rx > 0 else "Sin Señal")

    st.markdown(f"""
    <div style="display:flex; justify-content:space-between; align-items:flex-end; margin-bottom:1.2rem;">
      <div>
        <h1 style="margin:0; font-size:1.85rem; font-weight:700;">🛰️ Gemelo Digital: {sat_title}</h1>
        <p style="margin:0; color:{MUTED}; font-size:0.9rem;">Parámetros CelesTrak en Vivo · Resolución Kepler-Newton · Neville · Romberg</p>
      </div>
      <div style="padding:0.4rem 0.85rem; border:1px solid {LINE}; border-radius:99px; background:{SURFACE}; font-size:0.85rem;">
        <span style="color:{c_status}">●</span> <strong>{sat_title}</strong> · {t_status} · {total_rx} pkts
      </div>
    </div>
    """, unsafe_allow_html=True)

    # Métricas y KPIs
    loss_pct = (dropped / (total_rx + dropped) * 100) if (total_rx + dropped) > 0 else 0
    arc_len = 0.0
    romberg_res = None
    romberg_bounds = None

    if total_rx >= 4:
        xs = [p["x"] for p in packets]
        ys = [p["y"] for p in packets]
        ts = [p["timestamp"] for p in packets]

        # Longitud de arco acumulada por diferencias
        for idx in range(1, len(ts)):
            arc_len += math.sqrt((xs[idx] - xs[idx - 1])**2 + (ys[idx] - ys[idx - 1])**2)

        # Cuadratura de Romberg sobre perfil de velocidad estimada
        speeds = []
        speed_ts = []
        for idx in range(1, len(ts)):
            dt = ts[idx] - ts[idx - 1]
            if dt > 0:
                dist = math.sqrt((xs[idx] - xs[idx - 1])**2 + (ys[idx] - ys[idx - 1])**2)
                speeds.append(dist / dt)
                speed_ts.append((ts[idx] + ts[idx - 1]) / 2.0)

        if len(speeds) >= 3:
            def _eval_speed(t_query: float) -> float:
                if t_query <= speed_ts[0]:
                    return speeds[0]
                if t_query >= speed_ts[-1]:
                    return speeds[-1]
                for k in range(len(speed_ts) - 1):
                    if speed_ts[k] <= t_query <= speed_ts[k + 1]:
                        frac = (t_query - speed_ts[k]) / (speed_ts[k + 1] - speed_ts[k])
                        return speeds[k] + frac * (speeds[k + 1] - speeds[k])
                return speeds[-1]

            try:
                romberg_res = romberg_integration(_eval_speed, speed_ts[0], speed_ts[-1], max_order=4)
                romberg_bounds = (speed_ts[0], speed_ts[-1])
            except Exception:
                romberg_res = None

    k1, k2, k3, k4 = st.columns(4)
    with k1:
        st.markdown(kpi("Telemetría Recibida", str(total_rx), " pkts", f"Próximo ID: {next_id}", ICE), unsafe_allow_html=True)
    with k2:
        st.markdown(kpi("Paquetes Perdidos", str(dropped), " pkts", f"{loss_pct:.1f}% tasa de pérdida", ALERT if dropped > 0 else OK), unsafe_allow_html=True)
    with k3:
        err_txt = f"±{interps[-1]['error']:.2f} km" if interps else "Sin pérdidas activas"
        st.markdown(kpi("Interpolaciones Neville", str(len(interps)), " pts", err_txt, AMBER if interps else MUTED), unsafe_allow_html=True)
    with k4:
        arc_txt = f"{arc_len:.0f}" if arc_len > 0 else "—"
        romberg_txt = f"Romberg: {romberg_res.value:.0f} km" if romberg_res else "Calculando…"
        st.markdown(kpi("Distancia Orbital", arc_txt, " km", romberg_txt, OK), unsafe_allow_html=True)

    if total_rx > 0:
        render_direct_ticker(packets, interps)

    st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)

    # Visualización gráfica
    col_orbit, col_ts = st.columns([1.3, 1])

    with col_orbit:
        st.markdown(f"<p style='color:{MUTED};font-size:.85rem;margin-bottom:.3rem'>Órbita 3D Inercial (Geocéntrica)</p>", unsafe_allow_html=True)
        fig_3d = go.Figure()

        # Cargar Tierra desde caché
        for trace in make_earth_wireframe_cached():
            fig_3d.add_trace(trace)

        if total_rx > 0:
            coords = [orbital_to_3d(p["x"], p["y"]) for p in packets]
            ox = [c[0] for c in coords]
            oy = [c[1] for c in coords]
            oz = [c[2] for c in coords]
            p_ids = [p["id"] for p in packets]

            # Trayectoria 3D
            fig_3d.add_trace(go.Scatter3d(
                x=ox, y=oy, z=oz, mode="lines+markers",
                name="Telemetría Real",
                line=dict(color=ICE, width=4),
                marker=dict(size=2.5, color=p_ids, colorscale=[[0, "#1d4ed8"], [0.5, ICE], [1, "#dbeafe"]], showscale=False),
                customdata=p_ids,
                hovertemplate="PKT %{customdata}<br>X=%{x:.0f} km<br>Y=%{y:.0f} km<br>Z=%{z:.0f} km<extra></extra>",
            ))

            # --- MODELO 3D REAL DEL SATÉLITE (CUERPO + PANELES SOLARES) ---
            last_pt = coords[-1]
            prev_pt = coords[-2] if len(coords) >= 2 else None
            
            # Inyectamos el satélite 3D orientado según el vector de movimiento
            sat_traces = make_satellite_3d_mesh(pos=last_pt, prev_pos=prev_pt, scale=700.0)
            for st_trace in sat_traces:
                fig_3d.add_trace(st_trace)

            # Punto de anclaje para el hover con información de telemetría
            fig_3d.add_trace(go.Scatter3d(
                x=[last_pt[0]], y=[last_pt[1]], z=[last_pt[2]], mode="markers",
                name="Telemetría Activa",
                marker=dict(size=4, color="white", opacity=0.8),
                hovertemplate="🛰️ Satélite Activo<br>X=%{x:.0f} km<br>Y=%{y:.0f} km<br>Z=%{z:.0f} km<extra></extra>",
                showlegend=False
            ))

        if interps:
            nev_coords = [orbital_to_3d(p["x"], p["y"]) for p in interps]
            fig_3d.add_trace(go.Scatter3d(
                x=[c[0] for c in nev_coords],
                y=[c[1] for c in nev_coords],
                z=[c[2] for c in nev_coords],
                mode="markers",
                name="Neville (Reconstrucción)",
                marker=dict(size=6, color=ALERT, symbol="x", line=dict(width=1, color="white")),
                customdata=[f"±{p['error']:.2f}" for p in interps],
                hovertemplate="Reconstrucción Neville<br>X=%{x:.0f} km<br>Y=%{y:.0f} km<br>Z=%{z:.0f} km<br>Error: %{customdata} km<extra></extra>",
            ))

        axis_limit = 11000
        axis_cfg = dict(
            range=[-axis_limit, axis_limit], backgroundcolor=BG,
            gridcolor=LINE, showbackground=True, zerolinecolor=LINE,
            title_font=dict(size=10, family="IBM Plex Mono"),
        )

        fig_3d.update_layout(
            template="plotly_dark",
            paper_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=0, r=0, t=0, b=0),
            height=510,
            scene=dict(
                xaxis=dict(title="X (km)", **axis_cfg),
                yaxis=dict(title="Y (km)", **axis_cfg),
                zaxis=dict(title="Z (km)", **axis_cfg),
                aspectmode="cube",
                camera=dict(
                    eye=dict(x=1.5, y=1.1, z=0.8),
                    up=dict(x=0, y=0, z=1),
                ),
            ),
            legend=dict(
                orientation="h", yanchor="bottom", y=0.98,
                font=dict(size=10, color=MUTED, family="Sora"),
                bgcolor="rgba(0,0,0,0)",
            ),
            font=dict(family="IBM Plex Mono", color=MUTED),
            uirevision="orbit_3d_state",  # Preserva rotación y zoom del usuario entre refrescos
        )
        st.plotly_chart(fig_3d, key="chart_3d_view", width="stretch")
        

    with col_ts:
        st.markdown(f"<p style='color:{MUTED};font-size:.85rem;margin-bottom:.3rem'>Descomposición de Coordenadas 2D</p>", unsafe_allow_html=True)
        fig_ts = make_subplots(
            rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
            subplot_titles=("Posición X en Plano Orbital (km)", "Posición Y en Plano Orbital (km)")
        )

        if total_rx > 0:
            ts_vals = [p["timestamp"] for p in packets]
            fig_ts.add_trace(go.Scatter(
                x=ts_vals, y=[p["x"] for p in packets], mode="lines",
                name="X Nominal", line=dict(color=ICE, width=1.5)
            ), row=1, col=1)
            fig_ts.add_trace(go.Scatter(
                x=ts_vals, y=[p["y"] for p in packets], mode="lines",
                name="Y Nominal", line=dict(color=OK, width=1.5)
            ), row=2, col=1)

        if interps:
            i_ts = [p["timestamp"] for p in interps]
            fig_ts.add_trace(go.Scatter(
                x=i_ts, y=[p["x"] for p in interps], mode="markers",
                name="X Neville", marker=dict(size=6, color=ALERT, symbol="x")
            ), row=1, col=1)
            fig_ts.add_trace(go.Scatter(
                x=i_ts, y=[p["y"] for p in interps], mode="markers",
                name="Y Neville", marker=dict(size=6, color=ALERT, symbol="x")
            ), row=2, col=1)

        fig_ts.update_layout(
            template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor=SURFACE, margin=dict(l=35, r=15, t=25, b=25),
            height=510, showlegend=False,
            font=dict(family="IBM Plex Mono", color=MUTED, size=11),
            uirevision="ts_plot_state",
        )
        fig_ts.update_xaxes(gridcolor=LINE, zerolinecolor=LINE)
        fig_ts.update_yaxes(gridcolor=LINE, zerolinecolor=LINE)
        st.plotly_chart(fig_ts, key="chart_ts_view", width="stretch")

    if total_rx > 0:
        render_calculation_process(packets, interps, romberg_res, romberg_bounds)

    # Tabla de auditoría numérica
    if interps:
        st.markdown(f"<p style='color:{MUTED};font-size:.85rem;margin-top:0.8rem'>📊 Auditoría de Interpolación Polinómica (Algoritmo de Neville)</p>", unsafe_allow_html=True)
        # Normalizar claves para soportar datos previos en caché
        normalized_interps = []
        for item in interps:
            row = dict(item)
            if "id" not in row:
                row["id"] = "—"
            normalized_interps.append(row)

        df_interp = pd.DataFrame(normalized_interps)
        rename_map = {
            "id": "ID Paquete",
            "timestamp": "Timestamp (rad)",
            "x": "X Estimada (km)",
            "y": "Y Estimada (km)",
            "error": "Cota Residuo (km)"
        }
        df_interp = df_interp.rename(columns=rename_map)

        # Ordenar columnas presentes
        desired_order = ["ID Paquete", "Timestamp (rad)", "X Estimada (km)", "Y Estimada (km)", "Cota Residuo (km)"]
        ordered_cols = [c for c in desired_order if c in df_interp.columns]
        df_interp = df_interp[ordered_cols]

        formatters = {}
        if "ID Paquete" in df_interp.columns:
            formatters["ID Paquete"] = lambda val: f"PKT {int(val):02d}" if isinstance(val, (int, float)) else str(val)
        if "Timestamp (rad)" in df_interp.columns:
            formatters["Timestamp (rad)"] = "{:.4f}"
        if "X Estimada (km)" in df_interp.columns:
            formatters["X Estimada (km)"] = "{:.2f}"
        if "Y Estimada (km)" in df_interp.columns:
            formatters["Y Estimada (km)"] = "{:.2f}"
        if "Cota Residuo (km)" in df_interp.columns:
            formatters["Cota Residuo (km)"] = "{:.4f}"

        st.dataframe(df_interp.style.format(formatters), width="stretch")
    elif total_rx == 0:
        st.markdown('<div class="empty">⏳ Esperando enlace de telemetría UDP…<br><small>Presiona el botón superior para lanzar la misión.</small></div>', unsafe_allow_html=True)


# --------------------------------------------------------------------------- #
# Punto de Entrada
# --------------------------------------------------------------------------- #
st.markdown(CSS, unsafe_allow_html=True)
hub = get_hub()

if hub.bind_error:
    st.error(hub.bind_error)
    st.info("Asegúrate de no tener otra consola ejecutando un listener en el puerto UDP 5005.")
    st.stop()

launcher = get_launcher()

# Barra de control y lanzamiento (Eje estático fuera del loop reactivo)
btn_col, status_col = st.columns([1.2, 4])
with btn_col:
    if st.button("🚀 Lanzar Satélite", width="stretch", type="primary"):
        started = launcher.launch(hub)
        if not started:
            st.toast("⚠️ Misión ya en curso. Espera a que complete la órbita.", icon="🛰️")
        else:
            st.toast("🚀 ¡Satélite lanzado! Transmitiendo paquetes UDP...", icon="✅")

with status_col:
    if launcher.is_running:
        st.markdown(f"<p style='color:{OK};font-size:.85rem;margin-top:.55rem'>🛰️ Satélite en vuelo — transmitiendo coordenadas a :{UDP_PORT}…</p>", unsafe_allow_html=True)
    else:
        st.markdown(f"<p style='color:{MUTED};font-size:.85rem;margin-top:.55rem'>💤 Satélite en plataforma — pulsa para iniciar simulación orbital.</p>", unsafe_allow_html=True)

# Invocación del fragmento aislado (refresco a 2 Hz sin reconstruir el DOM general)
st.fragment(run_every=0.5)(telemetry_dashboard)(hub, launcher)