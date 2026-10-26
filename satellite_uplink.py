import json
import math
import random
import socket
import sys
import time
import urllib.request
from kepler_solver import compute_orbital_coordinates

# Configuración de salida UTF-8 para evitar errores de codificación en consolas Windows
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


def safe_print(*args, **kwargs):
    """Imprime mensajes protegiendo contra errores charmap en Windows."""
    try:
        print(*args, **kwargs)
    except UnicodeEncodeError:
        try:
            cleaned = [str(a).encode("ascii", "replace").decode("ascii") for a in args]
            print(*cleaned, **kwargs)
        except Exception:
            pass


# Catálogo NORAD CatNR 25544 = Estación Espacial Internacional (ISS)
CELESTRAK_ISS_URL = "https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=json"


def fetch_iss_orbital_elements() -> tuple[float, float, str]:
    """Consulta la API pública de CelesTrak para extraer los parámetros TLE reales de la ISS."""
    safe_print("🌐 [RED] Consultando API pública de CelesTrak (NORAD Satellite Catalog)...")
    try:
        req = urllib.request.Request(
            CELESTRAK_ISS_URL,
            headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
        )
        with urllib.request.urlopen(req, timeout=3.5) as response:
            data = json.loads(response.read().decode("utf-8"))
            target = data[0] if isinstance(data, list) else data

            # Excentricidad real y nombre oficial
            ecc = float(target.get("ECCENTRICITY", 0.0008))
            sat_name = target.get("OBJECT_NAME", "ISS (ZARYA)")

            # Cálculo astrodinámico del semieje mayor 'a' a partir del Movimiento Medio (rev/día)
            # Tercera Ley de Kepler: n = sqrt(mu / a^3) -> a = (mu * (T / 2pi)^2)^(1/3)
            mean_motion = float(target.get("MEAN_MOTION", 15.49))
            period_seconds = 86400.0 / mean_motion
            mu_earth = 398600.4418  # Parámetro gravitacional estándar terrestre (km^3/s^2)
            semi_major_axis = (mu_earth * (period_seconds / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)

            safe_print(f"✅ [ÉXITO] Elementos orbitales cargados para {sat_name}:")
            safe_print(f"   • Excentricidad real (e): {ecc:.6f}")
            safe_print(f"   • Semieje mayor real (a): {semi_major_axis:.2f} km")
            return semi_major_axis, ecc, sat_name

    except Exception as err:
        safe_print(f"⚠️ [AVISO] Conexión externa no disponible ({err}). Activando parámetros de respaldo LEO.")
        # Fallback nominal con parámetros de órbita baja terrestre equivalentes
        return 6788.0, 0.0008, "ISS (Respaldo Local)"


def run_satellite():
    UDP_IP = "127.0.0.1"
    PORT_SAT = 5006
    PORT_STATION = 5005

    # 1. Extracción de parámetros reales en vivo
    semi_axis, ecc, sat_name = fetch_iss_orbital_elements()

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    except Exception:
        pass

    try:
        sock.bind((UDP_IP, PORT_SAT))
    except Exception as exc:
        safe_print(f"Aviso de enlace en :{PORT_SAT} ({exc}), transmitiendo...")

    sock.settimeout(0.5)

    safe_print(f"🛰️ Telemetría de [{sat_name}] inicializada en socket local :{PORT_SAT}")
    safe_print(f"📡 Esperando señal de sincronización (GROUND_READY) en :{PORT_STATION}...")

    # 2. Handshake con la estación terrena (con fallback por tiempo)
    connected = False
    start_handshake = time.time()
    while not connected and (time.time() - start_handshake < 1.5):
        try:
            msg, addr = sock.recvfrom(1024)
            if msg == b"GROUND_READY":
                safe_print(f"🚀 ¡Enlace establecido con Estación Terrena {addr}! Comenzando telemetría...")
                connected = True
                break
        except (socket.timeout, ConnectionResetError, OSError):
            pass

    if not connected:
        safe_print("🚀 Enlace directo iniciado hacia Estación Terrena...")

    # 3. Transmisión continua de órbita con inyección de fallas físicas
    steps = 90  # Muestreo angular uniforme en 2*pi
    blackout_range = range(30, 48)  # Zona de desvanecimiento ionosférico inducido

    try:
        for i in range(steps):
            t = i * (2.0 * math.pi / steps)
            # Resuelve Kepler con Newton-Raphson sobre los parámetros reales
            pos = compute_orbital_coordinates(t, ecc, semi_axis)

            payload = {
                "id": i,
                "timestamp": round(t, 4),
                "x": round(pos.x, 2),
                "y": round(pos.y, 2),
                "sat_name": sat_name
            }

            # Simulación de pérdida en el canal físico UDP (80% de tasa de caída en blackout)
            if i in blackout_range and random.random() < 0.80:
                safe_print(f"[!] PKT {i:02d} | Desvanecimiento atmosférico (Paquete descartado en canal UDP)")
                time.sleep(0.12)
                continue

            sock.sendto(json.dumps(payload).encode("utf-8"), (UDP_IP, PORT_STATION))
            safe_print(f"[>] [{sat_name}] PKT {i:02d} -> X={payload['x']:.1f} km, Y={payload['y']:.1f} km")
            time.sleep(0.12)

        safe_print(f"\n🛰️ Transmisión de la órbita de {sat_name} completada. Socket liberado exitosamente.\n")
    finally:
        try:
            sock.close()
        except Exception:
            pass


if __name__ == "__main__":
    run_satellite()