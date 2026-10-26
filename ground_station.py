import socket
import json
import time
import sys
from reconstruction import neville_interpolate

# Configuración de salida UTF-8 para consolas Windows (evita UnicodeEncodeError con emojis)
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
    """Imprime mensajes protegiendo contra errores de codificación en terminales Windows."""
    try:
        print(*args, **kwargs)
    except UnicodeEncodeError:
        try:
            cleaned = [str(a).encode("ascii", "replace").decode("ascii") for a in args]
            print(*cleaned, **kwargs)
        except Exception:
            pass


def run_ground_station():
    UDP_IP = "127.0.0.1"
    PORT_STATION = 5005
    PORT_SAT = 5006

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    except Exception:
        pass
    sock.bind((UDP_IP, PORT_STATION))
    sock.settimeout(6.0)

    safe_print(f"📡 Estación Terrena activa en :{PORT_STATION}. Sincronizando con satélite...")
    # Mandamos señal de inicio al satelite
    try:
        sock.sendto(b"GROUND_READY", (UDP_IP, PORT_SAT))
    except Exception:
        pass

    telemetry_history = []
    expected_id = 0

    try:
        while True:
            try:
                data, addr = sock.recvfrom(2048)
            except ConnectionResetError:
                time.sleep(0.5)
                try:
                    sock.sendto(b"GROUND_READY", (UDP_IP, PORT_SAT))
                except Exception:
                    pass
                continue
            packet = json.loads(data.decode("utf-8"))
            packet_id = packet["id"]

            if packet_id > expected_id:
                missing_count = packet_id - expected_id
                safe_print(f"\n⚠️ [DROP] Pérdida de {missing_count} paquetes (IDs {expected_id} a {packet_id - 1})")

                if len(telemetry_history) >= 3:
                    safe_print(f"⚙️ Ejecutando Neville para reconstruir los {missing_count} paquetes perdidos...")
                    nodes_x = [(p["timestamp"], p["x"]) for p in telemetry_history[-3:]] + [(packet["timestamp"], packet["x"])]
                    nodes_y = [(p["timestamp"], p["y"]) for p in telemetry_history[-3:]] + [(packet["timestamp"], packet["y"])]

                    prev_pkt = telemetry_history[-1]
                    id_prev = expected_id - 1
                    t_prev = prev_pkt["timestamp"]

                    for lost_id in range(expected_id, packet_id):
                        frac = (lost_id - id_prev) / (packet_id - id_prev)
                        t_est = t_prev + frac * (packet["timestamp"] - t_prev)
                        res_x = neville_interpolate(nodes_x, t_est)
                        res_y = neville_interpolate(nodes_y, t_est)
                        err_cota = max(res_x.estimated_error, res_y.estimated_error)

                        safe_print(f"   🎯 PKT {lost_id:02d} estimado en t={t_est:.4f}: X={res_x.estimated_value:.1f} km, Y={res_y.estimated_value:.1f} km (±{err_cota:.2f} km)")

            safe_print(f"[<] PKT {packet_id:02d} | t={packet['timestamp']:.2f}s | Pos=({packet['x']:.1f}, {packet['y']:.1f})")
            telemetry_history.append(packet)
            expected_id = packet_id + 1

    except socket.timeout:
        safe_print("\n📡 Fin de ventana de visibilidad orbital. Estación en standby.")
    finally:
        try:
            sock.close()
        except Exception:
            pass


if __name__ == "__main__":
    run_ground_station()