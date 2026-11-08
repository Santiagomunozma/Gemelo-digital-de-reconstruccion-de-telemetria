import math
import random
from typing import List, Tuple
from kepler_solver import compute_orbital_coordinates
from reconstruction import neville_interpolate, lagrange_interpolate

def generate_telemetry(steps: int, semi_axis: float, ecc: float) -> List[Tuple[float, float, float]]:
    """Genera datos de telemetria base (t, x, y) a lo largo de la orbita."""
    telemetry = []
    for i in range(steps):
        # Asumimos que el tiempo t es proporcional a la anomalia media M para este test
        t = i * (2.0 * math.pi / steps)
        pos = compute_orbital_coordinates(t, ecc, semi_axis)
        telemetry.append((t, pos.x, pos.y))
    return telemetry

def run_blackout_simulation():
    # Parametros orbitales
    a = 7000.0  # Semieje mayor (km)
    e = 0.25    # Excentricidad
    total_samples = 24  # Muestras por orbita (telemetria cada 15 grados aprox)
    
    print("--- Simulacion de Apagon de Telemetria (Blackout) ---")
    
    # 1. Generar la telemetria real (Ground Truth)
    real_data = generate_telemetry(total_samples, a, e)
    
    # 2. Simular un corte de conexion. Elegimos un indice al azar para iniciar el corte.
    # Evitamos los bordes para tener datos de soporte (nodos) a ambos lados.
    blackout_start = random.randint(4, total_samples - 8)
    blackout_length = 5  # Perderemos 5 muestras consecutivas
    blackout_end = blackout_start + blackout_length
    
    print(f"[*] Enlace estable. Datos recibidos: indices 0 a {blackout_start-1}")
    print(f"[*] ¡ALERTA! Pérdida de paquetes detectada. Indices {blackout_start} a {blackout_end-1} perdidos.")
    print(f"[*] Enlace recuperado. Datos recibidos: indices {blackout_end} en adelante.\n")

    # 3. Preparar los datos de soporte (nodos) para la interpolacion
    # Tomamos 3 puntos antes del corte y 3 puntos despues
    support_nodes_x = []
    support_nodes_y = []
    
    # Nodos antes del corte
    for i in range(blackout_start - 3, blackout_start):
        support_nodes_x.append((real_data[i][0], real_data[i][1]))
        support_nodes_y.append((real_data[i][0], real_data[i][2]))
        
    # Nodos despues del corte
    for i in range(blackout_end, blackout_end + 3):
        support_nodes_x.append((real_data[i][0], real_data[i][1]))
        support_nodes_y.append((real_data[i][0], real_data[i][2]))

    # 4. Intentar reconstruir el punto medio del corte
    target_idx = blackout_start + (blackout_length // 2)
    target_t, real_x, real_y = real_data[target_idx]
    
    print(f"--- Reconstruyendo posicion para el tiempo t = {target_t:.3f} rad (Indice {target_idx}) ---")
    
    # Reconstruccion con Lagrange (solo X para demostracion)
    est_x_lagrange = lagrange_interpolate(support_nodes_x, target_t)
    
    # Reconstruccion con Neville (X e Y, que nos da cota de error)
    res_x_neville = neville_interpolate(support_nodes_x, target_t)
    res_y_neville = neville_interpolate(support_nodes_y, target_t)
    
    # 5. Analisis de resultados
    error_x_km = abs(real_x - res_x_neville.estimated_value)
    error_y_km = abs(real_y - res_y_neville.estimated_value)
    
    print(f"\nGround Truth (Kepler/Newton):")
    print(f"  X = {real_x:8.2f} km, Y = {real_y:8.2f} km")
    
    print(f"\nPrediccion Neville (Reconstruccion usando {len(support_nodes_x)} nodos):")
    print(f"  X = {res_x_neville.estimated_value:8.2f} km, Y = {res_y_neville.estimated_value:8.2f} km")
    
    print(f"\nAnalisis de Error Absoluto:")
    print(f"  Error en X: {error_x_km:.2f} km")
    print(f"  Error en Y: {error_y_km:.2f} km")
    
    # Validacion: Lagrange debe dar el mismo resultado teorico que el valor de Neville
    diff_methods = abs(est_x_lagrange - res_x_neville.estimated_value)
    print(f"\nComprobacion numerica: Diferencia entre polinomios de Lagrange y Neville: {diff_methods:.2e} (Debe ser cercano a 0)")

if __name__ == "__main__":
    # Fijamos la semilla temporalmente para que el corte ocurra en el mismo sitio
    # mientras validamos que funciona, luego puedes comentar esta linea.
    random.seed(42) 
    run_blackout_simulation()