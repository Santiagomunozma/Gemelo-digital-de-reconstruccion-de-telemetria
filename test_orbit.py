import math
from kepler_solver import compute_orbital_coordinates

def run_orbit_benchmark():
    # Parámetros: órbita elíptica LEO (semieje mayor 7000 km, e = 0.25)
    semi_axis_km = 7000.0
    eccentricity = 0.25
    steps = 8

    print(f"{'M (deg)':>7} | {'X (km)':>10} | {'Y (km)':>10} | {'True Anom (ν)':>13} | {'Iter':>4} | {'Residual':>10}")
    print("-" * 65)

    for i in range(steps + 1):
        deg = i * (360 // steps)
        m_rad = math.radians(deg)
        pos = compute_orbital_coordinates(m_rad, eccentricity, semi_axis_km)
        meta = pos.solver_metadata
        
        nu_deg = math.degrees(pos.true_anomaly) % 360
        print(f"{deg:7d}° | {pos.x:10.2f} | {pos.y:10.2f} | {nu_deg:12.2f}° | {meta.iterations:4d} | {meta.residual_error:10.2e}")

if __name__ == "__main__":
    run_orbit_benchmark()