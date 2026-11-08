import math
from orbital_metrics import composite_simpson, romberg_integration

def run_integration_benchmark():
    # Modelo de potencia solar recibida (W) en funcion de la anomalia media M in [0, 2pi]
    # Supongamos que el satelite pasa la mitad de la orbita al sol y la otra mitad en eclipse
    p_peak = 1200.0  # 1200 W de potencia pico de los paneles
    
    def solar_power(m: float) -> float:
        # Al sol entre 0 y pi, eclipse entre pi y 2pi
        illumination = math.cos(m - math.pi / 2.0)
        return p_peak * max(0.0, illumination)

    a_time, b_time = 0.0, 2.0 * math.pi

    print("--- Benchmark de Integración Numérica: Energía Solar Acumulada ---")
    
    # 1. Simpson con 16 subintervalos
    res_simp = composite_simpson(solar_power, a_time, b_time, n=16)
    
    # 2. Romberg con cota adaptativa
    res_romberg = romberg_integration(solar_power, a_time, b_time, max_order=6)

    # Valor analitico exacto de la integral de p_peak * cos(t) entre 0 y pi = 2 * p_peak = 2400.0
    exact_energy = 2400.0

    print(f"Valor analítico exacto: {exact_energy:.4f} Joules normalizados\n")
    print(f"[{res_simp.method}]")
    print(f"  Resultado:    {res_simp.value:.4f}")
    print(f"  Evaluaciones: {res_simp.function_evaluations}")
    print(f"  Error real:   {abs(res_simp.value - exact_energy):.2e}")

    print(f"\n[{res_romberg.method}]")
    print(f"  Resultado:    {res_romberg.value:.4f}")
    print(f"  Evaluaciones: {res_romberg.function_evaluations}")
    print(f"  Error real:   {abs(res_romberg.value - exact_energy):.2e}")

if __name__ == "__main__":
    run_integration_benchmark()