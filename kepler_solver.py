from dataclasses import dataclass
import math
from typing import Optional


@dataclass(frozen=True)
class NumericalResult:
    eccentric_anomaly: float
    iterations: int
    residual_error: float
    converged: bool


@dataclass(frozen=True)
class OrbitalPosition:
    x: float
    y: float
    true_anomaly: float
    solver_metadata: NumericalResult


def solve_kepler_newton(
    mean_anomaly: float,
    eccentricity: float,
    tolerance: float = 1e-10,
    max_iter: int = 100,
) -> NumericalResult:
    """Resuelve la ecuacion de Kepler M = E - e*sin(E) mediante Newton-Raphson.

    Garantiza convergencia cuadratica en la busqueda de la anomalia excentrica E.
    """
    if not (0.0 <= eccentricity < 1.0):
        raise ValueError("La excentricidad debe estar en el intervalo [0, 1).")

    # Reduccion de M al intervalo [0, 2*pi]
    m_norm = mean_anomaly % (2.0 * math.pi)

    # Estimacion inicial optima segun el valor de excentricidad
    if eccentricity < 0.8:
        e_k = m_norm
    else:
        e_k = math.pi if m_norm > math.pi else m_norm + eccentricity * math.sin(m_norm)

    for step in range(1, max_iter + 1):
        f_val = e_k - eccentricity * math.sin(e_k) - m_norm
        f_prime = 1.0 - eccentricity * math.cos(e_k)

        # Evitar divisiones por cero en casos patologicos
        if abs(f_prime) < 1e-12:
            break

        delta = f_val / f_prime
        e_next = e_k - delta

        if abs(delta) < tolerance:
            final_residual = abs(e_next - eccentricity * math.sin(e_next) - m_norm)
            return NumericalResult(
                eccentric_anomaly=e_next,
                iterations=step,
                residual_error=final_residual,
                converged=True,
            )

        e_k = e_next

    # Fallback si agota iteraciones
    final_residual = abs(e_k - eccentricity * math.sin(e_k) - m_norm)
    return NumericalResult(
        eccentric_anomaly=e_k,
        iterations=max_iter,
        residual_error=final_residual,
        converged=False,
    )


def compute_orbital_coordinates(
    mean_anomaly: float,
    eccentricity: float,
    semi_major_axis: float,
    tolerance: float = 1e-10,
) -> OrbitalPosition:
    """Calcula las coordenadas (x, y) en km a partir de la resolucion numerica."""
    res = solve_kepler_newton(mean_anomaly, eccentricity, tolerance)
    e_val = res.eccentric_anomaly

    # Posicion en el plano orbital 2D (plano del periapsis)
    x = semi_major_axis * (math.cos(e_val) - eccentricity)
    y = semi_major_axis * math.sqrt(1.0 - eccentricity**2) * math.sin(e_val)

    # Calculo del angulo real observado desde el foco (Tierra)
    true_anomaly = 2.0 * math.atan2(
        math.sqrt(1.0 + eccentricity) * math.sin(e_val / 2.0),
        math.sqrt(1.0 - eccentricity) * math.cos(e_val / 2.0),
    )

    return OrbitalPosition(
        x=x, y=y, true_anomaly=true_anomaly, solver_metadata=res
    )


if __name__ == "__main__":
    # Test de verificacion: Satelite en orbita eliptica tipo Molniya o LEO de prueba
    semi_axis_km = 7000.0  # ~600 km de altitud media sobre la Tierra
    ecc = 0.25  # Orbita visiblemente eliptica
    m_test = math.radians(45.0)

    pos = compute_orbital_coordinates(m_test, ecc, semi_axis_km)
    meta = pos.solver_metadata

    print(f"Estado de convergencia: {meta.converged}")
    print(f"Iteraciones requeridas: {meta.iterations}")
    print(f"Error residual |f(E)|: {meta.residual_error:.2e}")
    print(f"Coordenadas del satelite: X={pos.x:.3f} km, Y={pos.y:.3f} km")