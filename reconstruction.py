from dataclasses import dataclass
from typing import List, Tuple


@dataclass(frozen=True)
class InterpolationResult:
    estimated_value: float
    estimated_error: float


def lagrange_interpolate(points: List[Tuple[float, float]], target_t: float) -> float:
    """Interpola un valor en target_t usando el polinomio de Lagrange clasico.
    
    points: Lista de tuplas (tiempo, valor_posicion)
    """
    n = len(points)
    total = 0.0

    for i in range(n):
        ti, yi = points[i]
        basis = 1.0
        for j in range(n):
            if i != j:
                tj, _ = points[j]
                basis *= (target_t - tj) / (ti - tj)
        total += yi * basis

    return total


def neville_interpolate(points: List[Tuple[float, float]], target_t: float) -> InterpolationResult:
    """Evalua la interpolacion en target_t mediante el algoritmo recursivo de Neville.
    
    Retorna el valor estimado y una cota del error de aproximacion local.
    """
    n = len(points)
    if n == 0:
        raise ValueError("Se requiere al menos un punto de apoyo.")

    # Matriz/tabla para el algoritmo de Neville: Q[i][j]
    q = [[0.0] * n for _ in range(n)]

    for i in range(n):
        q[i][0] = points[i][1]  # Columna 0 con los valores f(t_i)

    # Llenado de la tabla triangular
    for j in range(1, n):
        for i in range(j, n):
            t_i = points[i][0]
            t_prev = points[i - j][0]
            
            numerator = (target_t - t_prev) * q[i][j - 1] - (target_t - t_i) * q[i - 1][j - 1]
            denominator = t_i - t_prev
            
            q[i][j] = numerator / denominator

    best_estimate = q[n - 1][n - 1]
    # La diferencia con la aproximacion del grado anterior sirve como indicador del residuo
    estimated_error = abs(q[n - 1][n - 1] - q[n - 1][n - 2]) if n > 1 else 0.0

    return InterpolationResult(estimated_value=best_estimate, estimated_error=estimated_error)