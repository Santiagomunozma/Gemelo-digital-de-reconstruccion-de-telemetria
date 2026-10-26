from dataclasses import dataclass
import math
from typing import Callable, List


@dataclass(frozen=True)
class IntegrationResult:
    method: str
    value: float
    function_evaluations: int
    estimated_error: float


def composite_trapezoidal(f: Callable[[float], float], a: float, b: float, n: int) -> float:
    """Aproximacion por regla del Trapecio Compuesta con n subintervalos."""
    if n <= 0:
        raise ValueError("El numero de subintervalos debe ser mayor a 0.")
    h = (b - a) / n
    total = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        total += f(a + i * h)
    return total * h


def composite_simpson(f: Callable[[float], float], a: float, b: float, n: int) -> IntegrationResult:
    """Regla de Simpson 1/3 Compuesta. n debe ser un numero par."""
    if n % 2 != 0:
        n += 1  # Forzar n par para cumplir la particion en pares de intervalos

    h = (b - a) / n
    total = f(a) + f(b)

    # Coeficientes 4 para indices impares y 2 para indices pares
    for i in range(1, n, 2):
        total += 4.0 * f(a + i * h)
    for i in range(2, n - 1, 2):
        total += 2.0 * f(a + i * h)

    res = (h / 3.0) * total
    return IntegrationResult(
        method="Simpson 1/3 Compuesto",
        value=res,
        function_evaluations=n + 1,
        estimated_error=0.0  # Simpson simple no estima error local intrinseco
    )

def romberg_integration(
    f: Callable[[float], float], 
    a: float, 
    b: float, 
    max_order: int = 6,  # Aumente max_order a 6 para asegurar profundidad
    tolerance: float = 1e-9
) -> IntegrationResult:
    """Metodo de Romberg mediante extrapolacion de Richardson sobre el Trapecio."""
    r = [[0.0] * (max_order + 1) for _ in range(max_order + 1)]
    eval_count = 0

    h = b - a
    r[0][0] = 0.5 * h * (f(a) + f(b))
    eval_count += 2

    for i in range(1, max_order + 1):
        n_sub = 2**i
        h = (b - a) / n_sub
        
        sum_new = sum(f(a + k * h) for k in range(1, n_sub, 2))
        eval_count += 2**(i - 1)
        r[i][0] = 0.5 * r[i - 1][0] + h * sum_new

        for j in range(1, i + 1):
            factor = 4.0**j
            r[i][j] = (factor * r[i][j - 1] - r[i - 1][j - 1]) / (factor - 1.0)

        error_est = abs(r[i][i] - r[i - 1][i - 1])
        # FIX Ciberseguro(?): Obligamos a profundizar al menos a orden 3 (8 subintervalos)
        # para evitar convergencia prematura si los bordes de f(x) son 0.
        if i >= 3 and error_est < tolerance:
            return IntegrationResult(
                method="Romberg (Richardson)",
                value=r[i][i],
                function_evaluations=eval_count,
                estimated_error=error_est
            )

    return IntegrationResult(
        method="Romberg (Richardson)",
        value=r[max_order][max_order],
        function_evaluations=eval_count,
        estimated_error=abs(r[max_order][max_order] - r[max_order - 1][max_order - 1])
    )