# 🛰️ Guía Maestra de Exposición · Gemelo Digital y Métodos Numéricos

> **Proyecto:** Gemelo Digital Orbital en Tiempo Real con Reconstrucción de Telemetría UDP  
> **Asignatura:** Métodos Numéricos / Análisis Numérico  
> **Caso de Estudio:** Estación Espacial Internacional (ISS - ZARYA, Catálogo NORAD #25544)

---

## 📋 Índice General de la Presentación

1. [Introducción y Objetivos del Proyecto](#1-introducción-y-objetivos-del-proyecto)
2. [Arquitectura del Gemelo Digital](#2-arquitectura-del-gemelo-digital)
3. [Método 1: Solución de la Ecuación de Kepler (Newton-Raphson)](#3-método-1-solución-de-la-ecuación-de-kepler-newton-raphson)
4. [Método 2: Reconstrucción Numérica con Algoritmo de Neville (NÚCLEO DE LA EXPOSICIÓN)](#4-método-2-reconstrucción-numérica-con-algoritmo-de-neville)
   - [¿Qué es y por qué lo usamos en lugar de Lagrange tradicional?](#por-qué-neville-y-no-lagrange-directo)
   - [Significado de cada variable y número](#significado-de-cada-variable-y-número)
   - [Construcción de la Tabla Triangular paso a paso](#construcción-de-la-tabla-triangular)
   - [Ejemplo Numérico Real con Coordenadas de la ISS](#ejemplo-numérico-resuelto-paso-a-paso)
   - [¿Qué significa la Cota de Residuo/Error?](#qué-significa-la-cota-de-residuo-o-error)
5. [Método 3: Cuadratura Numérica y Longitud de Arco (Romberg)](#5-método-3-cuadratura-numérica-romberg)
6. [Preguntas Típicas del Profesor y Respuestas Clave](#6-preguntas-típicas-del-profesor-defensa-del-proyecto)

---

## 1. Introducción y Objetivos del Proyecto

### ¿Qué es este proyecto?
Es un **Gemelo Digital (Digital Twin)** para el seguimiento orbital en tiempo real de la Estación Espacial Internacional (ISS). 

### ¿Qué problema real de ingeniería resuelve?
Cuando un satélite orbita la Tierra a 27,600 km/h, atraviesa zonas de sombra de antenas o perturbaciones en la ionosfera (como tormentas solares). En esos momentos, la señal de telemetría UDP **se degrada y se pierden paquetes consecutivos**.

La estación terrena no puede esperar a que el satélite vuelva a tener cobertura para saber dónde está: **debe reconstruir matemáticamente en tiempo real la posición física donde se encontraba el satélite cuando el paquete no llegó**, calculando además **con qué margen de incertidumbre (error)** se hace dicha estimación.

---

## 2. Arquitectura del Gemelo Digital

El software está dividido en tres capas concurrentes:

```
┌─────────────────────────────────────────────────────────────┐
│ 1. API CelesTrak (NORAD GP)                                 │
│    Descarga en vivo: semieje mayor (a) y excentricidad (e)   │
└──────────────────────────────┬──────────────────────────────┘
                               │ Parámetros orbitales reales
┌──────────────────────────────▼──────────────────────────────┐
│ 2. Emisor Satelital (satellite_uplink.py)                   │
│    - Resuelve Kepler vía Newton-Raphson                      │
│    - Transmite coordenadas X, Y por UDP en puerto 5005       │
│    - Simula corte de señal (Blackout / Tormenta solar)       │
└──────────────────────────────┬──────────────────────────────┘
                               │ Paquetes UDP (con pérdidas inducidas)
┌──────────────────────────────▼──────────────────────────────┐
│ 3. Estación Terrena y Dashboard (dashboard.py)              │
│    - Detección de saltos en ID de paquetes                   │
│    - Reconstrucción polinómica recursiva (Neville)          │
│    - Integración de longitud de arco (Romberg)               │
│    - Renderizado 3D inercial con malla geométrica            │
└─────────────────────────────────────────────────────────────┘
```

---

## 3. Método 1: Solución de la Ecuación de Kepler (Newton-Raphson)

Para saber dónde está el satélite en un instante angular $t \in [0, 2\pi]$, debemos resolver la **Ecuación Trascendente de Kepler**:

$$M = E - e \cdot \sin(E)$$

* **$M$ (Anomalía Media):** Representa el tiempo transcurrido escalado a ángulo uniforme.
* **$e$ (Excentricidad):** Qué tan achatada es la elipse ($e \approx 0.00069$ para la ISS, casi circular).
* **$E$ (Anomalía Excéntrica):** **La incógnita geométrica**. No se puede despejar algebraicamente porque $E$ está a la vez fuera y dentro del seno.

### Aplicación de Newton-Raphson
Definimos la función de raíz:
$$f(E) = E - e \cdot \sin(E) - M = 0$$
Su derivada respecto a $E$:
$$f'(E) = 1 - e \cdot \cos(E)$$

Fórmula iterativa:
$$E_{k+1} = E_k - \frac{E_k - e \cdot \sin(E_k) - M}{1 - e \cdot \cos(E_k)}$$

* **Criterio de parada:** Tolerancia $|E_{k+1} - E_k| < 10^{-10}$ rad.
* **Convergencia:** Cuadrática (alcanza precisión de máquina en solo 3 a 4 iteraciones).
* Una vez obtenido $E$, las coordenadas orbitales 2D se calculan como:
  $$x = a(\cos E - e), \quad y = a\sqrt{1 - e^2}\sin E$$

---

## 4. Método 2: Reconstrucción Numérica con Algoritmo de Neville

> **ESTE ES EL NÚCLEO MATEMÁTICO DEL PROYECTO.**

### ¿Por qué Neville y no Lagrange directo?

1. **Forma de Lagrange tradicional:**
   $$P(t) = \sum_{i=0}^n y_i \prod_{j \ne i} \frac{t - t_j}{t_i - t_j}$$
   * Desventaja: Si cambias el número de nodos o quieres evaluar un único punto, debes recalcular todas las bases multiplicativas desde cero. Es computacionalmente costoso ($O(n^2)$ multiplicaciones por evaluación) y propenso a errores de redondeo.
2. **Algoritmo de Neville:**
   * Es una **formulación recursiva** que genera una tabla triangular.
   * Permite evaluar el polinomio interpolante directamente en el punto deseado $t^*$ sin expandir coeficientes simbólicos.
   * **Ventaja insustituible:** Compara automáticamente aproximaciones de orden $n$ contra orden $n-1$, entregando **una cota intrínseca del error residual local** sin necesidad de conocer la derivada $(n+1)$-ésima de la función.

---

### Significado de cada variable y número

Para interpolar la posición en el tiempo $t^*$ de un paquete perdido, tomamos $n=4$ nodos conocidos: 3 anteriores a la falla y 1 posterior a la falla:

$$\{(t_0, y_0), (t_1, y_1), (t_2, y_2), (t_3, y_3)\}$$

* **$t_i$ (Nodos Temporales):** Instante en radianes de los paquetes que **sí llegaron** a la estación.
* **$y_i$ (Valores Nodales):** Coordenada observada ($X$ o $Y$ en km) en el tiempo $t_i$.
* **$t^*$ (Timestamp Objetivo):** El instante exacto en el que debió haber llegado el paquete perdido.
* **$Q_{i,j}$ (Elemento de la Tabla de Neville):**
  * El subíndice $i$ indica el **nodo base inferior**.
  * El subíndice $j$ indica el **grado del polinomio interpolante**.
  * Por ejemplo:
    * $Q_{0,0} = y_0$: Polinomio de grado 0 (constante en el nodo 0).
    * $Q_{1,1}$: Recta (polinomio grado 1) que interpola entre $t_0$ y $t_1$, evaluada en $t^*$.
    * $Q_{2,2}$: Parábola (polinomio grado 2) que interpola entre $t_0, t_1, t_2$, evaluada en $t^*$.
    * $Q_{3,3}$: Polinomio cúbico final (grado 3) que interpola los 4 nodos, evaluado en $t^*$.

---

### La Relación de Recurrencia de Neville

Para cualquier elemento de la tabla:

$$Q_{i,j} = \frac{(t^* - t_{i-j}) \cdot Q_{i,j-1} - (t^* - t_i) \cdot Q_{i-1,j-1}}{t_i - t_{i-j}}$$

#### ¿De dónde sale esta fórmula?
Es una combinación convexa ponderada entre dos polinomios de grado $(j-1)$: uno que cubre los primeros nodos ($Q_{i-1,j-1}$) y otro que cubre los últimos ($Q_{i,j-1}$).

---

### Construcción de la Tabla Triangular

Visualmente, la tabla de Neville se llena de izquierda a derecha:

| $i$ | Nodo $t_i$ | Grado 0 ($j=0$) | Grado 1 ($j=1$) | Grado 2 ($j=2$) | Grado 3 ($j=3$) [Resultado] |
|:---:|:---:|:---:|:---:|:---:|:---:|
| **0** | $t_0$ | $Q_{0,0} = y_0$ | | | |
| **1** | $t_1$ | $Q_{1,0} = y_1$ | $Q_{1,1} = \frac{(t^*-t_0)Q_{1,0} - (t^*-t_1)Q_{0,0}}{t_1 - t_0}$ | | |
| **2** | $t_2$ | $Q_{2,0} = y_2$ | $Q_{2,1} = \frac{(t^*-t_1)Q_{2,0} - (t^*-t_2)Q_{1,0}}{t_2 - t_1}$ | $Q_{2,2} = \frac{(t^*-t_0)Q_{2,1} - (t^*-t_2)Q_{1,1}}{t_2 - t_0}$ | |
| **3** | $t_3$ | $Q_{3,0} = y_3$ | $Q_{3,1} = \frac{(t^*-t_2)Q_{3,0} - (t^*-t_3)Q_{2,0}}{t_3 - t_2}$ | $Q_{3,2} = \frac{(t^*-t_1)Q_{3,1} - (t^*-t_3)Q_{2,1}}{t_3 - t_1}$ | **$Q_{3,3}$ (Valor Final Estimado)** |

---

### Ejemplo Numérico Resuelto Paso a Paso

Supongamos que durante la órbita de la ISS se pierde el paquete intermedio. Queremos reconstruir la coordenada **$X$** en el instante **$t^* = 2.50$ rad**:

#### Datos de entrada (4 paquetes conocidos):
* Nodo 0: $t_0 = 2.20$, $X_0 = -4000.0\text{ km}$
* Nodo 1: $t_1 = 2.40$, $X_1 = -3500.0\text{ km}$
* Nodo 2: $t_2 = 2.60$, $X_2 = -2900.0\text{ km}$
* Nodo 3: $t_3 = 2.80$, $X_3 = -2200.0\text{ km}$

Punto a interpolar: $t^* = 2.50\text{ rad}$.

#### Paso 1: Columna de Grado 0 ($j=0$)
Son simplemente las coordenadas originales:
* $Q_{0,0} = -4000.0$
* $Q_{1,0} = -3500.0$
* $Q_{2,0} = -2900.0$
* $Q_{3,0} = -2200.0$

#### Paso 2: Columna de Grado 1 ($j=1$, Rectas secantes)
* **Cálculo de $Q_{1,1}$ (entre nodos 0 y 1):**
  $$Q_{1,1} = \frac{(2.50 - 2.20)(-3500.0) - (2.50 - 2.40)(-4000.0)}{2.40 - 2.20} = \frac{(0.30)(-3500.0) - (0.10)(-4000.0)}{0.20} = \frac{-1050 + 400}{0.20} = -3250.0\text{ km}$$
* **Cálculo de $Q_{2,1}$ (entre nodos 1 y 2):**
  $$Q_{2,1} = \frac{(2.50 - 2.40)(-2900.0) - (2.50 - 2.60)(-3500.0)}{2.60 - 2.40} = \frac{(0.10)(-2900.0) - (-0.10)(-3500.0)}{0.20} = \frac{-290 - 350}{0.20} = -3200.0\text{ km}$$
* **Cálculo de $Q_{3,1}$ (entre nodos 2 y 3):**
  $$Q_{3,1} = \frac{(2.50 - 2.60)(-2200.0) - (2.50 - 2.80)(-2900.0)}{2.80 - 2.60} = \frac{(-0.10)(-2200.0) - (-0.30)(-2900.0)}{0.20} = \frac{220 - 870}{0.20} = -3250.0\text{ km}$$

#### Paso 3: Columna de Grado 2 ($j=2$, Parábolas)
* **Cálculo de $Q_{2,2}$ (usando $Q_{2,1}$ y $Q_{1,1}$):**
  $$Q_{2,2} = \frac{(2.50 - 2.20)(-3200.0) - (2.50 - 2.60)(-3250.0)}{2.60 - 2.20} = \frac{(0.30)(-3200.0) - (-0.10)(-3250.0)}{0.40} = \frac{-960 - 325}{0.40} = -3212.5\text{ km}$$
* **Cálculo de $Q_{3,2}$ (usando $Q_{3,1}$ y $Q_{2,1}$):**
  $$Q_{3,2} = \frac{(2.50 - 2.40)(-3250.0) - (2.50 - 2.80)(-3200.0)}{2.80 - 2.40} = \frac{(0.10)(-3250.0) - (-0.30)(-3200.0)}{0.40} = \frac{-325 - 960}{0.40} = -3212.5\text{ km}$$

#### Paso 4: Columna de Grado 3 ($j=3$, Polinomio Cúbico Final)
* **Cálculo de $Q_{3,3}$:**
  $$Q_{3,3} = \frac{(2.50 - 2.20)(-3212.5) - (2.50 - 2.80)(-3212.5)}{2.80 - 2.20} = \frac{(0.30)(-3212.5) + (0.30)(-3212.5)}{0.60} = -3212.5\text{ km}$$

---

### ¿Qué significa la Cota de Residuo o Error?

En la columna final de la tabla del dashboard se lee: **Cota Residuo (km)**.

$$\text{Error Estimado} = |Q_{3,3} - Q_{3,2}|$$

#### Explicación intuitiva para el profesor:
* $Q_{3,2}$ es la mejor aproximación usando un polinomio parabólico (grado 2).
* $Q_{3,3}$ es la aproximación refinada agregando un nodo más (polinomio cúbico grado 3).
* La diferencia entre ambos mide **cuánto varió la estimación al aumentar un grado de libertad**. Si la diferencia es de apenas $\pm 0.05\text{ km}$, significa que el polinomio ya convergía y la predicción es altamente confiable. Si la diferencia es grande (por ejemplo $\pm 150\text{ km}$), indica que la brecha de paquetes perdidos fue muy ancha y el error por extrapolación aumenta.

---

## 5. Método 3: Cuadratura Numérica (Romberg)

### ¿Para qué sirve en este sistema?
El satélite no viaja en línea recta; recorre una elipse alabeada en 3 dimensiones. Para calcular la **Distancia Orbital Acumulada** recorrida a lo largo del tiempo, se integra el perfil de velocidad instantánea $v(t)$:

$$S = \int_{t_0}^{t_f} v(t) \, dt = \int_{t_0}^{t_f} \sqrt{\left(\frac{dx}{dt}\right)^2 + \left(\frac{dy}{dt}\right)^2} \, dt$$

### ¿Cómo funciona la Cuadratura de Romberg?
Combina la **Regla del Trapecio compuesta** con la **Extrapolación de Richardson** para acelerar la convergencia:

1. **Nivel 1 (Trapecio):** Divide el intervalo en subintervalos de paso $h_k = \frac{b - a}{2^k}$.
   $$R_{k, 0} = \frac{1}{2} R_{k-1, 0} + h_k \sum_{i=1}^{2^{k-1}} f(a + (2i - 1)h_k)$$
2. **Nivel 2 (Extrapolación de Richardson):** Cancela los términos dominantes del error de orden $O(h^2), O(h^4), O(h^6)$:
   $$R_{j, k} = \frac{4^k \cdot R_{j, k-1} - R_{j-1, k-1}}{4^k - 1}$$

* **Ventaja:** Con solo 4 o 5 particiones se alcanza la misma precisión que Trapecio con miles de subintervalos, calculando la distancia orbital exacta en milisegundos.

---

## 6. Preguntas Típicas del Profesor (Defensa del Proyecto)

### Pregunta 1: *"¿Por qué usan exactamente 4 nodos para Neville y no usan 15 o todos los paquetes recibidos?"*
> **Respuesta:**  
> Por el **Fenómeno de Runge**. Cuando se incrementa excesivamente el grado de un polinomio de interpolación sobre nodos equiespaciados, los bordes del polinomio experimentan oscilaciones descontroladas de gran amplitud. Usar un grado cúbico ($n=4$ nodos) acota el error, evita oscilaciones espurias y garantiza estabilidad numérica en tiempo real con complejidad $O(n^2)$.

---

### Pregunta 2: *"¿Por qué la tabla triangular usa nodos antes y después de la brecha y no solo los anteriores?"*
> **Respuesta:**  
> Porque trabajar con nodos a ambos lados de la brecha convierte el problema en una **interpolación interior**, donde el error está formalmente acotado por el término del resto de Cauchy:
> $$R_n(t) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod_{i=0}^n (t - t_i)$$
> Si solo usáramos los nodos anteriores, estaríamos haciendo **extrapolación exterior**, donde el término $\prod (t - t_i)$ crece explosivamente y los errores se disparan.

---

### Pregunta 3: *"¿Qué diferencia hay entre la Anomalía Media $M$, la Excéntrica $E$ y la Verdadera $\nu$?"*
> **Respuesta:**  
> * **$M$ (Media):** Es un ángulo ficticio que avanza a velocidad angular constante como si la órbita fuera un círculo perfecto ($M = n \cdot t$).
> * **$E$ (Excéntrica):** Es el parámetro geométrico auxiliar proyectado sobre la circunferencia circunscrita que permite linealizar el área barrida. Se obtiene con Newton-Raphson.
> * **$\nu$ (Verdadera):** Es el ángulo físico real medido desde el foco (el centro de la Tierra) hacia el satélite.

---

### Pregunta 4: *"¿Cómo garantizan que Newton-Raphson siempre converja al resolver Kepler?"*
> **Respuesta:**  
> Porque la derivada de la función de Kepler es $f'(E) = 1 - e \cdot \cos(E)$. Como la excentricidad de satélites terrestres cerrados cumple $0 \le e < 1$, se cumple estrictamente que:
> $$f'(E) \ge 1 - e > 0 \quad \forall E \in \mathbb{R}$$
> La primera derivada jamás se anula y la función es estrictamente creciente. Al no haber puntos con pendiente cero, Newton-Raphson tiene convergencia garantizada sin caer en divisiones por cero ni ciclos infinitos.

---

### Pregunta 5: *"¿De dónde salieron los números reales que se ven en el dashboard?"*
> **Respuesta:**  
> El emisor consulta en vivo la API de **CelesTrak (NORAD 25544)** que provee los elementos orbitales de la ISS actualizados por el Comando de Defensa Aeroespacial de EE.UU. Mediante la **Tercera Ley de Kepler generalizada**:
> $$a = \sqrt[3]{\mu \cdot \left(\frac{T}{2\pi}\right)^2}$$
> con $\mu = 398600.4418\text{ km}^3/\text{s}^2$, se obtiene el semieje mayor real de $6797.64\text{ km}$, ubicando a la estación a $\approx 426\text{ km}$ sobre la Tierra con excentricidad casi circular $e = 0.00069$.

---

## 7. Resumen Ejecutivo para la Diapositiva Final

| Método Numérico | Aplicación en el Proyecto | Ventaja Clave |
|---|---|---|
| **Newton-Raphson** | Solución de la Ec. de Kepler $M = E - e\sin(E)$ | Convergencia cuadrática en $\le 4$ iteraciones |
| **Algoritmo de Neville** | Reconstrucción de telemetría perdida en ráfaga | Evaluación recursiva eficiente + cota de error instantánea |
| **Cuadratura de Romberg** | Longitud de arco / distancia orbital acumulada | Aceleración de Richardson sobre la regla del trapecio |
