# 🛰️ Gemelo Digital Orbital · Métodos Numéricos

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.37%2B-FF4B4B.svg)](https://streamlit.io/)
[![Plotly](https://img.shields.io/badge/Plotly-5.20%2B-3F4F75.svg)](https://plotly.com/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Simulador y Gemelo Digital en Tiempo Real** para seguimiento satelital, reconstrucción de telemetría bajo pérdida de paquetes y cálculo de métricas orbitales mediante métodos numéricos avanzados.

---

## 📌 Tabla de Contenidos

- [Descripción General](#-descripción-general)
- [Fundamentación Matemática y Métodos Numéricos](#-fundamentación-matemática-y-métodos-numéricos)
  - [1. Ecuación de Kepler y Newton-Raphson](#1-ecuación-de-kepler-y-newton-raphson)
  - [2. Reconstrucción de Telemetría (Algoritmo de Neville)](#2-reconstrucción-de-telemetría-algoritmo-de-neville)
  - [3. Cuadratura Numérica (Romberg y Simpson 1/3)](#3-cuadratura-numérica-romberg-y-simpson-13)
- [Arquitectura del Sistema](#-arquitectura-del-sistema)
- [Estructura del Proyecto](#-estructura-del-proyecto)
- [Instalación y Requisitos](#-instalación-y-requisitos)
- [Guía de Uso](#-guía-de-uso)
  - [Opción 1: Dashboard Integrado (Recomendada)](#opción-1-dashboard-integrado-recomendada)
  - [Opción 2: Modo Consola (Emisor / Receptor UDP Separados)](#opción-2-modo-consola-emisor--receptor-udp-separados)
- [Interfaz del Gemelo Digital](#-interfaz-del-gemelo-digital)
- [Créditos y Licencia](#-créditos-y-licencia)

---

## 🌌 Descripción General

Este proyecto implementa un **Gemelo Digital (Digital Twin)** para la simulación, recepción y análisis en tiempo real de la órbita de un satélite artificial en órbita baja terrestre (LEO) con características elípticas.

El sistema combina:
1. **Dinámica Orbital:** Simulación matemática exacta del movimiento kepleriano.
2. **Telecomunicaciones UDP:** Emisión y recepción de paquetes de telemetría con simulación de perturbaciones espaciales (tormentas solares y desvanecimiento de señal).
3. **Reconstrucción Numérica:** Estimación de posiciones perdidas en vuelo mediante interpolación polinómica recursiva con cotas de error en tiempo real.
4. **Visualización 3D Interactiva:** Interfaz en Streamlit y Plotly con proyección geocéntrica inercial (ECI), globo terrestre en malla 3D y KPIs reactivos actualizados a 2 Hz vía `st.fragment`.

---

## 📐 Fundamentación Matemática y Métodos Numéricos

El núcleo del software está sustentado en tres pilares del análisis numérico:

### 1. Ecuación de Kepler y Newton-Raphson
*Módulo:* [`kepler_solver.py`](kepler_solver.py)

Para determinar la posición de un cuerpo en una órbita elíptica en función del tiempo, se debe resolver la ecuación trascendente de Kepler:

$$M = E - e \cdot \sin(E)$$

donde:
- $M$: Anomalía media ($M = n(t - t_0)$).
- $E$: Anomalía excéntrica (incógnita).
- $e$: Excentricidad orbital ($0 \le e < 1$).

Definiendo la función raíz:
$$f(E) = E - e \cdot \sin(E) - M = 0$$
$$f'(E) = 1 - e \cdot \cos(E)$$

Se aplica el método de **Newton-Raphson** con convergencia cuadrática:
$$E_{k+1} = E_k - \frac{E_k - e \cdot \sin(E_k) - M}{1 - e \cdot \cos(E_k)}$$

**Estimación Inicial Óptima:**
- Para $e < 0.8$: $E_0 = M$
- Para $e \ge 0.8$: $E_0 = \pi$ si $M > \pi$, o $E_0 = M + e \cdot \sin(M)$

Una vez hallado $E$ con tolerancia de error residual $|f(E)| < 10^{-10}$, se obtienen las coordenadas orbitales en el plano del periapsis:
$$x = a(\cos E - e), \quad y = a\sqrt{1 - e^2}\sin E$$
$$\nu = 2 \arctan2\left(\sqrt{1+e}\sin\frac{E}{2}, \sqrt{1-e}\cos\frac{E}{2}\right)$$

---

### 2. Reconstrucción de Telemetría (Algoritmo de Neville)
*Módulo:* [`reconstruction.py`](reconstruction.py)

Durante ventanas de interferencia (ej. blackout atmosférico o radiación solar intensa), se pierden paquetes de telemetría enviados por UDP. Para estimar la posición del satélite en los instantes perdidos a partir de los últimos $n$ nodos válidos $(t_i, y_i)$, se utiliza el **Algoritmo de Neville**:

$$Q_{i, j} = \frac{(t - t_{i-j}) Q_{i, j-1} - (t - t_i) Q_{i-1, j-1}}{t_i - t_{i-j}}$$

**Ventajas clave:**
- Genera la tabla triangular de aproximaciones sucesivas de grado creciente de forma recursiva y computacionalmente eficiente ($O(n^2)$).
- **Cota intrínseca del error:** Proporciona una estimación directa del residuo local comparando el valor final con el del orden previo:
$$\epsilon \approx |Q_{n-1, n-1} - Q_{n-1, n-2}|$$

---

### 3. Cuadratura Numérica (Romberg y Simpson 1/3)
*Módulo:* [`orbital_metrics.py`](orbital_metrics.py)

Para calcular la longitud total del arco de órbita recorrido y la distancia acumulada a partir del perfil de velocidad $v(t)$:

$$s = \int_{t_0}^{t_f} v(t) \, dt$$

El proyecto implementa:
1. **Regla del Trapecio Compuesta:**
   $$T(h) = \frac{h}{2} \left[ f(a) + 2\sum_{i=1}^{n-1} f(a + ih) + f(b) \right]$$

2. **Regla de Simpson 1/3 Compuesta:**
   $$S(h) = \frac{h}{3} \left[ f(a) + 4\sum_{i=1,3,\dots}^{n-1} f(a + ih) + 2\sum_{i=2,4,\dots}^{n-2} f(a + ih) + f(b) \right]$$

3. **Integración de Romberg:**
   Aplica extrapolación sucesiva de Richardson a la regla del trapecio con duplicación de subintervalos:
   $$R_{i, j} = \frac{4^j R_{i, j-1} - R_{i-1, j-1}}{4^j - 1}$$
   Acelerando drásticamente el orden de convergencia de $O(h^2)$ a $O(h^{2(j+1)})$.

---

## 🏗️ Arquitectura del Sistema

El flujo de información y procesamiento opera según el siguiente diagrama:

```
                  ┌───────────────────────────────────────────────┐
                  │              satellite_uplink.py              │
                  │   - Genera trayectoria con kepler_solver.py   │
                  │   - Simula tormenta solar (pérdida UDP)       │
                  │   - Transmite JSON vía UDP Socket (:5006)     │
                  └───────────────────────┬───────────────────────┘
                                          │
                                UDP Packets (JSON)
                                (Puerto 5005)
                                          │
                  ┌───────────────────────▼───────────────────────┐
                  │              OrbitalTelemetryHub              │
                  │         (Thread en dashboard.py)              │
                  │   - Recepción continua no bloqueante          │
                  │   - Handshake bidireccional "GROUND_READY"    │
                  │   - Detección de pérdida de paquetes (gap ID) │
                  └───────────┬───────────────────────┬───────────┘
                              │                       │
           Si hay pérdida de paquetes            Métricas y Trayectoria
                              │                       │
                  ┌───────────▼───────────┐   ┌───────▼───────────────────┐
                  │   reconstruction.py   │   │    orbital_metrics.py     │
                  │  Algoritmo de Neville │   │  Cuadratura de Romberg    │
                  │  (Interpolación 2D)   │   │  (Longitud de arco/dist.) │
                  └───────────┬───────────┘   └───────┬───────────────────┘
                              │                       │
                              └───────────┬───────────┘
                                          │
                  ┌───────────────────────▼───────────────────────┐
                  │                 dashboard.py                  │
                  │  - Proyección Inercial 3D (i=51.6°, Ω=30°)    │
                  │  - Globo terráqueo wireframe + atmósfera     │
                  │  - KPIs en tiempo real (Pérdida, Error, Dist)│
                  │  - Refresco reactivo a 2 Hz (st.fragment)    │
                  └───────────────────────────────────────────────┘
```

---

## 📁 Estructura del Proyecto

```text
├── dashboard.py           # Dashboard principal (Streamlit + Plotly 3D + Hub UDP)
├── ground_station.py      # Estación terrena standalone en consola (listener UDP)
├── satellite_uplink.py    # Emisor satelital UDP (simulador de órbita y pérdidas)
├── kepler_solver.py       # Solucionador numérico de Kepler (Newton-Raphson)
├── reconstruction.py      # Interpolación de Neville y polinomios de Lagrange
├── orbital_metrics.py     # Cuadratura numérica (Trapecio, Simpson 1/3, Romberg)
├── requirements.txt       # Dependencias del proyecto
└── README.md              # Documentación técnica del proyecto
```

---

## ⚙️ Instalación y Requisitos

### Prerrequisitos
- **Python 3.10** o superior instalado en el sistema.
- Gestor de paquetes `pip`.

### 1. Clonar el repositorio
```bash
git clone https://github.com/tu-usuario/proyecto-metodos-numericos.git
cd "proyecto-metodos-numericos"
```

### 2. Crear y activar entorno virtual (Recomendado)
- **En Windows:**
  ```powershell
  python -m venv venv
  .\venv\Scripts\activate
  ```
- **En Linux / macOS:**
  ```bash
  python3 -m venv venv
  source venv/bin/activate
  ```

### 3. Instalar dependencias
```bash
pip install -r requirements.txt
```

*Contenido de `requirements.txt`:*
- `streamlit>=1.37.0`
- `plotly>=5.20.0`
- `pandas>=2.2.0`

---

## 🚀 Guía de Uso

### Opción 1: Dashboard Integrado (Recomendada)

La forma más sencilla e interactiva de utilizar el sistema es ejecutando directamente el dashboard:

```bash
streamlit run dashboard.py
```

1. Se abrirá automáticamente en tu navegador web (por defecto en `http://localhost:8501`).
2. El hub de telemetría interno se enlazará automáticamente al puerto UDP `5005`.
3. Haz clic en el botón superior **🚀 Lanzar Satélite**:
   - Se iniciará el emisor en un hilo en segundo plano.
   - Observarás la trayectoria orbital trazarse en 3D alrededor de la Tierra.
   - Al cruzar la zona de interferencia solar, el sistema detectará los paquetes perdidos y calculará automáticamente las coordenadas reconstruidas mediante **Neville**, mostrándolas como marcas rojas `×` junto con la tabla de auditoría de errores.

---

### Opción 2: Modo Consola (Emisor / Receptor UDP Separados)

Si deseas probar la comunicación de red distribuida en dos consolas independientes:

1. **Terminal 1 (Estación Terrena en consola):**
   ```bash
   python ground_station.py
   ```
   *Escuchará en el puerto 5005 e imprimirá los paquetes recibidos y las interpolaciones calculadas.*

2. **Terminal 2 (Satélite emisor):**
   ```bash
   python satellite_uplink.py
   ```
   *Establecerá el handshake con la estación y comenzará la transmisión de la órbita.*

---

## 🖥️ Interfaz del Gemelo Digital

El panel de control incluye:

- **Indicador de Enlace en Vivo:** Monitoreo dinámico del estado del socket UDP (`Enlace Activo`, `Standby` o `Sin Señal`).
- **4 Tarjetas de Métricas (KPIs):**
  - **Telemetría Recibida:** Total de paquetes recibidos e ID del próximo paquete esperado.
  - **Paquetes Perdidos:** Conteo acumulado y porcentaje de pérdida de paquetes.
  - **Interpolaciones Neville:** Total de puntos recuperados numéricamente y cota de residuo calculada ($\pm \text{km}$).
  - **Distancia Orbital:** Longitud de arco acumulada y resultado de la cuadratura de Romberg.
- **Visualizador 3D Geocéntrico Inercial:**
  - Modelo wireframe y malla translúcida de la Tierra ($R = 6371\text{ km}$).
  - Inclinación orbital configurada a $51.6^\circ$ (equivalente a la ISS) con rotaciones por RAAN y argumento de periapsis.
  - Trazado de telemetría real vs. puntos reconstruidos por Neville.
  - Persistencia de cámara (zoom y rotación conservados entre refrescos vía `uirevision`).
- **Descomposición 2D en Plano Orbital:** Gráficos temporales de coordenadas $X(t)$ e $Y(t)$.
- **Tabla de Auditoría Numérica:** Registro tabular de cada punto interpolado, con timestamps y cotas de error formalmente calculadas.

---

## 👥 Créditos y Licencia

Desarrollado como proyecto aplicativo para la materia de **Métodos Numéricos**. Demuestra la aplicación rigurosa de métodos de raíces, interpolación y cuadratura numérica a problemas reales de ingeniería aeroespacial y telecomunicaciones.

Distribuido bajo la licencia [MIT](LICENSE).
