# Molecular Modelling

El **modelado molecular** es un término general que describe un conjunto de técnicas computacionales utilizadas para simular y experimentar con las estructuras físicas y propiedades químicas de las moléculas y materiales.

Los pasos generales para configurar el entorno de los ejemplos que usan ASE están en la [guía de instalación](docs/installation/readme.md).

## Alcance computacional del repositorio

Aunque este README presenta distintas áreas del modelado molecular, los ejemplos computacionales disponibles se concentran principalmente en cálculos de estructura electrónica de moléculas aisladas: optimización geométrica y análisis de frecuencias con DFT. La mayoría de los scripts Python usan ASE como interfaz y NWChem como calculador; también hay un ejemplo para dióxido de carbono con Quantum ESPRESSO y herramientas de análisis/postprocesamiento.

Los ejemplos son scripts independientes, no una aplicación unificada ni un flujo automatizado común. Su ejecución requiere Python, ASE y las bibliotecas que importe cada script, además del programa externo correspondiente (por ejemplo, NWChem o Quantum ESPRESSO). También hay que ejecutar cada script desde la carpeta que contiene sus archivos de entrada, ya que varios usan rutas relativas.

Por tanto, la dinámica molecular de líquidos, las propiedades de transporte y los métodos multiescala descritos más abajo son contexto científico, no capacidades cubiertas de forma general por los scripts actuales. [examples/butane.v2.md](examples/butane.v2.md) propone un futuro protocolo GROMACS/OPLS-AA para butano líquido; es una guía de trabajo y no un flujo ejecutable validado incluido en los ejemplos actuales.

## ¿En qué consiste el modelado molecular?

El objetivo principal es representar y visualizar las **formas tridimensionales** de las moléculas y simular cómo se mueven e interactúan, ya que son demasiado pequeñas para ser vistas con microscopía clásica.

Generalmente, el proceso en simulaciones computarizadas sigue tres pasos:

1. **Definición de las Características:** Se establecen las características preliminares de la molécula, incluyendo el número de átomos y la forma en que los elementos se unen o enlazan.
2. **Representación 3D:** Se definen las posiciones relativas de los átomos constituyentes para crear una representación tridimensional de la molécula.
3. **Simulación de Interacciones:** Se simulan las interacciones químicas y moleculares para predecir reactividad, estabilidad y cambios estructurales.

## Tipos de Técnicas Principales

El modelado molecular abarca varios métodos, que varían en la forma en que tratan a los átomos:

* **Mecánica Molecular (MM):**
  * Trata a los **átomos** como las unidades más pequeñas (enfoque de "bola y resorte").
  * Utiliza **campos de fuerza** (expresiones matemáticas empíricas) para estimar la energía potencial del sistema.
  * Es útil para sistemas grandes (como proteínas o polímeros) porque es menos intensivo computacionalmente.

* **Química Cuántica (QC) / Métodos *Ab Initio*:**
  * Modela explícitamente **protones, neutrones y electrones**.
  * Utiliza los principios de la mecánica cuántica para calcular propiedades con mayor precisión.
  * Es computacionalmente más exigente y se limita a sistemas más pequeños.

## El Puente Intermedio: Métodos Multiescala

Entre la Mecánica Molecular (átomos explícitos) y la Química Cuántica (electrones explícitos) existe una categoría adicional que el modelado molecular clásico suele dejar fuera del mapa, pero que es clave cuando el sistema de interés es demasiado grande incluso para MM, o cuando se necesita conectar la escala atómica con fenómenos macroscópicos:

* **Modelos de Grano Grueso (*Coarse-Grained*):**
  * Agrupan varios átomos en una sola "pseudo-partícula" o *bead*, reduciendo drásticamente los grados de libertad del sistema.
  * Sacrifican detalle atómico a cambio de acceder a escalas espaciales y temporales mucho mayores (miles de moléculas, microsegundos en vez de picosegundos).

* **Dinámica de Partículas Disipativas (*DPD - Dissipative Particle Dynamics*):**
  * Extiende la idea de grano grueso incorporando fuerzas conservativas, disipativas y aleatorias entre partículas, preservando el comportamiento hidrodinámico del sistema.
  * Es especialmente útil para modelar mesoescalas de materiales blandos, interfaces y procesos de agregación o difusión que la Dinámica Molecular clásica no puede alcanzar por su costo computacional.

* **Método de Elementos Discretos (*DEM - Discrete Element Method*):**
  * Trata las unidades del sistema (partículas, granos, agregados) como cuerpos rígidos o semi-rígidos que interactúan mecánicamente (contacto, fricción, cohesión).
  * Conecta la escala mesoscópica con el comportamiento estructural y mecánico observable a nivel macroscópico, típico en simulaciones de materiales granulares o compuestos.

Estos métodos no reemplazan a MM ni a QC, sino que actúan como **puente entre escalas**: permiten propagar información desde interacciones atómicas (obtenidas vía MD o incluso QC) hacia la respuesta mecánica de un material a escalas de ingeniería, en un enfoque típicamente llamado **modelado multiescala**.

## Técnicas de Simulación Clave

Dentro del modelado molecular, hay métodos clave para investigar el comportamiento dinámico o estático de un sistema:

* **Minimización de Energía:**
  * Busca la **estructura molecular más estable** (el estado de energía potencial más bajo).
  * Proporciona una imagen estática útil para comparar estados de sistemas similares.

* **Dinámica Molecular (DM):**
  * Calcula el **comportamiento del sistema en función del tiempo**.
  * Resuelve las leyes del movimiento de Newton para obtener las trayectorias de los átomos en el espacio y el tiempo.
  * Proporciona información sobre procesos dinámicos e incluye los efectos de la **temperatura**.

* **Acoplamiento Molecular (*Molecular Docking*):**
  * Predice la **orientación preferida de una molécula pequeña** (ligando) cuando se une a una segunda molécula (receptor, como una proteína).
  * Es fundamental en el diseño de fármacos para predecir la afinidad de unión.

## Aplicaciones

El modelado molecular es una herramienta rutinaria en diversas áreas científicas:

* **Investigación Farmacéutica (Diseño de Fármacos):**
  * Diseño de nuevos fármacos (ej. el desarrollo de Saquinavir, uno de los primeros medicamentos anti-VIH diseñados con programas bioinformáticos).
  * Refinamiento y optimización de estructuras de medicamentos.

* **Biología Computacional y Bioinformática:**
  * Investigación del **plegamiento de proteínas**, catálisis enzimática, estabilidad proteica y reconocimiento molecular.

* **Ciencia de Materiales:**
  * Descubrimiento y desarrollo de nuevas sustancias sintéticas.

* **Modelado Multiescala de Composites y Materiales de Construcción:**
  * Combinación de MD (interacciones atómicas), DPD (mesoescala) y DEM (respuesta mecánica macroscópica) para predecir el comportamiento de materiales compuestos complejos, como matrices cementicias reforzadas con polímeros.
