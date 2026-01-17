# Mecánica Molecular (MM), Dinámica Molecular (DM) y Materiales Compuestos Colombianos

El análisis de la energía de materiales complejos mediante **campos de fuerza (CF)** es el núcleo de la **Mecánica Molecular (MM)** y la **Dinámica Molecular (DM)**. Pensando en un contexto de investigación de **materiales compuestos colombianos**, como aquellos reforzados con fibras naturales como el **fique** o la **guadua**, el uso de CF se vuelve crucial para evaluar su estabilidad y propiedades. 

---

## La Función de la Energía en Materiales Complejos

En un sistema molecular, la energía total potencial ($E_{potencial}$) se calcula utilizando el **Campo de fuerza**. Esta energía es la suma de términos que describen las interacciones entre los átomos:

$E_{potencial} = E_{enlazada} + E_{no\ enlazada}$

Donde:

1.  **$E_{enlazada}$** (Intramolecular): Contribución de átomos que están directamente conectados por enlaces químicos.
2.  **$E_{no\ enlazada}$** (Intermolecular): Contribución de átomos que están espacialmente cerca, pero no conectados por enlaces.

Al evaluar esta energía para diferentes configuraciones (geometrías) de un material, podemos predecir su **estabilidad (Minimización de Energía)**, su **comportamiento dinámico (Dinámica Molecular)** y cómo se integran sus componentes.

---

## Campos de Fuerza Típicos para Materiales (Considerando Colombia)

Dada la naturaleza de los **materiales compuestos** (polímeros, fibras naturales/sintéticas y sus interfaces), se requieren CF que manejen sistemas grandes y heterogéneos.

| Tipo de Campo de Fuerza               | Enfoque                                                      | Aplicación en Materiales Compuestos                          | Ejemplos Comunes                   |
| :------------------------------------ | :----------------------------------------------------------- | :----------------------------------------------------------- | :--------------------------------- |
| **Clásico / *All-Atom***              | Describe cada átomo de forma explícita. Ofrecen alta precisión. | Interacciones detalladas entre la matriz polimérica y la superficie de las fibras (ej. **celulosa** de fique o guadua) o el estudio de defectos cristalinos. | **AMBER**, **CHARMM**, **OPLS-AA** |
| **Reactivo (*Reactive*)**             | Permite la **ruptura y formación de enlaces** durante la simulación. | Modelado de la **degradación térmica** de la matriz, la **interacción química** en la interfaz o la **catálisis** en procesos de síntesis. | **ReaxFF**, **EVB**                |
| **Grano Grueso (*Coarse-Grained*)**   | Agrupan varios átomos en una sola "esfera" (partícula). Sacrifican detalle por **eficiencia computacional**. | Simulación del **comportamiento a escala macro** (ej. grandes agregados, viscosidad, autoensamblaje) o la interacción de polímeros de alto peso molecular. | **MARTINI**, **SIRAH**             |
| **Basado en Machine Learning (MLFF)** | Utilizan el aprendizaje automático para reproducir la precisión de la Química Cuántica a la velocidad de la Mecánica Molecular. Ideal para sistemas novedosos y complejos donde no existen parámetros pre-calibrados, ofreciendo alta precisión en la energía de la interfaz. | **ANI**, **SchNet**                |

---

## Aplicación en Problemas Reales de Materiales en Colombia

Para un material compuesto hecho, por ejemplo, de **Matriz de Polietileno (HDPE)** reforzada con **Fibras de Fique/Guadua** (como se investiga en varias universidades colombianas):

1.  **Evaluación de la Interfase:**
    * **CF: CHARMM/OPLS-AA (Clásico)**
    * **Objetivo:** Simular la energía de adhesión entre la matriz polimérica (HDPE) y la celulosa/lignina de la fibra natural. Una **menor energía de interfaz** indica una unión más fuerte y, por ende, mejores propiedades mecánicas del compuesto final.
2.  **Propiedades Mecánicas (Nanoindentación):**
    * **CF: ReaxFF (Reactivo) o LAMMPS (Plataforma)**
    * **Objetivo:** Simular cómo reacciona el material a una fuerza localizada (como una nanoindentación). Esto implica modelar la deformación y potencialmente la ruptura de enlaces, para obtener propiedades como la dureza y el módulo de elasticidad.
3.  **Comportamiento a Largo Plazo (DM):**
    * **CF: MARTINI (Grano Grueso)**
    * **Objetivo:** Estudiar cómo la estructura general del polímero cambia o se degrada con el tiempo o la temperatura (ej. envejecimiento acelerado), algo vital para la durabilidad de un producto en el clima tropical.

El estudio titulado "Interfaz Gráfica de Usuario para la Simulación por Dinámica Molecular de Películas Delgadas" es un ejemplo de investigación colombiana que implementa la Dinámica Molecular para simular la dureza de materiales.

