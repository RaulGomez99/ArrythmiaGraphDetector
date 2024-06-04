# Arrythmia Graph Detection

## Descripción

Esta librería es de un Trabajo de Final de Grado de la Universitat de València para el Grado en Ciencia de Datos. El objetivo de este trabajo es la detección de arritmias mediante la detección de reentradas haciendo un analisis anatómicos, modelizando esta mediante grafos y aplicando algoritmos de detección de ciclos en grafos.

## Instalación

Instalar mediante pip y el enlace de GitHub:

```bash
pip install git+https://github.com/RaulGomez99/ArrythmiaGraphDetector
```

## Uso

La librería hace uso de archivos OBJ para la representación de los modelos anatómicos (en cado de querer detectar reentradas funcionales hacerlo con una resolución aproximada de 2.3 milímetros de distancia entre nodos para que su uso sea eficiente y no superior de limite superior de distancia de 20 milímetros). También es necesario tener un CSV con la información de los nodos (su dirección de fibras, material y modelo). Y por último es necesario un archivo CSV con las velocidades de conducción en función de material y modelo.

El diccionario de indices es el siguiente:

 Material:
 - 1 -> Right Atrium
- 2 -> Crista Terminalis
- 3 -> Pulmonary Veins
- 4 -> Bachmann's Bundle
- 5 -> Interatrial Septum
- 6 -> Sinoatrial Node
- 7 -> Left Fibrous Ring
- 8 -> Coronary Sinus
- 9 -> Mitral Valve/Left Atrial Appendage
- 10 -> Fossa Ovalis

Model:
- 191 -> Right Atrium
- 192 -> Bachmann's Bundle
- 193 -> Right Atrial Appendage
- 194 -> Tricuspid Valve
- 195 -> Left Atrium
- 196 -> Pulmonary Veins
- 197 -> Left Atrial Appendage
- 198 -> Mitral Valve


## Funciones

- `crear_grafo`: Función de lectura de un archivo obj para obtener los datos necesarios

- `crear_grafo_vtk`: Función para crear un grafo a partir de un archivo obj y un archivo csv del VTK con los datos de cada punto

- `pintar_puntos`: Función para pintar los puntos de un grafo y mostrar reentradas anatómicas

- `pintar_puntos_rotores`: Función para pintar los puntos de un grafo y mostrar reentradas funcionales con mapa de calor en función del tiempo del camino

- `pintar_puntos_rotores_binario`: Función para pintar los puntos de un grafo y mostrar reentradas funcionales con mapa de calor binario en función si hay una posible reentrada o no

- `calcular_tiempo_camino`: Función para calcular el tiempo de un camino en función de las velocidades de conducción

- `calcular_distancia_camino`: Función para calcular la distancia de un camino

- `calcular_tiempo_maximo_punto`: Función para calcular el tiempo máximo en recorrer cualquier camino que pasa por ese punto

- `detectar_rotores`: Función para detectar reentradas funcionales en un grafo, mediante topes de tiempo y distancia

- `guardar_caminos`: Función para guardar los caminos en un archivo CSV

- `cargar_caminos`: Función para cargar los caminos de un archivo CSV

- `obtener_velocidades_csv`: Función para obtener las velocidades de conducción en función de material y modelo

## Ejemplo

```python
import ArrythmiaGraphDetector as agd
import numpy as np

[puntos, triangulos, _, matriz_adyacencia, reentradas_anatomicas] = agd.crear_grafo("auricula.obj")
agd.pintar_puntos(puntos, matriz_adyacencia, caminos = reentradas_anatomicas, show_index = False)

datos_puntos = agd.crear_grafo_vtk("auricula.obj", "datos_auricula.csv")[4]

puntos = np.array(puntos) # Importante porque si no, no se puede calcular velocidades
velocidades = agd.obtener_velocidades_csv("datos_velocidades.csv", datos_puntos, 1)

reentradas_funcionales = agd.detectar_rotores(matriz_adyacencia, puntos, velocidades, limites_espacio = [12, 20], limites_tiempo = [200, 99999])

tiempos = [int(round(agd.calcular_tiempo_camino(camino, velocidades, puntos)*1000)) for camino in reentradas_funcionales]
tiempo_puntos = agd.calcular_tiempo_maximo_punto(puntos, tiempos, reentradas_funcionales)

agd.pintar_puntos_rotores(puntos, tiempo_puntos, triangulos)

agd.guardar_caminos(reentradas_funcionales, "caminos.csv")
caminos = agd.cargar_caminos("caminos.csv")

```
## Autor

- Raúl Gómez López

## Contacto

- Raúl Gómez López: gomezlopezraul1999@gmail.com

## Licencia

 No License