#from .Functions import pintar_puntos, pintar_puntos_rotores, pintar_puntos_rotores_binario, calcular_tiempo_camino, calcular_distancia_camino, calcular_tiempo_maximo_punto, detectar_rotores, crear_grafo, crear_grafo_vtk, guardar_caminos, cargar_caminos, obtener_velocidades_csv

#%%
from math import sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import plotly.graph_objects as go
import networkx as nx
import ast

# Constantes utilizadas
LIMITE_DIST_INF = 8
LIMITE_DIST_SUP = 20

LIMITE_TMP_INF = 100
LIMITE_TMP_SUP = 99999

VALOR_REENTRADA = 200
VALOR_MED_REENTRADA = 120

DICCIONARIO_MATERIALS = {"1": "RA", "2": "CT", "3": "PVS", "4": "BB", "5": "IST", "6": "SAN", "7": "LFO", "8": "CS", "9":"MV/LAA", "10": "FO"}
DICCIONARIO_MODELS = {"191": "RA", "192": "BB", "193": "RAA", "194": "TV", "195": "LA", "196": "PVS", "197":"LAA", "198":"MV"}


# Funciones internas
# ----------------------------------------------------------------------------------
def calcular_distancia(punto1:list, punto2:list)->float:
    """
    Función interna para calcular distancias entre dos puntos
    
    Args:
        punto1 (list): Punto de origen
        punto2 (list): Punto de destino

    Returns:
        float: Distancia euclidea entre los puntos
    """
    return sqrt((punto1[0] - punto2[0])**2+(punto1[1] - punto2[1])**2)


def transformar_nombres(puntos:list)->list:
    """
    Función interna para transformar los nombres de los puntos a strings para crear el grafo

    Args:
        puntos (list): Lista de los puntos

    Returns:
        list: Lista con los nombre de los indices de los puntos pasados a string
    """
    nombres = []
    for i in range(len(puntos)):
        nombres.append(str(i))
    return nombres

def transformar_relaciones(matriz_adyacencia:list)->list:
    """
    Función interna para transformar la matriz de adyacencia en una lista de relaciones

    Args:
        matriz_adyacencia (list): Matriz de adyacencia

    Returns:
        list: Relaciones en formato grafo
    """
    relaciones = []
    for i in range(matriz_adyacencia.shape[0]):
        for j in range(matriz_adyacencia.shape[1]):
            if(matriz_adyacencia[i,j] != 0):
                relaciones.append((str(i), str(j), matriz_adyacencia[i,j]))
    return relaciones

def transformar_coordenadas(puntos:list)->dict:
    """
    Función interna para transformar las coordenadas de los puntos en un diccionario
    
    Args:
        puntos (list): Lista de los puntos
        
    Returns:
        dict: Diccionario con las coordenadas de los puntos
    """
    coordenadas = dict()
    for i in range(len(puntos)):
        coordenadas[str(i)] = (puntos[i][0], puntos[i][1])
    return coordenadas

def son_perpendiculares(v:list, w:list)->float:
    """
    Función interna para calcular el grado de perpendicularidad
    
    Args:
        v (list): Vector 1
        w (list): Vector 2
        
    Returns:
        float: Grado de perpendicularidad (0, 90)
    """
    producto_punto = np.dot(v, w)
    norma_v = np.linalg.norm(v)
    norma_w = np.linalg.norm(w)
    angulo = np.degrees(np.arccos(producto_punto / (norma_v * norma_w)))
    if angulo > 90: angulo = 180 - angulo
    return angulo/90

def detectar_camino_rotor(matriz_adyacencia:list, puntos:list, punto:int, limites_espacio:list, limites_tiempo:list, camino:list = [], unic:list = [], ant:list = None, caminos:list = [], velocidades:list = [])->None:
    """
    Función interna para detectar los caminos de los rotores

    Args:
        matriz_adyacencia (list): Matriz de adyacencia
        puntos (list): Lista de puntos
        punto (int): Punto en el que se encuentra
        limites_espacio (list): Límites de espacio [limite_inferior, limite_superior]
        limites_tiempo (list): Límites de tiempo [limite_inferior, limite_superior]
        camino (list, optional): Camino actual. Defaults to []
        unic (list, optional): Caminos únicos. Defaults to []
        ant ([type], optional): Punto anterior. Defaults to None
        caminos (list, optional): Caminos de rotores encontrados. Defaults to []
        velocidades: Velocidades de los puntos [[vector], velocidad (mm/s), penalización]
        
    Returns:
        None
    """
    if camino == []:
        camino = [punto]
    
    dist = calcular_distancia_camino(camino, puntos)
    if dist > limites_espacio[1]: return
    if dist + calcular_distancia(puntos[punto], puntos[camino[0]]) > limites_espacio[1]: return
    
    if velocidades != []:
        temp = calcular_tiempo_camino(camino, velocidades, puntos)
        temp = int(round(temp*1000))
        if temp > limites_tiempo[1]: return
    
    
    conexiones = [i for i, x in enumerate(matriz_adyacencia[punto]) if x != 0]
    for i in conexiones:
        if i == ant: continue
        if i in camino:
            if i != camino[0]: continue
            camino2 = list(camino)
            camino2.append(i)
            dist2 = calcular_distancia_camino(camino2, puntos)
            if velocidades != []:
                temp2 = calcular_tiempo_camino(camino, velocidades, puntos)
                temp2 = int(round(temp2*1000))
            else: temp2 = limites_tiempo[0]
            
            if (limites_espacio[0] <= dist2 <= limites_espacio[1]) and (limites_tiempo[0] <= temp2 <= limites_tiempo[1]):
                camino2 = camino2[camino2.index(i):]
                camino2 = camino2[:-1]
                camino2.sort()
                if camino2 not in unic:
                    unic.append(camino2)
                    caminos.append(camino + [i])
                continue
            else: continue
        camino2 = list(camino)
        camino2.append(i)
        detectar_camino_rotor(matriz_adyacencia, puntos, punto = i, limites_espacio = limites_espacio, limites_tiempo = limites_tiempo, camino = camino2, unic = unic, ant = punto, caminos = caminos, velocidades = velocidades)
     
def retorna_color_relacion(tiempos:list, valor_max:int, valor_med:int)->str:
    """
    Función interna para retornar el color de una relación en base a los tiempos

    Args:
        tiempos (list): Tiempos de los nodos implicados
        valor_max (int): Tiempo (ms) que se considera rojo en el mapa de calor
        valor_med (int): Tiempo (ms) en el que se considera amarillo en el mapa de calor

    Returns:
        str: Color formato HTML
    """
    return valor_a_color_html(max(tiempos), valor_max, valor_med)

def retorna_color_poligono(tiempos:list, max_tiempo:int, med_tiempo:int)->str:
    """
    Función interna para retornar el color de un polígono en base a los tiempos
    
    Args:
        tiempos (list): Lista de tiempos de los nodos
        max_tiempo (int): Tiempo (ms) que se considera rojo en el mapa de calor
        med_tiempo (int): Tiempo (ms) en el que se considera amarillo en el mapa de calor
    Returns:
        str: Color en formato HTML
    """
    return valor_a_color_html(max(tiempos), max_tiempo, med_tiempo)

def valor_a_color_html(valor:int, valor_max:int = VALOR_REENTRADA, valor_med:int = VALOR_MED_REENTRADA):
    """
    Función interna para transformar un valor en un color en formato HTML

    Args:
        valor (int): Tiempo a transformar
        valor_max (int, optional): Tiempo (ms) que se considera rojo. Defaults to VALOR_REENTRADA.
        valor_med (int, optional): Tiempo (ms) que se considera amarillo. Defaults to VALOR_MED_REENTRADA.

    Returns:
        _type_: _description_
    """
    
    cmap_dict = {
        (0, valor_med): ('#0000FF', '#FCFF0E'),           # Azul a amarillo
        (valor_med, valor_max): ('#FCFF0E', '#FF0000'),   # Amarillo a rojo 
    }
    if valor >= valor_max: return '#FF0000'
    for rango, colores in cmap_dict.items():
        if rango[0] <= valor < rango[1]:
            inicio, fin = colores
            inicio_rgb = mcolors.hex2color(inicio)
            fin_rgb = mcolors.hex2color(fin)
            fraccion = (valor - rango[0]) / (rango[1] - rango[0])
            color_interpolado = tuple((1 - fraccion) * inicio_rgb[i] + fraccion * fin_rgb[i] for i in range(3))
            return mcolors.rgb2hex(color_interpolado)
        
def punto_in_rotor(punto:int, rotores:list)->bool:
    """
    Función interna para saber si un punto está en algún rotor

    Args:
        punto (int): Punto a comprobar
        rotores (list): Lista de rotores

    Returns:
        bool: Devuelve si está en algún rotor
    """
    for rotor in rotores:
        if punto in rotor: return True
    return False

def filtrar_rotores(caminos:list, puntos:list, tmp_min:int, tmp_max:int, dist_min:float, dist_max:float, velocidades:list)->list:
    """
    Función interna para filtrar los caminos en base a los tiempos y distancias (solo filtra no encuentra nuevos caminos)

    Args:
        caminos (list): Lista de caminos
        puntos (list): Lista de puntos
        tmp_min (int): Tiempo en milisegundos minimo para el filtrado
        tmp_max (int): Tiempo en milisegundos máximo para el filtrado
        dist_min (float): Distancia mínima en milímetros para el filtrado
        dist_max (float): Distancia máxima en milímetros para el filtrado
        velocidades (list): Velocidades de los puntos [[vector], velocidad (mm/s), penalización]

    Returns:
        list: Lista de rotores filtrados
    """
    distancias = [calcular_distancia_camino(camino, puntos) for camino in caminos]
    tiempos = [int(round(calcular_tiempo_camino(camino, velocidades, puntos)*1000)) for camino in caminos]
    caminos_filtrados = []
    for i in range(len(caminos)):
        if dist_min <= distancias[i] <= dist_max and tmp_min <= tiempos[i] <= tmp_max:
            caminos_filtrados.append(caminos[i])
            
    return caminos_filtrados
   
# Funciones de pintar
# ----------------------------------------------------------------------------------       

def pintar_puntos(puntos:list, matriz_adyacencia:list, caminos:list = None, show_index:bool = False, tiempos:list = []):
    """
    Función para pintar los puntos y las relaciones entre ellos y los caminos si se quiere
    
    Args:
        puntos (list): Lista de los puntos
        matriz_adyacencia (list): Matriz de adyacencia
        caminos (list, optional): Lista de caminos. Defaults to None.
        show_index (bool, optional): Mostrar los indices de los puntos. Defaults to False.
        tiempos (list, optional): Lista de tiempos de los caminos. Defaults to [].
        
    Returns:
        None
    """
    grafo = nx.Graph()
    nombres = transformar_nombres(puntos)
    relaciones = transformar_relaciones(matriz_adyacencia)
    coordenadas = transformar_coordenadas(puntos)
    
    grafo.add_nodes_from(nombres)
    grafo.add_weighted_edges_from(relaciones)
    
    colors = ["black" for _ in range(len(relaciones))]
    edges_widht = [1 for _ in range(len(relaciones))]
    node_colors = ["blue" for _ in range(len(puntos))]
    
    if(caminos != None):
        for c, camino in enumerate(caminos):
            rels = list(grafo.edges)
            for i in range(len(camino)-1):
                if((str(camino[i]), str(camino[i+1])) in rels):
                    pos = rels.index((str(camino[i]), str(camino[i+1])))
                else:
                    pos = rels.index((str(camino[i+1]), str(camino[i])))
                edges_widht[pos] = 5
                if tiempos == []: colors[pos] = "r"
                else:
                    if tiempos[c] > VALOR_REENTRADA: colors[pos] = "r"
                    elif tiempos[c] > VALOR_MED_REENTRADA: colors[pos] = "y"
                    else: colors[pos] = "g"
        
    nx.draw(grafo, pos = coordenadas, node_color = node_colors, with_labels= show_index, edge_color = colors, node_size = 5, font_size = 8, width = edges_widht)
    
def pintar_puntos_rotores(puntos:list, tiempos:list, triangulos:list, max_tiempo:int = VALOR_REENTRADA, med_tiempo:int = VALOR_MED_REENTRADA)->None:
    """
    Función para pintar un mapa de calor del corazón según lo propenso que se sea a los rotores

    Args:
        puntos (list): Lista de puntos
        tiempos (list): Lista del máximo tiempo de cada punto
        triangulos (list): Triangulos del obj
        max_tiempo (int, optional): Valor (ms) para que sea rojo. Defaults to VALOR_REENTRADA.
        med_tiempo (int, optional): Valor (ms) para que sea amarillo. Defaults to VALOR_MED_REENTRADA.

    Returns:
        None
    """
    fig = go.Figure()
    x_lim = [min([punto[0] for punto in puntos]), max([punto[0] for punto in puntos])]
    y_lim = [min([punto[1] for punto in puntos]), max([punto[1] for punto in puntos])]
    z_lim = [min([punto[2] for punto in puntos]), max([punto[2] for punto in puntos])]

    x = [punto[0] for punto in puntos]
    y = [punto[1] for punto in puntos]
    z = [punto[2] for punto in puntos]
    
    
    i2 = []
    j = []
    k = []
    triangle_colors = [
       [0.0, '#0000ff'],  # Azul
       [med_tiempo/max_tiempo, '#fcff0e'],  # Amarillo
       [1.0, '#ff0000']   # Rojo
        
    ]
    intensidades = [min(t/200, 1) for t in tiempos]
    intensidades.append(1)
    intensidades.append(0)
    for i in range(len(triangulos)):
        i2.append(triangulos[i][0]-1)
        j.append(triangulos[i][1]-1)
        k.append(triangulos[i][2]-1)
        
        
    fig.add_trace(go.Mesh3d(x=x, y=y, z=z, i=i2, j=j, k=k, 
                            colorscale = triangle_colors, opacity=1, 
                            intensity = intensidades, name="Rotores",
                            colorbar=dict(
                                title="Tiempo de rotor",
                                tickmode="array",
                                tickvals=[0.001, med_tiempo/max_tiempo, 1],  
                                ticktext=["0 ms", str(med_tiempo)+" ms", str(max_tiempo)+" ms"] 
                            )))
    fig.update_layout(scene=dict(
        xaxis=dict(range = x_lim, visible=False),
        yaxis=dict(range = y_lim, visible=False),
        zaxis=dict(range = z_lim, visible=False),
    ))
    fig.show()
    
def pintar_puntos_rotores_binario(puntos:list, puntos_in_rotor:list, triangulos:list)->None:
    """
    Función para pintar un mapa de calor del corazón según si hay posibles rotores o no
    
    Args:
        puntos (list): Lista de puntos
        puntos_in_rotor (list): Lista de puntos que están en un rotor
        triangulos (list): Triangulos del obj
        max_tiempo (int, optional): Valor (ms) para que sea rojo. Defaults to VALOR_REENTRADA.
        med_tiempo (int, optional): Valor (ms) para que sea amarillo. Defaults to VALOR_MED_REENTRADA.

    Returns:
        None
    """
    fig = go.Figure()
    x_lim = [min([punto[0] for punto in puntos]), max([punto[0] for punto in puntos])]
    y_lim = [min([punto[1] for punto in puntos]), max([punto[1] for punto in puntos])]
    z_lim = [min([punto[2] for punto in puntos]), max([punto[2] for punto in puntos])]

    x = [punto[0] for punto in puntos]
    y = [punto[1] for punto in puntos]
    z = [punto[2] for punto in puntos]
    
    
    i2 = []
    j = []
    k = []
    triangle_colors = [
       [0.0, '#22E103'],  # Verde
       [1.0, '#ff0000']   # Rojo
        
    ]
    intensidades = [0 for _ in range(len(puntos))]
    intensidades.append(1)
    intensidades.append(0)
    for i in puntos_in_rotor:
        intensidades[i] = 1
    for i in range(len(triangulos)):
        i2.append(triangulos[i][0]-1)
        j.append(triangulos[i][1]-1)
        k.append(triangulos[i][2]-1)
        
        
    fig.add_trace(go.Mesh3d(x=x, y=y, z=z, i=i2, j=j, k=k, 
                            colorscale = triangle_colors, opacity=1, 
                            intensity = intensidades, name="Rotores",
                            colorbar=dict(
                                title="Salud del corazón",
                                tickmode="array",
                                tickvals=[0, 0.001, 1],  
                                ticktext=["", "Sano", "Proarritmico"]  
                            )))
    fig.update_layout(scene=dict(
        xaxis=dict(range = x_lim, visible=False),
        yaxis=dict(range = y_lim, visible=False),
        zaxis=dict(range = z_lim, visible=False),
    ))
    
    # fig.add_trace(go.Scatter3d(
    #     x=[0, 0],  # coordenadas x para la línea
    #     y=[-20, 20],  # rango de coordenadas y para la línea
    #     z=[0, 0],  # rango de coordenadas z para la línea
    #     mode='lines',
    #     line=dict(color='black', width=8),
    #     name='Línea vertical en x=0'
    # ))
    fig.show()

# Funciones de calculos de caminos
# ----------------------------------------------------------------------------------
def calcular_tiempo_camino(camino:list, velocidades:list, puntos:list)->float:
    """
    Función para calcular el tiempo de un camino en base a las velocidades de los puntos
    
    Args:
        camino (list): Lista de puntos del camino
        velocidades (list): Lista de velocidades de los puntos [[vector], velocidad (mm/s), penalización]
        puntos (list): Lista de puntos
        
    Returns:
        float: Tiempo en segundos del camino
    """
    tiempo = 0
    for i in range(len(camino)-1):
        punto = puntos[camino[i]]
        punto2 = puntos[camino[i+1]]
        porc_pen = son_perpendiculares(punto2-punto, velocidades[camino[i]][0])
        pen = porc_pen*(velocidades[camino[i]][2]-1)+1
        veloc = pen*velocidades[camino[i]][1]
        dist = calcular_distancia(punto, punto2)
        tiempo += dist/veloc
    return tiempo

def calcular_distancia_camino(camino:list, puntos:list)->float:
    """
    Función interna para calcular la distancia de un camino

    Args:
        camino (list): Lista de puntos del camino
        puntos (list): Lista de puntos

    Returns:
        float: Distancia del camino en milímetros
    """
    distancia = 0
    for i in range(len(camino)-1):
        distancia += calcular_distancia(puntos[camino[i]], puntos[camino[i+1]])
    return distancia

def calcular_tiempo_maximo_punto(puntos:list, tiempos:list, rotores:list)->float:
    """
    Función para calcular el tiempo máximo de un punto en base a los caminos

    Args:
        puntos (list): Lista de puntos
        tiempos (list): Lista de tiempos de los puntos en milisegundos
        rotores (list): Lista de reentradas funcionales

    Returns:
        float: Lista con los tiempos máximod de un punto en milisegundos
    """
    tiempos_punto = []
    for i in range(len(puntos)):
        rotores_in_punto = [i for camino in rotores if i in camino]
        tiempo_punto = 0
        for rotor in rotores_in_punto:
            tiempo_punto = max(tiempo_punto, tiempos[rotor])
        
        tiempos_punto.append(tiempo_punto)
    return tiempos_punto


# Funciones de rotores
# ----------------------------------------------------------------------------------
def detectar_rotores(matriz_adyacencia:list, puntos:list, velocidades:list, limites_espacio:list = [LIMITE_DIST_INF, LIMITE_DIST_SUP], limites_tiempo:list = [LIMITE_TMP_INF, LIMITE_TMP_SUP])->list:
    """
    Función para detectar los rotores en un grafo con su matriz de adyacencia
    
    Args:
        matriz_adyacencia (list): Matriz de adyacencia
        puntos (list): Lista de puntos
        velocidades (list): Lista de velocidades de los puntos [[vector], velocidad (mm/s), penalización]
        limites_espacio (list, optional): Límites de espacio [limite_inferior, limite_superior]. Defaults to [LIMITE_DIST_INF, LIMITE_DIST_SUP].
        limites_tiempo (list, optional): Límites de tiempo [limite_inferior, limite_superior]. Defaults to [LIMITE_TMP_INF, LIMITE_TMP_SUP].
    
    Returns:
        list: Lista de caminos que forman los rotores
    """
    caminos = []
    unic = []
    puntos = np.array(puntos)
    for i in range(len(matriz_adyacencia[0])):
        print(str(i) + " de " + str(len(matriz_adyacencia[0])) + " => " + str(round(i/len(matriz_adyacencia[0])*100, 2)) + "%")
        detectar_camino_rotor(matriz_adyacencia, puntos, punto = i, limites_espacio = limites_espacio, limites_tiempo = limites_tiempo, camino = [], unic = unic, caminos = caminos, velocidades = velocidades)
    return caminos


# Funciones de lectura de objetos
# ----------------------------------------------------------------------------------
def crear_grafo(nombre_archivo:str)->list:
    """
    Función de lectura de un archivo obj para obtener los datos necesarios

    Args:
        nombre_archivo (str): Nombre del archivo obj

    Returns:
        list: Te de devuelve una lista con los siguientes datos:
            Puntos: Lista de puntos
            Conexiones: Lista de triangulos del archivo obj
            Normales: Lista de normales
            Matriz de adyacencia: Matriz de adyacencia creada a partir de las relaciones
            Caminos: Reentradas anatómicas encontradas 
    """
    [puntos, conexiones, normales, matriz_adyacencia, triangulos] = leer_obj(nombre_archivo)
    caminos = busca_caminos(matriz_adyacencia, triangulos)
    caminos = filtrar_caminos(caminos)
    return [puntos, conexiones, normales, matriz_adyacencia, caminos]

def busca_caminos(matriz_adyacencia:list, triangulos:list)->list:
    """
    Función para buscar las reentradas anatómicas

    Args:
        matriz_adyacencia (list): Matriz de adyacencia del grafo
        triangulos (list): Triangulos del obj

    Returns:
        list: Lista de caminos que forman las reentradas anatómicas
    """
    caminos = []
    for i in range(len(matriz_adyacencia[0])):
        obtener_camino(matriz_adyacencia, triangulos, i, caminos = caminos)
    return caminos

def filtrar_caminos(caminos:list)->list:
    """
    Función para filtrar los caminos repetidos para que solo haya una vez cada reentrada anatómica

    Args:
        caminos (list): Lista de reentradas anatómicas

    Returns:
        list: Lista de reentradas anatómicas únicas
    """
    caminos_filtrados = []
    unic = []
    for camino in caminos:
        camino2 = list(camino)
        camino2 = camino2[:-1]
        camino2.sort()
        if camino2 not in unic:
            unic.append(camino2)
            caminos_filtrados.append(camino)
    return caminos_filtrados

def obtener_camino(matriz_adyacencia:list, triangulos:list, punto:int, anterior:int = -1, camino:list = [], caminos:list = [])->None:
    """
    Función para obtener los caminos de las reentradas anatómicas mirando la cantidad de triangulos que lo conectan (si es menor que dos significa que hay algún agujero)

    Args:
        matriz_adyacencia (list): Matriz de adyacencia del grafo
        triangulos (list): Lista de triangulos del obj
        punto (int): Punto en el que se encuentra
        anterior (int, optional): Punto del que viene (-1 si empieza el camino). Defaults to -1.
        camino (list, optional): Camino actual. Defaults to [].
        caminos (list, optional): Lista de caminos. Defaults to [].

    Returns:
        None
    """
    if camino == []:
        camino = [punto]
    conexiones = [i for i, x in enumerate(matriz_adyacencia[punto]) if x == 1]
    triangulos_punto = [x for x in triangulos if punto in x]

    for i in conexiones:
        if i == anterior: continue
        triangulos_punto_conexion = [x for x in triangulos_punto if i in x]
        if len(triangulos_punto_conexion) < 2:
            if i in camino:
                camino2 = list(camino)
                camino2.append(i)
                return caminos.append(camino2)
            else:
                camino2 = list(camino)
                camino2.append(i)
                obtener_camino(matriz_adyacencia, triangulos, i, punto, camino2, caminos)

def leer_obj(nombre_archivo:str)->list:
    """
    Función para leer un archivo obj y obtener los datos necesarios

    Args:
        nombre_archivo (str): Nombre del archivo obj

    Returns:
        list: Lista con los siguientes datos:
            Puntos: Lista de puntos
            Conexiones: Lista de conexiones
            Normales: Lista de normales
            Matriz de adyacencia: Matriz de adyacencia creada a partir de las relaciones
            Triangulos: Lista de triangulos del obj 
    """
    puntos = []
    conexiones = []
    normales = []
    triangulos = []
    matriz_adyacencia = None
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            if linea.startswith('v '):
                puntos.append([float(x) for x in linea.split()[1:]])
            elif linea.startswith('vn '):
                normales.append([float(x) for x in linea.split()[1:]])
            elif linea.startswith('f '):
                if matriz_adyacencia is None:
                    matriz_adyacencia = np.zeros((len(puntos), len(puntos)))
                conects = [int(x.split('/')[0]) for x in linea.split()[1:]]
                conexiones.append(conects)
                matriz_adyacencia[conects[0]-1, conects[1]-1] = 1
                matriz_adyacencia[conects[1]-1, conects[0]-1] = 1
                matriz_adyacencia[conects[1]-1, conects[2]-1] = 1
                matriz_adyacencia[conects[2]-1, conects[1]-1] = 1
                matriz_adyacencia[conects[0]-1, conects[2]-1] = 1
                matriz_adyacencia[conects[2]-1, conects[0]-1] = 1
                
                # if conects[0] == 383 or conects[1] == 383 or conects[2] == 383:
                #     print(conects)
                triangulos.append([x -1 for x in conects])
    return [puntos, conexiones, normales, matriz_adyacencia, triangulos]

def crear_matriz(puntos:list, conexiones:list)->list:
    """
    Función para crear una matriz de adyacencia a partir de los puntos y las conexiones

    Args:
        puntos (list): Lista de puntos
        conexiones (list): Lista de conexiones

    Returns:
        list: Matriz de adyacencia
    """
    matriz = []
    for i in range(len(puntos)):
        matriz.append([0]*len(puntos))
    for conexion in conexiones:
        matriz[conexion[0]-1][conexion[1]-1] = 1
        matriz[conexion[1]-1][conexion[0]-1] = 1
        matriz[conexion[1]-1][conexion[2]-1] = 1
        matriz[conexion[2]-1][conexion[1]-1] = 1
        matriz[conexion[0]-1][conexion[2]-1] = 1
        matriz[conexion[2]-1][conexion[0]-1] = 1
    return matriz


# Funciones para reentradas funcionales
# ----------------------------------------------------------------------------------
def punto_mas_cercano(punto:int, puntos:list)->int:
    """
    Función para encontrar el punto más cercano a un punto (usarla para ver el punto más cercano de cada uno de los puntos del obj al VTK)

    Args:
        punto (int): Punto de origen
        puntos (list): Lista de puntos del VTK

    Returns:
        int: Punto más cercano del obj al VTK
    """
    punto = np.array(punto)
    distancias = [np.linalg.norm(punto - p) for p in puntos]
    return np.argmin(distancias)

def crear_grafo_vtk(nombre_archivo_obj:str, nombre_archivo_csv_vtk:str)->list:
    """
    Función para crear un grafo a partir de un archivo obj y un archivo csv del VTK con los datos de cada punto

    Tarda tiempo en ejecutar, recomendable guardar en csv el datos_puntos
    Args:
        nombre_archivo_obj (str): Nombre del archivo obj con la estructura del corazón
        nombre_archivo_csv_vtk (str): Nombre del archivo csv con los datos del VTK

    Returns:
        list: Lista con los siguientes datos:
            Puntos: Lista de puntos
            Conexiones: Lista de conexiones
            Normales: Lista de normales
            Matriz de adyacencia: Matriz de adyacencia creada a partir de las relaciones
            Triangulos: Lista de triangulos del obj
            Datos de los puntos: Datos de los puntos del VTK
    """
    [puntos, conexiones, normales, matriz_adyacencia, triangulos] = leer_obj(nombre_archivo_obj)
    datos_vtk = pd.read_csv(nombre_archivo_csv_vtk)
    puntos_vtk = np.array(datos_vtk[["Points:0",  "Points:1",  "Points:2"]].values)
    datos_puntos = pd.DataFrame(puntos, columns = ["x", "y", "z"])
    datos_puntos["punto_mas_cercano"] = [punto_mas_cercano(p, puntos_vtk) for p in puntos]
    datos_puntos["material"] = datos_vtk["material"].values[datos_puntos["punto_mas_cercano"]]
    datos_puntos["model"] = datos_vtk["model"].values[datos_puntos["punto_mas_cercano"]]
    datos_puntos["f_x"] = datos_vtk["fibers:0"].values[datos_puntos["punto_mas_cercano"]]
    datos_puntos["f_y"] = datos_vtk["fibers:1"].values[datos_puntos["punto_mas_cercano"]]
    datos_puntos["f_z"] = datos_vtk["fibers:2"].values[datos_puntos["punto_mas_cercano"]]
    return [puntos, conexiones, normales, matriz_adyacencia, datos_puntos]


# Funciones de lectura guardadas de anteriores ejecuciones
# ----------------------------------------------------------------------------------
def guardar_caminos(caminos:list, nombre_archivo:str)->None:
    """
    Función para guardar los caminos en un archivo csv

    Args:
        caminos (list): Lista de caminos a guardar (Tanto reentradas anatómicas, como funcionales)
        nombre_archivo (str): Nombre del archivo csv donde se guardará
    """
    pd.DataFrame(np.array(caminos)).to_csv(nombre_archivo, index = False)
    
def cargar_caminos(nombre_archivo:str)->list:
    """
    Función para cargar los caminos guardados en un archivo csv

    Args:
        nombre_archivo (str): Nombre del archivo csv donde se guardaron los caminos

    Returns:
        list: Lista de caminos guardados
    """
    caminos = pd.read_csv(nombre_archivo)
    return [ast.literal_eval(s) for s in np.array(caminos["0"]).tolist()]
    
def obtener_velocidades_csv(nombre_archivo:str, puntos_detallados:pd.DataFrame, penalizador:float)->list:
    """
    Función para obtener las velocidades de los puntos de un archivo csv

    Args:
        nombre_archivo (str): Nombre del archivo csv con las velocidades según material y model (Tiene que llamarse velocidad y anisotropía los otros datos y velocidad tiene que ser en cm/s y anisotropía en % (se puede hacer dividiendo la velocidad en la traspuesta entre la velocidad normal))
        puntos_detallados (pd.DataFrame): DataFrame con los puntos detallados con los datos del VTK
        penalizador (float): Penalizador para las velocidades (Número entre 0 y 1)

    Returns:
        list: Lista de velocidades de los puntos
            Estructura de las velocidades: [[vector], velocidad (mm/s), penalización]
    """
    DICCIONARIO_MATERIALS = {"1": "RA", "2": "CT", "3": "PVS", "4": "BB", "5": "IST", "6": "SAN", "7": "LFO", "8": "CS", "9":"MV/LAA", "10": "FO"}
    DICCIONARIO_MODELS = {"191": "RA", "192": "BB", "193": "RAA", "194": "TV", "195": "LA", "196": "PVS", "197":"LAA", "198":"MV"}
    model_material = np.array([puntos_detallados["material"].values, puntos_detallados["model"].values]).T
    dataframe_velocidades = pd.DataFrame(np.unique(model_material, axis = 0), columns = ["material", "model"])
    
    dataframe_velocidades["material_name"] = dataframe_velocidades["material"].apply(lambda x: DICCIONARIO_MATERIALS[str(x)])
    dataframe_velocidades["model_name"] = dataframe_velocidades["model"].apply(lambda x: DICCIONARIO_MODELS[str(x)])
    
    csv_velocidades = pd.read_csv(nombre_archivo)
    
    df_vels = pd.merge(dataframe_velocidades, csv_velocidades, on=["material", "model"], how='left')
    df_vels["velocidades"] = df_vels["velocidades"] * 10 * penalizador # Pasar a mm/s y penalizadas para probar corazones menos sanos
    
    puntos_con_velocidades = pd.merge(puntos_detallados, df_vels, on=["material", "model"], how='left')
    
    velocidades_formateada = []
    for indice, fila in puntos_con_velocidades.iterrows():
        coords = [fila["f_x"], fila["f_y"], fila["f_z"]]
        dato = [coords, fila["velocidades"], fila["anisotropia"]]
        velocidades_formateada.append(dato)
    return velocidades_formateada