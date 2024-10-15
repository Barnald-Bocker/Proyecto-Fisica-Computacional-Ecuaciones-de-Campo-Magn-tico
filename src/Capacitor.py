#!/usr/bin/env python3

# Importa las bibliotecas necesarias para los diferentes cálculos y para graficar
# los resultados obtenidos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def grid_generator(factor = 1, V1 = 1.0, V2 = -1.0, boundary = 0.0):
    """
    Es la función que genera la grilla en la cual se va a trabajar para obtener la solución a la ecuación diferencial. El tamaño de la grilla va a ser de 10*factor filas y 10*factor+1 columnas. La grilla está hecha para manejar condiciones de Dirichlet. El propósito es modelar un capacitor, por eso se van a generar dos placas, una con voltaje V1 ubicada en x = 2 y la otra con voltaje V2 ubicada en x = 8. Asimismo, el borde del capacitor va a tener un voltaje en defecto igual a 0V. El problema trata con una región cuadrada de 10x10 cm, por lo tanto para obtener más precisión se usa el factor de escala ya mencionado.

    Examples:
        factor = 1
        V1 = 1.0
        V2 = -1.0
        boundary = 0.0

        >>> print(grid_generator(1.0,1.0,-1.0,0.0))
         [[ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
          [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
          [ 0.  0.  1.  0.  0.  0.  0.  0. -1.  0.  0.]
          [ 0.  0.  1.  0.  0.  0.  0.  0. -1.  0.  0.]
          [ 0.  0.  1.  0.  0.  0.  0.  0. -1.  0.  0.]
          [ 0.  0.  1.  0.  0.  0.  0.  0. -1.  0.  0.]
          [ 0.  0.  1.  0.  0.  0.  0.  0. -1.  0.  0.]
          [ 0.  0.  1.  0.  0.  0.  0.  0. -1.  0.  0.]
          [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
          [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]]

    Args:
        factor (int): Primer argumento, es el factor de escala para obtener una grilla más detallada. Tiene que ser un entero
        V1 (float): Segundo argumento, es el voltaje que corresponde a la primera placa del capacitor.
        V2 (float): Tercer argumento, es el voltaje que corresponde a la segunda placa del capacitor.
        boundary (float): Cuarto argumento, es el voltaje que corresponde a la frontera que delimita el capacitor
    Returns:
        phi (Numpy array): La grilla 10x10cm escalada por el factor para que sea un numpy array de tamaño (10*factor, 10*factor+1) para modelar el capacitor. El +1 en las columnas es debido a que una grilla que tiene una cantidad M de divisiones, contiene M+1 puntos.
    """
    # Crea el numpy array lleno de ceros para empezar a trabajar el problema.
    # Se agregar el "dtype = float" para guardar los espacios con números flotantes desde un principio.
    # factor tiene que ser un int, o si no va a tirar error
    factor = int(factor)
    phi = np.zeros((10*factor,10*factor+1), dtype = float)
    # Se incluyen las condiciones de frontera, como se escala con el factor, mantiene las dimensiones del problema original pero se tiene que multiplicar por este mismo. Se hace uso de "list slicing" para evitar el uso de for loops. Se hace desde 2*factor hasta 8*factor para no afectar las condiciones de frontera, pero aún así incluir los 6*factor puntos correspondientes a la grilla. Se aplica en las columnas correspondientes a 2*factor y 8*factor.
    phi[2*factor:8*factor, 2*factor] = V1
    phi[2*factor:8*factor, 8*factor] = V2
    # Se devuelve la matriz ya construida
    return phi

print(grid_generator(1.1,1.0,-1.0,0.0))
