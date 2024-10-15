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

def gauss_seidel(phi, omega, factor, tol):
    """
    Esta es la función que resuelve la ecuación diferencial parcial de Laplace de manera iterativa. En este caso la solución se guarda y se trabaja sobre la misma matriz, lo que es la principal diferencia de los métodos de Jacobi (modificado y simple). Una propiedad que es sumamente importante, es que siempre es estable.

    Examples:
    factor = 1
    V1 = 1.0
    V2 = -1.0
    omega = 0.9
    tol = 1e-5

    Examples:
         >>> result, iterations = gauss_seidel( 1, 0.9, 1, 1e-5)
         >>> print(result)
         [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
          [ 0.00000000e+00  1.94984598e-01  3.61769011e-01  2.52087439e-01 1.21533562e-01 -1.41308295e-06 -1.21531760e-01 -2.52085676e-01 -3.61767842e-01 -1.94984609e-01  0.00000000e+00]
          [ 0.00000000e+00  4.18168790e-01  1.00000000e+00  5.25040044e-01  2.34042744e-01  1.03625792e-06 -2.34042966e-01 -5.25038493e-01 -1.00000000e+00 -4.18167388e-01  0.00000000e+00]
          [ 0.00000000e+00  4.77685703e-01  1.00000000e+00  6.14026208e-01  2.89605985e-01 -5.83595094e-08 -2.89606016e-01 -6.14026132e-01 -1.00000000e+00 -4.77681031e-01  0.00000000e+00]
          [ 0.00000000e+00  4.92561212e-01  1.00000000e+00  6.41460677e-01  3.10354450e-01  1.52099287e-06 -3.10355266e-01 -6.41457315e-01 -1.00000000e+00 -4.92557557e-01  0.00000000e+00]
          [ 0.00000000e+00  4.92563136e-01  1.00000000e+00  6.41461553e-01  3.10355031e-01  3.98493823e-07 -3.10355102e-01 -6.41460292e-01 -1.00000000e+00 -4.92560068e-01  0.00000000e+00]
          [ 0.00000000e+00  4.77682608e-01  1.00000000e+00  6.14024658e-01  2.89606329e-01  2.47307453e-06 -2.89603028e-01 -6.14023807e-01 -1.00000000e+00 -4.77680087e-01  0.00000000e+00]
          [ 0.00000000e+00  4.18167581e-01  1.00000000e+00  5.25039348e-01  2.34044672e-01  5.26638008e-06 -2.34039885e-01 -5.25034828e-01 -1.00000000e+00 -4.18164196e-01  0.00000000e+00]
          [ 0.00000000e+00  1.94982205e-01  3.61767506e-01  2.52082957e-01  1.21534324e-01  3.15599565e-06 -1.21528280e-01 -2.52083071e-01 -3.61768031e-01 -1.94982326e-01  0.00000000e+00]
          [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]]

    Args:
        phi (Numpy array): Primer argumento, la grilla 10x10cm escalada por el factor para que sea un numpy array de tamaño (10*factor, 10*factor+1) para modelar el capacitor.
        omega (float): Segundo argumento, parámetro que depende del problema, dependiendo de su escogencia puede acelerar el cálculo.
        factor (int): Tercer argumento, es el factor de escala para obtener una grilla más detallada. Tiene que ser un entero.
        tol (float): Cuarto argumento, es la tolerancia aceptada para la precisión del problema, se obtiene obteniendo la diferencia máxima entre phi y phi_copy

    Returns:
        phi (Numpy array): La grilla 10x10cm escalada por el factor para que sea un numpy array de tamaño (10*factor, 10*factor+1) para modelar el capacitor.
        its (int): Iteraciones que fueron necesarias para llegar a la tolerancia especificada.
    """
    # Se asegura que factor sea un int
    factor = int(factor)
    # Se va a guardar una copia para evaluar el error
    phi_copy = phi.copy()
    # Desarrollo del método Gauss-Seidel
    # Se inicializa la variable que se usa para el condicional
    delta = 1.0
    # Se inicializa el contador de iteraciones
    its = 0
    # La escogencia de los parámetros de este for loop tienen un razonamiento. Se puede notar que la matriz especificada tiene tamaño (10*factor, 10*factor+1). Al empezar el for loop en 1 y terminarlo en el tamaño-1 (10*factor-1, 10*factor respectivamente) se evitan las condiciones de frontera que se encuentran en los límites de estas matrices, que serían la primer y última fila y columna, de esta forma se evita el uso de un condicional dentro del ciclo
    while delta > tol:
        for i in range(1,10*factor-1):
            for j in range(1,10*factor):
                # Condicional que evita las condiciones de frontera
                if 2*factor <= i < 8*factor and (j == 2*factor or j == 8*factor):
                    pass
                else:
                    # Para cualquier otro caso que no sean las condiciones de frontera, se evalua el algoritmo
                    phi[i,j] = ((1.0+omega)*0.25)*(phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1]) - omega*phi[i,j]
        # Se calcula la diferencia máxima con respecto a los valores que se obtuvieron en la iteración anterior
        delta = np.max(np.abs(phi - phi_copy))
        # Se copian los valores de la nueva iteración, para seguir aplicando el delta y mejorar la precisión
        phi_copy = phi.copy()
    # En el momento que se cumple la condición en la que delta <= tol, se rompe el ciclo y se devuelve la matriz con la solución phi y las iteraciones que fueron requeridas para llegar a la tolerancia mencionada
    # Por naturaleza de python, se devuelve en forma de tupla
    return phi, its

factor = 10
omega = 0.9
tol = 1e-5
main_grid = grid_generator(factor)
result, iterations = gauss_seidel(main_grid, omega, factor, tol)
plt.imshow(result, cmap = "gray")
plt.tight_layout()
plt.title(f"Gráfico para omega = {0.9}")
plt.show()
