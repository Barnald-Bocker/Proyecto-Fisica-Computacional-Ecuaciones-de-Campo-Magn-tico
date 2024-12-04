#Métodos numéricos usados en el proyecto 

##Método de relajación de Jacobi 

Este es un método que utiliza diferencias finitas para definir las derivadas del siguiente problema el cual corresponde a la ley de Gauss en dos dimensiones:

\begin{equation}
   \frac{\partial^{2}\phi}{\partial{x^{2}}} + \frac{\partial^{2}\phi}{\partial{y^{2}}} = 0
\end{equation}

en donde las derivadas vistas se pueden aproximar usando diferencias centrales para así poder expresarlas de la siguiente manera

\begin{equation}
 \frac{\partial^{2}\phi}{\partial{x^{2}}} \approx \frac{\phi (x +a, y) - 2\phi (x,y) + \phi (x -a,y)}{a^{2}}
\end{equation}

\begin{equation}
 \frac{\partial^{2}\phi}{\partial{y^{2}}} \approx \frac{\phi (x,y + a) - 2\phi (x,y) + \phi (x,y - a)}{a^{2}}
\end{equation}

De esta manera, la ecuación diferencial dada se puede aproximar de la siguiente manera:
\begin{equation}
 \phi (x +a, y) - 2\phi (x,y) + \phi (x -a,y) + \phi (x,y + a) - 2\phi (x,y) + \phi (x,y - a) \approx  0
\end{equation}

\begin{equation}
 \phi (x +a, y) + \phi (x -a,y) + \phi (x,y + a) + \phi (x,y - a) -4\phi (x,y) \approx 0
\end{equation}

De esta manera se consigue una nueva forma de la ley de Gauss.
\begin{equation}
 \phi (x,y) = \frac{1}{4} \cdot \phi (x +a, y) + \phi (x -a,y) + \phi (x,y + a) + \phi (x,y - a)
\end{equation}
Como se puede observar, cada punto $\phi(x,y)$ se calcula al usar los valores próximos a dicho punto.
De esta manera se encuentran una serie de ecuaciones que corresponden a un sistema de ecuaciones a resolver. De este mismo sistema se puede extraer una matriz, a este arreglo se le va a nombrar la matriz A con entradas $a_{ij}$ dentro de esta matriz recae la condicion de convergencia de este método la conforman los criterios de la diagonal dominante, los cuales son los siguiente

condición necesaria: $|a_{ii}| > |a_{ij}|$

-condición suficiente: $|a_{ii}| > \sum|a_{ij}|$

La condición necesaria se debe cumplir obligatoriamente para todos los elementos de la diagonal de la matriz.
La condición suficiente, una vez cumplida, se puede asegurar que se llegara a la convergencia. Puede haber casos que no cumplan esta condicion y aun asi se puede dar la convergencia
Es importante aclarar que este método solo proporciona una aproximación a la solución, esta aproximación puede ser precisa dependiendo del valor de tolerancia que se logre alcanzar.
##Metodo de Jacobi modificado(overrelaxation)
El metodo de sobrerelajación de Jacobi forma parte de los métodos de relajación acelerada en el cual se logra acelerar la convergencia gracias a un parámetro $\omega$, dando lugar a una nueva forma de definir $\phi_{prime}$ la cual es la siguiente:
\begin{equation}
 \phi_{prime} = (1+\omega)(\phi (x +a, y) + \phi (x -a,y) + \phi (x,y + a) + \phi (x,y - a)) - \omega \phi (x,y)
\end{equation}
A nivel de matemáticas se tiene que considerar el siguiente vector el cual se va a llamar el vector residuo $r = b- A\bar{x}$ siendo $\bar{x}$ una aproximación a la solución del problema a resolver de forma $Ax = b$. Para cada iteración de procesos como Jacobi se da una aproximación a la solución por cada iteración y por lo tanto un vector residual. El objetivo es que este vector residual converga a 0 de manera rapida al usar el parámetro $\omega$. 
De este parámetro dependera la convergencia del metodo. $\omega$ se escoge dependiendo del problema a resolver. 

si $0 < \omega < 1$ se le conoce como under-relaxation.
 
si $1 < \omega$ se le conoce como over-relaxation.

Si se hace buen uso de este método, la convergencia se logra de buena manera, sin embargo este método no siempre logra ser convergente, esto puede ser tanto por la escogencia del parámetro $\omega$ o por la naturaleza del problema en si.    

##Metodo de Gauss-Eidel

Este método es una forma optimizada para el metodo de relajación de Jacobi. Si se aplica correctamente, se requieren menos iteraciones para que se llegue a la convergencia. A nivel computacional lo que sucede es que ya no se necesita solicitar un nuevo arreglo de datos $\phi_{prime}$, los datos se iran actualizando dentro del mismo $\phi$. De esta manera no se tiene que solicitar mas memoria.
Los criterios de convergencia son los mismos que los del método de relajacion de Jacobi.
A nivel matemático, el metodo de jacobi asume un vector solución  inicial $x_{0} = 0$ y se consiguen las nuevas actualizaciones del vector. Mientras que Gauss Eidel no se espera hasta conseguir el vector completo, los valores se van actualizando conforme se consiguen, sin esperar a que se consiga el vector solución completo correspondiente. De esta manera se consigue una convergencia acelerada y mayor eficiencia computacional.
