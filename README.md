# Explicación de métodos de paralelización: Memoria compartida y Distribuida
## Memoria Compartida: Método de Wavelength

Este método se ve utilizado para calculos reiterativos, para el caso del capacitor simple se comienza en el primer elemento de la grilla, ubicado en (0,0), y se procesan de diagonal en diagonal, en este caso se utliza una diagonal creciente la cual se encuentra fijada por el punto que pertenece a la columna 0; este método al tratar con datos de las diagonales asegura una no interdependencia.

Este método presenta multiples ventajas, al ser datos independiente se da una paralelización altamente eficiente, asimismo esta metodología presenta una execlente escalabilidad donde es altamente eficiente para grillas grandes y multiples procesadores. Asimismo presenta ciertas desventajas como el hecho de que es aplicable siempre y cuando sea organizable en diagonales, de igual manera puede presentar carga desbalanceada, donde segun el número de hilos que se esten usando en las primeras iteraciones al tener pocos datos puede tener un subaprovechamiento.

## Memoría Distribuida: Método de Celdas Fantasmas

Este método se utiliza cuando al dividir un grilla en diferentes secciones asignadas cada una a un proceso diferente los datos presentan interdependencia a los valores de un proceso vecino, de esta manera surge la necesidad de una comunicación entre ellos y esto se resuelve añadiendo celdas que almacenaran los datos necesarios de otro proceso. Dicha metodología presenta una serie de ventajas como lo son la reducción de la reducción de la comunicación  entre procesos, de igual manera presenta una excelente escalabilidad conforme se añaden nodos lo que hace que sea estable.

Por otro lado tambien se presentan desventajas tales como que las celdas fantasmas aumentan la cantidad de almacenamiento utlizado debido a que se necesitan datos adicionales, asimimsmo su implementacion puede implicar un nivel de dificultad alto, asimismo la sincronización de datos puede tornarse problemática.
