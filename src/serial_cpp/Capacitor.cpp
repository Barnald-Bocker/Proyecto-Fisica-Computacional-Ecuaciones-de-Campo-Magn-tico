#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> grid_generator(int factor, double V1, double V2, double boundary = 0.0){

  // Manejo de excepciones para evitar factores de escala no válidos.
  if(factor < 1){
    std::cerr << "El factor de escala debe ser un número entero mayor que cero!" << std::endl;
    exit(1);
    }

  // Se define la variable 'column' para facilitar la visualización del código, al ser un valor sumamente utilizado durante el proceso.
  int column = 10*factor+1;
  // Se genera una grilla de dimensiones (10*factor)x(10*factor+1), empleando un vector para que sea contiguo en memoria.
  std::vector<double> phi(10*factor*column);
  // Generar las condiciones de frotenra, correspondientes a las líneas verticales de potencial. 
  for(int i=2*factor; i < 8*factor; ++i){
    phi[i*column+2*factor] = V1;
    phi[i*column+8*factor] = V2;
  }
  // Definimos como valor de retorno la grilla.
  return phi;
}

int gauss_seidel(std::vector<double>& phi, double omega, int factor, double tolerance){

  // Manejo de excepciones para asegurar que el factor de escala ingresado a la función sea el mismo valor utilizado en la generación de la grilla.
  if(phi.size() != (10*factor)*(10*factor+1)){
    std::cerr << "El factor de escala ingresado es diferente al utilizado para generar la grilla!" << std::endl;
    exit(1);
  }

  // Manejo de excepciones para asegurar un valor de tolerancia válido.
  if(tolerance <= 0){
    std::cerr << "El valor de tolerancia no puede ser nulo ni negativo!" << std::endl;
    exit(1);
    }

  // Definimos una copia de la grilla, servirá para establecer la diferencia de la grilla después de cada iteración.
  std::vector<double> phiCopy = phi;
  // Definimos una variable que irá almacenando la diferencia máxima entre grillas, con un valor inicial de '1.0'.
  double delta = 1.0;
  // La variable 'diff' se empleará para ir calculando la diferencia entrada por entrada entre las grillas. Si el valor de 'diff' calculado en una iteracion es mayor al almacenado en 'maxDelta', se reestablece el valor de maxDelta.  
  double diff;
  double maxDelta;
  // La variable 'its' hara la funcion de contador del numero de iteraciones.
  int its = 0;
  // Nuevamente, se define la variable 'column' para simplificar la lectura del codigo.
  int column = 10*factor+1;

  // Ciclo encargado de iterar hasta haber alcanzado la tolerancia.
  while(delta > tolerance){
  // Ciclo 'for' que corresponde a las filas de la grilla.
  for(int i=1; i < 10*factor-1; ++i){
    // Ciclo 'for' que corresponde a las columnas de la grilla.
    for(int j=1; j < column-1; ++j){
      // Condiciones de frontera. En caso de encontrarse con las lineas de potencial, no se calcula el promedio sobre la entrada.
      if((2*factor <= i && i < 8*factor) && (j == 2*factor || j == 8*factor)){
        continue;
	// En caso de analizar una entrada diferente a las condiciones de frontera, se calcula el promedio de las cuatro entradas subyacentes.
      } else{
	phi[i*column + j] = ((1.0+omega)*0.25)*(phi[(i+1)*column + j] + phi[(i-1)*column + j] + phi[i*column + (j+1)] + phi[i*column + (j-1)]) - omega*phi[i*column + j];
      }
    }
  }

  // Se establece un valor inicial de 0.0 para la maxima diferencia entre grillas.  
  maxDelta = 0.0;
  // Se calcula el valor absoluto (std::abs) de la diferencia entre grillas entrada por entrada, almacenando en 'maxDelta' el mayor valor encontrado.
  for(int i=1; i < 10*factor-1; ++i){
    for(int j=1; j < column-1; ++j){
      diff = std::abs(phi[i*column + j] - phiCopy[i*column + j]);
      // Si se encuentra una diferencia entre entradas mayor a la maxima diferencia, se actualiza 'maxDelta'.
      if(diff > maxDelta){
        maxDelta = diff;
      }
    }
  }

  // Se asigna la maxima diferencia encontrada a la variable 'delta'. En caso de ser mayor a la tolerancia, sigue iterando. Caso contrario, se detiene la condicion.
  delta = maxDelta;
  // Establecer una copia de la grilla.
  phiCopy = phi;
  // Sumar una iteracion.
  its += 1;
  }
  std::cout << "Diferencia máxima: " << maxDelta << std::endl;
  return its;
}

void print_grid(std::vector<double> &matrix, int factor){
  std::cout.precision(4);
  for (int i = 0; i < 10*factor; i++){
    for (int j = 0; j < 10*factor+1; j++){
      if (j<10*factor){
        std::cout<< matrix[i*(10*factor+1)+j] << ";";
      }
      else{
        std::cout<<matrix[i*(10*factor+1)+j]<< std::endl;
      }
      }
      }
 }

int main(){
  int factor = 10;
  double volt1 = 1.0;
  double volt2 = -1.0;
  std::vector<double> grid = grid_generator(factor, volt1, volt2);
  double omega = 0.9;
  double tolerance = 1e-4;
  int iterations = gauss_seidel(grid,omega, factor, tolerance);
  print_grid(grid,factor);
}
