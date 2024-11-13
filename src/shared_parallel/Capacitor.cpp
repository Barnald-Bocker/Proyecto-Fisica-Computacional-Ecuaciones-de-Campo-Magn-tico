#include <iostream>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

// Función generadora de la grilla inicial.
std::vector<double> grid_generator(int factor, double V1, double V2, double boundary = 0.0){

	// Manejo de excepciones para evitar factores de escala no válidos.
	if(factor < 1){
	    std::cerr << "El factor de escala debe ser un número entero mayor que cero!" << std::endl;
	    exit(1);
	}

	// Generar la grilla con las dimensiones (10*factor x 10*factor+1) y las líneas de potencial con V1, V2
	int column = 10*factor+1;
	std::vector<double> phi(10*factor*column); 
	for(int i=2*factor; i < 8*factor; ++i){
		phi[i*column+2*factor] = V1;
		phi[i*column+8*factor] = V2;
	}

	return phi;
}


// Función que regresa la medición del tiempo actual.
double seconds(){
	struct timeval tmp;
	double sec;
	gettimeofday( &tmp, (struct timezone *)0 );
	sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;

	return sec;

}

// Función que se encarga de las iteraciones del método de Gauss-Seidel.
// Como variables de entrada se encuentra una referencia a la grilla, el valor 'omega', el factor de escala y la tolerancia con la que se desea evaluar.
int gauss_seidel(std::vector<double>& phi, double omega, int factor, double tolerance){

	// Manejo de excepciones, evitar que el factor de escala ingresado en este método sea diferente al utilizado al momento de generar la grilla.
	if(phi.size() != (10*factor)*(10*factor+1)){
		std::cerr << "El factor de escala ingresado es diferente al utilizado para generar la grilla!" << std::endl;
		exit(1);
	}

	// Manejo de excepciones para asegurar que el valor de tolerancia sea un número válido.
	if(tolerance <= 0){
		std::cerr << "El valor de tolerancia no puede ser nulo ni negativo!" << std::endl;
		exit(1);
	}

	// Se realiza un deep copy de la grilla. Nótese que la sintaxis es de una copia superficial, sin embargo, la biblioteca 'vector' se encarga de gestionar esta operación como una copia profunda.
	std::vector<double> phiCopy = phi;
	double delta = 1.0;
	double diff;
	double maxDelta;
	int its = 0;
	// Se define la variable 'column' para simplificar la apariencia del código.
	int column = 10*factor+1;
	int initValue;

	// Se inicia el contador de tiempo antes de iniciar las iteraciones para poder medir la escalabilidad de la paralelización.
	double time_1 = seconds();
	int num_procs;

	// Vamos a recorrer la grilla por líneas diagonales, empezando por un primer ciclo encargado de ir variando las diagonales de tal forma que las entradas iniciales de cada diagonal se vayan desplazando a lo largo de la primera columna. Una vez se alcance la última fila de la primera columna, se realizará un desplazamiento de la entrada inicial de cada diagonal a lo largo de toda la fila mencionada, que corresponde al segundo ciclo.
	while(delta > tolerance){
		
		// Se establece 'initValue = 1' para el primer ciclo. Corresponderá a la columna que posee las entradas iniciales de las diagonales.
		initValue = 1;
		// Primer ciclo. initValue se encargará de iniciar las diagonales en la segunda columna de la grilla. Se evita la primera columna por las condiciones de frontera.
		// Inicializamos el ambiente OMP.
		#pragma omp parallel
		{
		num_procs = omp_get_num_threads();
		
		for(int i=1; i < 10*factor-1; ++i){
			// Se paraleliza el recorrido para cada línea diagonal. (Fork_1)
			#pragma omp for
			for(int d=0; d < i; ++d){ // 'd' corresponde al ciclo que recorre las líneas diagonales.
			// Nótese que ahora las entradas de la grilla corresponden a '[i-d,j+d]'.
			// Para el primer ciclo, j->initValue. Por lo tanto: '[i-d,initValue+d]'.
			// El valor 'd' se resta a la fila y suma a la columna porque el objetivo es generar líneas diagonales que se recorran con monotonía creciente.
				if((2*factor <= (i-d) && (i-d) < 8*factor) && ((initValue+d) == 2*factor || (initValue+d) == 8*factor)){
					continue;
				} else{
					phi[(i-d)*column + (initValue+d)] = ((1.0+omega)*0.25)*(phi[(i-d+1)*column + (initValue+d)] + phi[(i-d-1)*column + (initValue+d)] + phi[(i-d)*column + (initValue+d+1)] + phi[(i-d)*column + (initValue+d-1)]) - omega*phi[(i-d)*column + (initValue+d)];
				}
			}
		}
		} // Se cierra el ambiente OMP para el primer ciclo. (Join_1)

		// Segundo ciclo, initValue ahora se encargará de iniciar las diagonales en la penúltima fila de la grilla. Se evita la última fila por las condiciones de frontera.
		// Se establece 'initValue = 10*factor-2' para el segundo ciclo. Corresponderá a la fila que poseerá las entradas iniciales de las diagonales.
		initValue = 10*factor-2;
		// Se inicializa el ambiente OMP para el segundo ciclo. (Fork_2)
		#pragma omp parallel
		{
		for(int j=2; j < 10*factor; ++j){
			// Se paraleliza el recorrido para cada línea diagonal.
			#pragma omp for
			for(int d=0; d < 10*factor-j; ++d){
			// Nótese que, nuevamente, las entradas de la grilla corresponden a '[i-d,j+d]'.
			// Para el segundo ciclo, i->initValue. Por lo tanto: '[initValue-d,j+d]'.
				if((2*factor <= (initValue-d) && (initValue-d) < 8*factor) && ((j+d) == 2*factor || (j+d) == 8*factor)){
                                        continue;
                                } else{
                                        phi[(initValue-d)*column + (j+d)] = ((1.0+omega)*0.25)*(phi[(initValue-d+1)*column + (j+d)] + phi[(initValue-d-1)*column + (j+d)] + phi[(initValue-d)*column + (j+d+1)] + phi[(initValue-d)*column + (j+d-1)]) - omega*phi[(initValue-d)*column + (j+d)];
				}
			}
		}
		} // Se cierra el ambiente OMP para el segundo ciclo. (Join_2)
		
		maxDelta = 0.0;
		for(int i=1; i < 10*factor-1; ++i){
			for(int j=1; j < column-1; ++j){
				diff = std::abs(phi[i*column + j] - phiCopy[i*column + j]);
				if(diff > maxDelta){
					maxDelta = diff;
				}
			}
		}

		delta = maxDelta;
		phiCopy = phi;
		its += 1;
	}

	double time_2 = seconds();
	std::cout << "Número de procesos: " << num_procs << std::endl;
	std::cout << "Tiempo transcurrido: " << time_2 - time_1 << std::endl;
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
  int iterations = gauss_seidel(grid, omega, factor, tolerance);
  print_grid(grid, factor);
}
