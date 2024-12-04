#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

void gauss_seidel(double V1, double V2, double omega, int factor, double tolerance){

  if(tolerance <= 0){
    std::cerr << "El valor de tolerancia no puede ser nulo ni negativo!" << std::endl;
    exit(1);
  }

  double delta = 1.0;
  double diff;
  double local_maxDelta;
  double new_value;
  int its = 0;
  int column = 10*factor+1;

  int size, rank;
  int row = 10*factor;
  double time_1 = MPI_Wtime();
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nlocal = row/size;
  int rest = row%size;

  if (rest && rank<rest){
    nlocal ++;
  }

  int start = nlocal*rank;
  if (rest && rest <= rank){
    start += rest;
  }

  int end;

  // En este caso, 'end' permite delimitar en qué punto de la 'grilla local' (local_section) se debe detener la iteración, para evitar tanto las condiciones de frontera como las filas fantasma.
  if(rank == 0 || rank == size-1){
    end = nlocal;
  } else{
    end = nlocal+1;
  }

  if(size == 1){
    end = nlocal-1;
  }

  std::vector<double> local_section;

  // IMPORTANTE: Celdas fantasma.
  // Se asigna una fila fantasma al primer y último proceso.
  // Se asginan dos filas fantasma (superior e inferior) a los procesos intermedios.
  if(size > 1){
    if(rank == 0 || rank == size - 1){
      local_section.resize((nlocal+1)*column, 0.0); // MUY IMPORTANTE: GRILLAS FANTASMA
    } else{
      local_section.resize((nlocal+2)*column, 0.0);
    }
  } else{
    local_section.resize(nlocal*column, 0.0);
  }

  // Establecer las condiciones de frontera.
  // Es necesario considerar la posición absoluta de cada entrada, por lo que es necesario un ajuste para contrarrestar el desplazamiento de las filas por la incorporación de filas fantasma.
  for(int i=1; i < end; ++i){
    int absolute_position = start+i-1;
    if(rank == 0){
      absolute_position = start+i;
    }
    for(int j=1; j < column-1; ++j){
      if ((2*factor <= absolute_position && absolute_position < 8*factor) && (j == 2*factor)){
        local_section[i*column + j] = V1;
      }
      if ((2*factor <= absolute_position && absolute_position < 8*factor) && (j== 8*factor)){
        local_section[i*column + j] = V2;
      }
    }
  }
   
  // En caso de evaluarse el código con un sólo proceso, no se necesita intercomunicación para filas fantasma.
  // Por lo tanto, esta sección se ejecutará únicamente para dos o más procesos.
  if(size > 1){
    // Se establece un 'request' por cada tipo de proceso.
    // Envío [0] y recepción [1] de fila superior.
    // Envío [2] y recepción [3] de fila inferior.
    MPI_Request request[4];
    // Transferencia de las filas fantasma superiores.
    if(rank > 0){
      // Todos los procesos, exceptuando el primero, envían su segunda fila (omitiendo la fila fantasma) al proceso anterior.
      // A su vez, reciben una fila del proceso anterior, almacenándola en la fila fantasma superior. 
      MPI_Isend(&local_section[column], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[0]);
      MPI_Irecv(&local_section[0], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    }

    // Transferencia de las filas fantasma inferiores.
    if(rank == 0){
      // El primer proceso envía la penúltima fila al siguiente proceso.
      // A su vez, recibe una fila del siguiente proceso, almacenándola en la fila fantasma.
      MPI_Isend(&local_section[(nlocal-1)*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
      MPI_Irecv(&local_section[nlocal*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
    } else if(rank < size - 1 && rank > 0){
      // Las filas intermedias repiten el proceso anterior, realizando un ajuste para 'nlocal'.
      MPI_Isend(&local_section[nlocal*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
      MPI_Irecv(&local_section[(nlocal+1)*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    }

    // Esperar a que se completen las transferencias de las filas fantasma superiores.
    if(rank > 0){
      // Envío de fila superior.
      MPI_Wait(&request[0], MPI_STATUS_IGNORE);
      // Recepción de fila superior.
      MPI_Wait(&request[1], MPI_STATUS_IGNORE);
    }

    if(rank < size - 1){
      // Envío de fila inferior.
      MPI_Wait(&request[2], MPI_STATUS_IGNORE);
      // Recepción de fila inferior.
      MPI_Wait(&request[3], MPI_STATUS_IGNORE);
    }
  // Finaliza la intercomunicación de procesos.
  }

  while(delta > tolerance){

    for(int i=1; i < end; ++i){
      int absolute_position = start+i-1;
      if (rank == 0){
        absolute_position = start+i;
      }

      // Establecemos un valor inicial para la mayor diferencia.
      local_maxDelta = 0.0;

      // Iteraciones del método de Gauss-Seidel.
      for(int j=1; j < column-1; ++j){
        if((2*factor <= absolute_position && absolute_position < 8*factor) && (j == 2*factor || j == 8*factor)){
          continue;
        } else{
	  // Cálculo del nuevo valor que recibirá la celda.
	  // Se define 'new_value' para almacenar temporalmente el valor calculado con la iteración de Gauss-Seidel.
	  // Esta estrategia permite desarrollar la metodología de Gauss-Seidel sin necesitar una copia de la grilla local.
	  new_value = ((1.0+omega)*0.25)*(local_section[(i+1)*column + j] + local_section[(i-1)*column + j] + local_section[i*column + (j+1)] + local_section[i*column + (j-1)]) - omega*local_section[i*column + j];
	  // Se compara la diferencia entre el nuevo valor calculado y el valor de la celda con el 'local_maxDelta'.
	  // En caso de que la diferencia sea mayor a la almacenada en 'local_maxDelta', se actualiza el valor máximo.
	  local_maxDelta = std::max(local_maxDelta, std::abs(local_section[i*column + j] - new_value));
	  // Finalmente, se asigna el valor calculado a la entrada correspondiente.
	  local_section[i*column + j] = new_value;
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Sumamos una iteración al contador.
    its++;
    // Volvemos a intercomunicar las filas fantasma, con el mismo proceso anterior al ciclo.
    if(size > 1){
      MPI_Request request[4];
      if(rank > 0){
        MPI_Isend(&local_section[column], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[0]);
        MPI_Irecv(&local_section[0], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
      }

      if(rank == 0){
        MPI_Isend(&local_section[(nlocal-1)*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
        MPI_Irecv(&local_section[nlocal*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
      } else if(rank < size - 1 && rank > 0){
        MPI_Isend(&local_section[nlocal*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
        MPI_Irecv(&local_section[(nlocal+1)*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
      }

      if(rank > 0){
        MPI_Wait(&request[0], MPI_STATUS_IGNORE);
        MPI_Wait(&request[1], MPI_STATUS_IGNORE);
      }

      if(rank < size - 1){
        MPI_Wait(&request[2], MPI_STATUS_IGNORE);
        MPI_Wait(&request[3], MPI_STATUS_IGNORE);
      }
    }

    // Comunicar el delta.
    double global_delta = 0.0;
    MPI_Allreduce(&local_maxDelta, &global_delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    delta = global_delta;

    if (delta > 1e10) {
      std::cerr << "WARNING: Delta is increasing uncontrollably, consider changing omega value." << std::endl;
      break;
    }

  
  }
  // Una vez finalizadas las iteraciones, se detiene el conteo de tiempo.
  // Considerar el tiempo que se dure imprimiendo sesgaría la verdadera evaluación de la escalabilidad de las iteraciones.
  double time_2 = MPI_Wtime();

  // Imprimir la matriz resultante.
  MPI_Request print_request;
  int signal = 1;

  // Si el código se ejecuta con únicamente un proceso, simplemente se imprime la matriz.
  if(size == 1){ 
    std::cout << "Número de procesos: " << size << std::endl;
    std::cout << "Tiémpo transcurrido: " << time_2 - time_1 << std::endl;
    for(int i = 0; i < row; i++){
      for(int j = 0; j < column; j++){
	if(j < column - 1){
        std::cout << local_section[i*column+j] << ";";
	} else{
	  std::cout << local_section[i*column+j] << std::endl;
	}
      }
    }
  // En caso de que se empleen dos o más procesos, es necesario coordinar el orden de impresión de los procesos.
  // Podemos emplear una variable 'signal' que se encargue de notificar a un proceso cuando el proceso anterior ha terminado de imprimir.
  } else{
    // El proceso 0 imprime de primero, al ser quien posee las primeras entradas de la grilla.
    // Es importante considerar que, para el proceso 0, local_grid contiene una fila fantasma al final, por lo que debemos omitirla.
    if(rank == 0){
      std::cout << "Número de procesos: " << size << std::endl;
      std::cout << "Tiempo transcurrido: " << time_2 - time_1 << std::endl;

      for(int i = 0; i < end; i++){
        for(int j = 0; j < column; j++){
	  if(j < column - 1){
	    std::cout << local_section[i*column+j] << ";";
	  } else{
	    std::cout << local_section[i*column+j] << std::endl;
	  }
        }
      }

      // Una vez el primer proceso finaliza la impresión, notifica al segundo proceso mediante el envío de la variable 'signal'.
      MPI_Isend(&signal, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &print_request);
      MPI_Wait(&print_request, MPI_STATUS_IGNORE);
    } else{
      // Esta sección del código la ejecutarán todos los procesos a partir del proceso 1.
      // Primeramente, se recibe la señal del proceso anterior, lo que inicializará el proceso una vez recibido 'signal'.
      MPI_Irecv(&signal, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &print_request);
      MPI_Wait(&print_request, MPI_STATUS_IGNORE);

      // Los procesos 'intermedios' (0 < rank < size - 1) contienen dos filas fantasmas, superior e inferior.
      // Por lo tanto, debemos omitir la primera y la última fila de local_section al imprimir.
      if(rank < size-1){
        for(int i = 1; i < end; i++){
          for(int j = 0; j < column; j++){
            if(j < column - 1){
              std::cout << local_section[i*column+j] << ";";
            } else{
              std::cout << local_section[i*column+j] << std::endl;
            }
          }
        }
	// Una vez se termina la impresión en terminal, se notifica al siguiente proceso.
        MPI_Isend(&signal, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &print_request);
        MPI_Wait(&print_request, MPI_STATUS_IGNORE);
      // Nótese que la siguiente sección del código únicamente será accesada por el último proceso.
      // En este caso, local_grid sólo contiene una fila fantasma superior, por lo que omitimos dicha fila al imprimir.
      // Además, mediante el 'if' establecido para los procesos intermedios, evitamos que el último proceso trate de enviar una señal al siguiente proceso.
      } else{
        for(int i = 1; i < end+1; i++){
          for(int j = 0; j < column; j++){
	    if(j < column-1){
	      std::cout << local_section[i*column+j] << ";";
	    } else{
	      std::cout << local_section[i*column+j] << std::endl;
	    }
	  }
        }
      }
    }
  }

  // Establecemos una barrera para asegurar que todos los procesos hayan finalizado antes de cerrar el ambiente MPI.
  MPI_Barrier(MPI_COMM_WORLD);
  // Finalización del ambiente MPI.
  MPI_Finalize();
}

int main(){
  int factor = 10;
  double volt1 = 1.0;
  double volt2 = -1.0;
  double omega = 0.8;
  double tolerance = 1e-4;
  gauss_seidel(volt1, volt2, omega, factor, tolerance);
  return 0;
}
