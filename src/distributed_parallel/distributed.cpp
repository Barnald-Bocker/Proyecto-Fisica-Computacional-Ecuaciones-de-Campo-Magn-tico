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

  //int end = start +  nlocal;
  if(rank == 0 || rank == size-1){
    end = nlocal;
  } else{
    end = nlocal+1;
  }

  std::vector<double> local_section;

  if(rank == 0 || rank == size - 1){
    local_section.resize((nlocal+1)*column, 0.0); // MUY IMPORTANTE: GRILLAS FANTASMA
  } else{
    local_section.resize((nlocal+2)*column, 0.0);
  }

  //Establecer condiciones de frontera

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

  std::vector<double> local_sectionCopy = local_section;

  while(delta > tolerance){
    // Se establece un 'request' por cada tipo de proceso.
    // Envío [0] y recepción [1] de fila superior.
    // Envío [2] y recepción [3] de fila inferior.
    MPI_Request request[4];
    // Transferencia de las filas fantasma superiores.
    if(rank > 0){
      MPI_Isend(&local_section[column], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[0]);
      MPI_Irecv(&local_section[0], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    }

    // Transferencia de las filas fantasma inferiores.
    if(rank == 0){
      MPI_Isend(&local_section[(nlocal-1)*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
      MPI_Irecv(&local_section[nlocal*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
    } else if(rank < size - 1 && rank > 0){
      MPI_Isend(&local_section[nlocal*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
      MPI_Irecv(&local_section[(nlocal+1)*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
    }

    // Esperar a que se completen las comunicaciones entre procesos.
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

    for(int i=1; i < end; ++i){
      int absolute_position = start+i-1;
      if (rank == 0){
        absolute_position = start + i;
      }

      local_maxDelta = 0.0;

      // Iteraciones del método de Gauss-Seidel.
      for(int j=1; j < column-1; ++j){
        if((2*factor <= absolute_position && absolute_position < 8*factor) && (j == 2*factor || j == 8*factor)){
          continue;
        } else{
	  // Cálculo del nuevo valor que recibirá la celda.
          local_section[i*column + j] = ((1.0+omega)*0.25)*(local_section[(i+1)*column + j] + local_section[(i-1)*column + j] + local_section[i*column + (j+1)] + local_section[i*column + (j-1)]) - omega*local_section[i*column + j];
	  // Comparar el valor obtenido con el de la iteración anterior.
	  local_maxDelta = std::max(local_maxDelta, std::abs(local_section[i*column + j] - local_sectionCopy[i*column + j]));
        }
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
    its++;
    // Copia
    local_sectionCopy = local_section;   
  }
  double time_2 = MPI_Wtime();

  // Imprimir la matriz resultante.
  MPI_Request print_request;
  int signal = 1;
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

    MPI_Isend(&signal, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &print_request);
    MPI_Wait(&print_request, MPI_STATUS_IGNORE);
  } else{
    MPI_Irecv(&signal, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &print_request);
    MPI_Wait(&print_request, MPI_STATUS_IGNORE);

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
      MPI_Isend(&signal, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &print_request);
      MPI_Wait(&print_request, MPI_STATUS_IGNORE);
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

  MPI_Barrier(MPI_COMM_WORLD);
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
