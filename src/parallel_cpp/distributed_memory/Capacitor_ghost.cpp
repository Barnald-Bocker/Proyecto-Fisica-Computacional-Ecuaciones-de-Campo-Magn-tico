#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

std::vector<double> gauss_seidel(double V1, double V2, double omega, int factor, double tolerance){


	if(tolerance <= 0){
		std::cerr << "El valor de tolerancia no puede ser nulo ni negativo!" << std::endl;
		exit(1);
	}
	double delta = 1.0;
	double diff;
	double maxDelta;
	int its = 0;
	int column = 10*factor+1;
  int size, rank;
  int row = 10*factor;
  MPI_Init(NULL, NULL);
  double time_1 = MPI_Wtime();
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
  //int end = start +  nlocal;
  int end = nlocal+1;
 
  std::vector<double> local_section;

  if(rank == 0 || rank == size - 1){
    local_section.resize((nlocal+1)*column); // MUY IMPORTANTE: GRILLAS FANTASMA 
    end = nlocal;
  }
  else{
    local_section.resize((nlocal+2)*column);
  }
  
  if (rank == size-1){
    end = nlocal-1;
  }
  
  //Establecer condiciones de frontera
  
  for(int i=1; i < end; ++i){
    int absolute_position = start+i-1;
    if (rank == 0){
      absolute_position = start + i;
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
  
  // Envio de las filas fantasma superiores.
  if (rank > 0){
    MPI_Request request_send, request_recv;
    MPI_Isend(&local_section[column], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request_send);
    MPI_Irecv(&local_section[end*column], column, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request_recv);
    MPI_Wait(&request_send, MPI_STATUS_IGNORE);  // Wait for send to complete
    MPI_Wait(&request_recv, MPI_STATUS_IGNORE);  // Wait for receive to complete
  }

  if (rank + 1 < size){
    MPI_Request request_send, request_recv;
    MPI_Isend(&local_section[end*column], column, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &request_send);
    MPI_Irecv(&local_section[0], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request_recv);
    MPI_Wait(&request_send, MPI_STATUS_IGNORE);  // Wait for send to complete
    MPI_Wait(&request_recv, MPI_STATUS_IGNORE);  // Wait for receive to complete
  }
  // Copia
  std::vector<double> local_sectionCopy = local_section;
  
  while(delta > tolerance){
    for(int i=1; i < end; ++i){
      int absolute_position = start+i-1;
      if (rank == 0){
        absolute_position = start + i;
      }
      for(int j=1; j < column-1; ++j){
        if((2*factor <= absolute_position && absolute_position < 8*factor) && (j == 2*factor || j == 8*factor)){
  			  continue;
  			} else{
  					local_section[i*column + j] = ((1.0+omega)*0.25)*(local_section[(i+1)*column + j] + local_section[(i-1)*column + j] + local_section[i*column + (j+1)] + local_section[i*column + (j-1)]) - omega*local_section[i*column + j];
  				}
  			}
  		}	
    maxDelta = 0.0;
    for(int i=1; i < end; ++i){
       int absolute_position = start+i-1;
       if (rank == 0){
         absolute_position = start + i;
       }
       for(int j=1; j < column-1; ++j){
         if((2*factor <= absolute_position && absolute_position < 8*factor) && (j == 2*factor || j == 8*factor)){
           continue;
         } else{
           diff = std::abs(local_section[i*column+j] - local_sectionCopy[i*column+j]);
           if(diff > maxDelta){
             maxDelta = diff;
           }

          }
        }
    }
    double global_delta;
    MPI_Allreduce(&maxDelta, &global_delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    delta = global_delta;

      // Comunicar el delta

    // Envio de las filas fantasma superiores.
    // Envio de las filas fantasma superiores.
    if (rank > 0){
      MPI_Request request_send, request_recv;
      MPI_Isend(&local_section[column], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request_send);
      MPI_Irecv(&local_section[end*column], column, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request_recv);
      MPI_Wait(&request_send, MPI_STATUS_IGNORE);  // Wait for send to complete
      MPI_Wait(&request_recv, MPI_STATUS_IGNORE);  // Wait for receive to complete
    }
    if (rank == 0){
  
    }
    if (rank + 1 < size){
      MPI_Request request_send, request_recv;
      MPI_Isend(&local_section[end*column], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request_send);
      MPI_Irecv(&local_section[0], column, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request_recv);
      MPI_Wait(&request_send, MPI_STATUS_IGNORE);  // Wait for send to complete
      MPI_Wait(&request_recv, MPI_STATUS_IGNORE);  // Wait for receive to complete
    }
    // Copia
    local_sectionCopy = local_section;
    its += 1;
  }
  //std::cout << delta << std::endl;
  switch(rank){
    case 0:
      {
      std::vector<double> grid(10*factor*(10*factor+1));
      int displacements[size];
      int counts[size];
      displacements[0] = 0; 
      for (int i=1; i<size-1;i++){
        displacements[i] =((row/size) + rest*( (i<(row%size)) ? 1 : 0))*i + rest*( (rest <= i) ? 1:0 );
        counts[size] = ((row/size) + rest*( (i<(row%size)) ? 1 : 0))*column;
      }
      MPI_Gatherv(&local_section,nlocal*column, MPI_DOUBLE, &grid, counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      double time_2 = MPI_Wtime();
      MPI_Finalize(); 
      //std::cout << "Numero de iteraciones: " << its << std::endl;
      return grid;
      break;
      }
    default:
      MPI_Gatherv(&local_section, nlocal*column, MPI_INT, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      break;
  }
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
  int factor = 1;
  double volt1 = 1.0;
  double volt2 = -1.0;
  double omega = 0.9;
  double tolerance = 1e-4;
  std::vector<double> grid = gauss_seidel(volt1,volt2,omega, factor, tolerance);
  print_grid(grid,factor);
}
