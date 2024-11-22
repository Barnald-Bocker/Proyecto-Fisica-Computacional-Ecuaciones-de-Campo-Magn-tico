#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> grid_generator(int factor, double V1, double V2, double boundary = 0.0){

	if(factor < 1){
	    std::cerr << "El factor de escala debe ser un número entero mayor que cero!" << std::endl;
	    exit(1);
	}

	int column = 10*factor+1;
	std::vector<double> phi(10*factor*column); 
	for(int i=2*factor; i < 8*factor; ++i){
		phi[i*column+2*factor] = V1;
		phi[i*column+8*factor] = V2;
	}

	return phi;
}

int gauss_seidel(std::vector<double>& phi, double omega, int factor, double tolerance){

	if(phi.size() != (10*factor)*(10*factor+1)){
		std::cerr << "El factor de escala ingresado es diferente al utilizado para generar la grilla!" << std::endl;
		exit(1);
	}

	if(tolerance <= 0){
		std::cerr << "El valor de tolerancia no puede ser nulo ni negativo!" << std::endl;
		exit(1);
	}

	std::vector<double> phiCopy = phi;
	double delta = 1.0;
	double diff;
	double maxDelta;
	int its = 0;
	int column = 10*factor+1;
	
	while(delta > tolerance){
		for(int i=1; i < 10*factor-1; ++i){
			for(int j=1; j < column-1; ++j){
				if((2*factor <= i && i < 8*factor) && (j == 2*factor || j == 8*factor)){
					continue;
				} else{
					phi[i*column + j] = ((1.0+omega)*0.25)*(phi[(i+1)*column + j] + phi[(i-1)*column + j] + phi[i*column + (j+1)] + phi[i*column + (j-1)]) - omega*phi[i*column + j];
				}
			}
		}
		
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
  int factor = 1;
  double volt1 = 1.0;
  double volt2 = -1.0;
  std::vector<double> grid = grid_generator(factor, volt1, volt2);
  double omega = 0.9;
  double tolerance = 1e-4;
  int iterations = gauss_seidel(grid,omega, factor, tolerance);
  print_grid(grid,factor);
}
