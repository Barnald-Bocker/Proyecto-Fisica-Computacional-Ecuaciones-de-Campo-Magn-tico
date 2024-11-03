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
				if(2*factor <= i && i < 8*factor && (j == 2*factor || j == 8*factor)){
					;
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


int jacobi_relaxation(std::vector<double>& phi, int factor, double tolerance){

	 if(phi.size() != (10*factor)*(10*factor+1)){
                std::cerr << "El factor de escala ingresado es diferente al utilizado para generar la grilla!" << std::endl;
                exit(1);
        }

        if(tolerance <= 0){
                std::cerr << "El valor de tolerancia no puede ser nulo ni negativo!" << std::endl;
                exit(1);
        }

	int column = 10*factor+1;
	std::vector<double> phiPrime = phi;
	std::vector<double> temp(10*factor*column);
        double delta = 1.0;
        double diff;
        double maxDelta;
        int its = 0;

        while(delta > tolerance){
                for(int i=1; i < 10*factor-1; ++i){
                        for(int j=1; j < column-1; ++j){
                                if(2*factor <= i && i < 8*factor && (j == 2*factor || j == 8*factor)){
                                        ;
                                } else{
                                        phiPrime[i*column + j] = 0.25*(phi[(i+1)*column + j] + phi[(i-1)*column + j] + phi[i*column + (j+1)] + phi[i*column + (j-1)]);
                                }
                        }
                }

		maxDelta = 0.0;
                for(int i=1; i < 10*factor-1; ++i){
                        for(int j=1; j < column-1; ++j){
                                diff = std::abs(phi[i*column + j] - phiPrime[i*column + j]);
                                if(diff > maxDelta){
                                        maxDelta = diff;
                                }
                        }
                }

		temp = phi;
		phi = phiPrime;
		phiPrime = temp;
                delta = maxDelta;
                its += 1;
        }

	phi = phiPrime;
        return its;
}


int SOR(std::vector<double>& phi, int factor, double omega, double tolerance){

         if(phi.size() != (10*factor)*(10*factor+1)){
                std::cerr << "El factor de escala ingresado es diferente al utilizado para generar la grilla!" << std::endl;
                exit(1);
        }

        if(tolerance <= 0){
                std::cerr << "El valor de tolerancia no puede ser nulo ni negativo!" << std::endl;
                exit(1);
        }

	int column = 10*factor+1;
        std::vector<double> phiPrime = phi;
        std::vector<double> temp(10*factor*column);
        double delta = 1.0;
        double diff;
        double maxDelta;
        int its = 0;

        while(delta > tolerance){
                for(int i=1; i < 10*factor-1; ++i){
                        for(int j=1; j < column-1; ++j){
                                if(2*factor <= i && i < 8*factor && (j == 2*factor || j == 8*factor)){
                                        ;
                                } else{
                                        phiPrime[i*column + j] = (1.0+omega)*(phi[(i+1)*column + j] + phi[(i-1)*column + j] + phi[i*column + (j+1)] + phi[i*column + (j-1)])*0.25 - omega*phi[i*column + j];
                                }
                        }
                }

                maxDelta = 0.0;
                for(int i=1; i < 10*factor-1; ++i){
                        for(int j=1; j < column-1; ++j){
                                diff = std::abs(phi[i*column + j] - phiPrime[i*column + j]);
                                if(diff > maxDelta){
                                        maxDelta = diff;
                                }
                        }
                }

                temp = phi;
                phi = phiPrime;
                phiPrime = temp;
                delta = maxDelta;
                its += 1;
        }

        phi = phiPrime;
        return its;
}
