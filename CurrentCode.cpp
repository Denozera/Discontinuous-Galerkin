#include <iostream>
#include <fstream>
#include <cmath>
#include "mesh.h"


using namespace std;
// Interval for the problem is assumed to be [-1 , 1] discritized with N cells
// Degree of polynomial is taken to be 2.
// Forgot to calculate delta of mass matrix 

// For the time being, lets say that we want to do the simulation for 80 cells with CFL number 0.2 (just to compare the results)
double zeta[3][N+2];					// Chebyshev Nodes

double basis[3][3][N+2];				// Basis Psi by Lagrange interpolation 
double basisI[3][3][N+2];
double basisD[3][3][N+2];
double mass[3][3][N+2];
double stiffness[3][3][N+2];
double inverseM[3][3][N+2];
double psiInterfaceValueL[3][N+3];
double psiInterfaceValueR[3][N+3];
//double interfaceVelocityFieldL[N];
double interfaceVelocityFieldR[N+2];
double fCap[3][N+2];
double ka[3][N+2];				// The Multiplication of stiffness matrix and time dependent Psi Coefficients
double rhs[3][N+2];				// Final evaluation of RHS
double deltaM[N+2];
double delta[N+2];
//double flux[N+1];

double integral[9];
double a[3][N+2];					// The time dependant part of U (Initial condition)
//double psiFlux[N+1];			// Because since there are N cells between -1 to 1, we need to calculate flux at N+1 points (Though total points are N+3)
double deltaT = 0.0004;			// For CFL to be 0.5

int i = 0;
int j = 0;
int l = 0;
int m = 0;
int k = 0;

double x, t, xL, xR;

void chebyshevNodes();
void basisCoefficient();
//double velocityField(int, int);
void initialization();
void polyIntegration();
void massMatrix();
void stiffnessMatrix();
void inverseMatrix();
void psiInterface();
void interfaceVelocityField();
void FCapCalc();
void EvaluatingRSH();
void Euler();
void print();

int main(){

		chebyshevNodes();

		initialization();

		basisCoefficient();

		massMatrix();

		stiffnessMatrix();

		psiInterface();

	//for(t=0; t<0.2; t += deltaT){	

		interfaceVelocityField();

		FCapCalc();

		EvaluatingRSH();

		Euler();

	//}


	print();

	//cout << u << endl;
	/*
	chebyshevNodes();
	cout << zeta[0][0] << "	" << zeta[1][0] << "	" << zeta[2][0] <<endl;
	*/
}


void chebyshevNodes(){
	//int n = 3;

	for(i=0; i<N+2; i++){
			
		xL = -1 + (double)2/N*(i-1);
		xR = xL + (double)2/N;

		for(j=0; j<3; j++){

			if((double)(2*j+1)/3 == 1){
				zeta[j][i] = (xL+xR)*0.5;
			}
			else{
				zeta[j][i] = (xL+xR)*0.5 - 0.5*(xR-xL)*cos((double)(2*j+1)/3*0.5*M_PI);
			}
		}
	}

}

void basisCoefficient(){

	/*for(j=0; j<N+2; j++){

		//delta[j] = (zeta[1][j]-zeta[0][j])*(zeta[2][j]-zeta[0][j])*(zeta[2][j]-zeta[1][j]);

		for(i=0; i<3; i++){
			basis[0][i][j] = zeta[(i+1)%3][j]*zeta[(i+2)%3][j]*(zeta[(i+2)%3][j]-zeta[(i+1)%3][j]);
		}

		for(i=0; i<3; i++){
			basis[1][i][j] = (zeta[(i+1)%3][j]+zeta[(i+2)%3][j])*(zeta[(i+2)%3][j]-zeta[(i+1)%3][j]);
		}

		for(i=0; i<3; i++){
			basis[2][i][j] = zeta[(i+2)%3][j]-zeta[(i+1)%3][j];
		}

	}*/

	for(j=0; j<N+2; j++){
		for(i=0; i<3; i++){
			for(l=0; l<3; l++){
				basisI[i][l][j] = 1 * pow(zeta[i][j], l);
			}
		}
	}

	for(k=0; k<N+2; k++){	
		delta[k] = 0;
		for(i = 0; i < 3; i++){
			delta[k] = delta[k] + (basisI[0][i][k] * (basisI[1][(i+1)%3][k] * basisI[2][(i+2)%3][k] - basisI[1][(i+2)%3][k] * basisI[2][(i+1)%3][k]));
		}
	}

	for(k=0; k<N+2; k++){
		for(i=0; i<3; i++){
	        for(j=0; j<3; j++){
	           	basis[i][j][k] = ((basisI[(j+1)%3][(i+1)%3][k]*basisI[(j+2)%3][(i+2)%3][k]) - (basisI[(j+1)%3][(i+2)%3][k]*basisI[(j+2)%3][(i+1)%3][k]))/delta[k];
	        }
	    }
	}

	/*for(j=0; j<N+2; j++){
		for(i=0; i<3; i++){
			for(l=0; l<3; l++){
				basis[i][l][j] = (double)basis[i][l][j]/1000000;
			}
		}
	}*/
}

/*double velocityField(i, j){
	//initinally we want to make a matrix such that each row contains value of psi evaluated at that zeta.

	//psiZeta[0][0] = basis[0][0] + basis[1][0]*zeta[0] + basis[2][0]*zeta[0]*zeta[0];

	if(abs(zeta[j][i]) < (double)1/3){
		u = 2;
	}
	else{
		u = 1;
	}


	return u;
	

}*/

void initialization(){

	//double cons = (double)4/(sqrt(2*M_PI));

	for(i=0; i<N+2; i++){
		for(j=0; j<3; j++){

			a[j][i] =  (double)8/sqrt(2*M_PI) * exp(-32 * pow(zeta[j][i], 2));
		}
	}
		
}

void polyIntegration(){

	for(i=0; i<9; i++){

		integral[i] = (double)1/(i+1) * (pow(xR,i+1) - pow(xL,i+1));
	}
}

void massMatrix(){

	double deltaM[N+2];
/*	
	mass[0][0] = basis[0][0]*basis[0][0]*integral[0] + basis[0][0]*basis[1][0]*integral[1] + basis[0][0]*basis[2][0]*integral[2];
	mass[0][0] = mass[0][0] + basis[1][0]*basis[0][0]*integral[1] + basis[1][0]*basis[1][0]*integral[2] + basis[1][0]*basis[2][0]*integral[3];
	mass[0][0] = mass[0][0] + basis[2][0]*basis[0][0]*integral[2] + basis[2][0]*basis[1][0]*integral[3] + basis[2][0]*basis[2][0]*integral[4];

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			mass[i][j] = basis[0][i]*basis[0][j]*integral[0] + basis[0][i]*basis[1][j]*integral[1] + basis[0][i]*basis[2][j]*integral[2];
			mass[i][j] = mass[i][j] + basis[1][i]*basis[0][j]*integral[1] + basis[1][i]*basis[1][j]*integral[2] + basis[1][i]*basis[2][j]*integral[3];
			mass[i][j] = mass[i][j] + basis[2][i]*basis[0][j]*integral[2] + basis[2][i]*basis[1][j]*integral[3] + basis[2][i]*basis[2][j]*integral[4];
		}
	}

*/

	

	for(k=0; k<N+2; k++){

		deltaM[k] = 0;

		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				mass[i][j][k] = 0;

			}
		}
	}

	for(k=0; k<N+2; k++){

		xL = -1 + (double)2/N*(k-1);
		xR = xL + (double)2/N;

		for(i=0; i<9; i++){

		integral[i] = (double)1/(i+1) * (pow(xR,i+1) - pow(xL,i+1));

		}
		
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				for(l=0; l<3; l++){
					for(m=0; m<3; m++){
						mass[i][j][k] = mass[i][j][k] + basis[l][i][k]*basis[m][j][k]*integral[l+m];         //	Caution, See if xL and xR are explicitly specified in the main loop for each iteration so that we don't have to use multi-dimentional array for Intergral. 
					}
				}
			}
		}
	}

	for(k=0; k<N+2; k++){	
		for(i = 0; i < 3; i++){
			deltaM[k] = deltaM[k] + (mass[0][i][k] * (mass[1][(i+1)%3][k] * mass[2][(i+2)%3][k] - mass[1][(i+2)%3][k] * mass[2][(i+1)%3][k]));
		}
	}

	for(k=0; k<N+2; k++){
		for(i=0; i<3; i++){
	        for(j=0; j<3; j++){
	           	inverseM[i][j][k] = ((mass[(j+1)%3][(i+1)%3][k]*mass[(j+2)%3][(i+2)%3][k]) - (mass[(j+1)%3][(i+2)%3][k]*mass[(j+2)%3][(i+1)%3][k]))/deltaM[k];
	        }
	    }
	}
}


void stiffnessMatrix(){

	for(k=0; k<N+2; k++){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				stiffness[i][j][k] = 0;
			}
		}
	}

	for(k=0; k<N+2; k++){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				basisD[i][j][k] = i*basis[i][j][k];
			}
		}
	}

	for(k=0; k<N+2; k++){
		xL = -1 + (double)2/N*(k-1);
		xR = xL + (double)2/N;

		for(i=0; i<9; i++){

		integral[i] = (double)1/(i+1) * (pow(xR,i+1) - pow(xL,i+1));

		}
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				for(l=0; l<3; l++){
					for(m=1; m<3; m++){
							stiffness[i][j][k] = stiffness[i][j][k] + basis[l][i][k]*basisD[m][j][k]*integral[l+m-1];	
					}
				}
			}
		}
	}

}

/*
void inverseMatrix(){

	for(k=0; k<N; k++){
		for(i=0; i<3; i++){
	        for(j=0; j<3; j++){
	           	inverseM[i][j][k] = ((mass[(j+1)%3][(i+1)%3][k]*mass[(j+2)%3][(i+2)%3][k]) - (mass[(j+1)%3][(i+2)%3][k]*mass[(j+2)%3][(i+1)%3][k]))/deltaM[k];
	        }
	    }
	}

}
*/

void psiInterface(){

	for(j=0; j<3; j++){
		for(i=0; i<N+3; i++){
			psiInterfaceValueL[j][i] = basis[0][j][i] + basis[1][j][i]*(-1 + (double)2/N*(i-1)) + basis[2][j][i]*(-1 + (double)2/N*(i-1))*(-1 + (double)2/N*(i-1));    // First row will represent all the values of psi 0 and so on
			psiInterfaceValueR[j][i] = basis[0][j][i] + basis[1][j][i]*(-1 + (double)2/N*i) + basis[2][j][i]*(-1 + (double)2/N*i)*(-1 + (double)2/N*i);
		}
	}
	
}

void interfaceVelocityField(){

	for(i=0; i<N+2; i++){
		interfaceVelocityFieldR[i] = 0;
	}

	for(i=0; i<N+2; i++){
		for(j=0; j<3; j++){
			//interfaceVelocityFieldL[i] = a[j][i]*psiInterfaceValue[j][i];
			interfaceVelocityFieldR[i] = interfaceVelocityFieldR[i] + a[j][i]*psiInterfaceValueR[j][i];
		}
	}
}


void FCapCalc(){

	for(i=1; i<N+1; i++){
		for(j=0; j<3; j++){
			fCap[j][i] =  interfaceVelocityFieldR[i]*psiInterfaceValueR[j][i] - interfaceVelocityFieldR[i-1]*psiInterfaceValueL[j][i];  // Here, in order not to get negetive term in the square brackets, we need to add ghost cells.
		}
	}
}


void EvaluatingRSH(){

	for(i=1; i<N+1; i++){
		for(j=0; j<3; j++){
			ka[j][i] = stiffness[j][0][i]*a[0][i] + stiffness[j][1][i]*a[1][i] + stiffness[j][2][i]*a[2][i];  //This will give us product of Stiffness Matrix and a's
		}
	}

	//Now we have inverse of Mass Matrix, ka and fCap.

	for(i=1; i<N+1; i++){
		for(j=0; j<3; j++){
			rhs[j][i] = inverseM[j][0][i]*(ka[0][i]-fCap[0][i]) + inverseM[j][1][i]*(ka[1][i]-fCap[1][i]) + inverseM[j][2][i]*(ka[2][i]-fCap[2][i]);  //This will give us the derivative of a
		}
	}
}

void Euler(){

	for(i=1; i<N+1; i++){
		for(j=0; j<3; j++){
			a[j][i] = a[j][i] + deltaT*rhs[j][i];
		}
	}

	for(j=0; j<3; j++){
		a[j][0] = a[j][N-1];
		a[j][N+1] = a[j][2];
	}	
}
   
void print(){
	ofstream output;

        output.open("IC.csv");
        for(i=1; i<N+1; i++){
            for(j=0; j<3; j++){
            	output << zeta[j][i] << " " << a[j][i] << endl;
            }
        }
        output.close();


        /*output.open("sec2.csv");
        for(i=1; i<N+1; i++){
        	output << zeta[0][i] << " " << zeta[1][i] << " " << zeta[2][i] << endl;
        	output << endl;
            for(j=0; j<3; j++){
            	output << stiffness[j][0][i] << " " << stiffness[j][1][i] << " " << stiffness[j][2][i] << endl;;
            }
            output << endl;
        }
        output.close();*/
    
		output.open("sec4.csv");
        for(i=1; i<N+1; i++){
    		output << ka[0][i] << " " << ka[1][i] << " " <<	ka[2][i] << endl;
			output << interfaceVelocityFieldR[i] << endl;
			for(j=0; j<3; j++){
            	output << psiInterfaceValueR[j][i] << " " ;
            }
            output << endl;

            output << endl;

        }
        output.close();
        
}   