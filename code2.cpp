#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;
// Interval for the problem is assumed to be [-1 , 1] discritized with N cells
// Degree of polynomial is taken to be 2.

double zeta[3];					// Chebyshev Nodes

double basis[3][3];				// Basis Psi by Lagrange interpolation 
double basisD[3][3];
double mass[3][3];
double stiffness[3][3];
double inverseM[3][3];
double psiInterfaceValue[3][N+1];
double fCap[3][N];

double integral[9];
double a[j];					// The time dependant part of U
double psiFlux[N+1];			// Because since there are N cells between -1 to 1, we need to calculate flux at N+1 points (Though total points are N+3)
double delta = 0;

int i = 0;
int j = 0;
int l = 0;
int m = 0;

double x, xL, xR;

void chebyshevNodes();
void basisCoefficient();
void polyIntegration();
void massMatrix();
void stiffnessMatrix();
void inverseMatrix();
void psiInterface();

int main(){

	chebyshevNodes();

	basisCoefficient();

	for(i=0; i<3; i++){

		cout << basis[0][i] << "	" << basis[1][i] << "	" << basis[2][i] <<endl;
	}

	cout << zeta[0] << "	" << zeta[1] << "	" << zeta[2] <<endl;

}

void chebyshevNodes(){
	//int n = 3;

	for(i=0; i<3; i++){

		if((double)(2*i+1)/3 == 1){
			zeta[i] = 0;
		}
		else{
			zeta[i] = -cos((double)(2*i+1)/3*0.5*M_PI);
		}
		
	}

}

void basisCoefficient(){

	delta = (zeta[1]-zeta[0])*(zeta[2]-zeta[0])*(zeta[2]-zeta[1]);

	for(i=0; i<3; i++){
		basis[0][i] = zeta[(i+1)%3]*zeta[(i+2)%3]*(zeta[(i+2)%3]-zeta[(i+1)%3]);
	}

	for(i=0; i<3; i++){
		basis[1][i] = (zeta[(i+1)%3]+zeta[(i+2)%3])*(zeta[(i+2)%3]-zeta[(i+1)%3]);
	}

	for(i=0; i<3; i++){
		basis[2][i] = zeta[(i+2)%3]-zeta[(i+1)%3];
	}

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			basis[i][j] = (double)basis[i][j]/delta;
		}
	}
}

void polyIntegration(){

	for(i=0; i<9; i++){

		integral[i] = (double)1/(i+1) * (pow(xR,i+1) - pow(xL,i+1));
	}
}

void massMatrix(){
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

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			mass[i][j] = 0;
		}
	}


	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			for(l=0; l<3; l++){
				for(m=0; m<3; m++){
					mass[i][j] = mass[i][j] + basis[l][i]*basis[m][j]*integral[l+m];
				}
			}
		}
	}
}


void stiffnessMatrix(){

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			stiffness[i][j] = 0;
		}
	}

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			basisD[i][j] = i*basis[i][j];
		}
	}

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			for(l=0; l<3; l++){
				for(m=0; m<3; m++){

					if(m != 0){
						stiffness[i][j] = stiffness[i][j] + basis[l][i]*basisD[m][j]*integral[l+m-1];
					}
					
				}
			}
		}
	}

}


void inverseMatrix(){

	for(i=0; i<3; i++){
        for(j=0; j<3; j++){

           	inverseM[i][j] = ((mass[(j+1)%3][(i+1)%3]*mass[(j+2)%3][(i+2)%3]) - (mass[(j+1)%3][(i+2)%3]*mass[(j+2)%3][(i+1)%3]))/delta;

        }
    }

}
   
void psiInterface(){

	for(j=0; j<3; j++){
		for(i=0; i<N+1; i++){
			psiInterfaceValue[j][i] = basis[0][j] + basis[1][j]*(-1 + (double)2/N*i) + basis[2][j]*(-1 + (double)2/N*i)*(-1 + (double)2/N*i);    // First row will represent all the values of psi 0 and so on
		}
	}
	
}


void fCapCalc(){

	for(j=0; j<3; j++){
		for(i=0; i<N; i++){
			fCap[j][i] = psiInterfaceValue[j][i+1]*psiInterfaceValue[j][i+1]*a[j][i+1] - psiInterfaceValue[j][i]*psiInterfaceValue[j][i]*a[j][i];  
		}
	}
}

