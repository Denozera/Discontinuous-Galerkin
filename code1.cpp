#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;
// Interval for the problem is assumed to be [1 , 1]
// Degree of polynomial is taken to be 2.

double zeta[3];					// Chebyshev Nodes

double basis[3][3];				// Basis Psi by Lagrange interpolation 

int i = 0;
int j = 0;

void chebyshevNodes();
void basisCoefficient();

int main(){

	chebyshevNodes();

	basisCoefficient();

	for(i=0; i<3; i++){

		cout << basis[0][i] << "	" << basis[1][i] << "	" << basis[2][i] <<endl;
	}

	cout << zeta[0] << "	" << zeta[1] << "	" << zeta[2] <<endl;

}

void chebyshevNodes(){
	int n = 3;

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

	double delta;

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

