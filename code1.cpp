#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;

double zeta[3];					// Three nodes for polynomial of 2nd order

double basis[3][3];

double basis1[3];
double basis2[3];
double basis3[3];

int i = 0;
int j = 0;

int main(){
	return 0;
}

void basisCoefficient(){

	double delta;

	delta = (zeta[1]-zeta[0])*(zeta[2]-zeta[0])*(zeta[2]-zeta[1]);

	for(i=0; i<3; i++){
		basis1[0][i] = zeta[(i+1)%3]*zeta[(i+2)%3]*(zeta[(i+2)%3]-zeta[(i+1)%3]);
	}

	for(i=0; i<3; i++){
		basis2[1][i] = (zeta[(i+1)%3]+zeta[(i+2)%3])*(zeta[(i+2)%3]-zeta[(i+1)%3]);
	}

	for(i=0; i<3; i++){
		basis3[2][i] = zeta[(i+2)%3]-zeta[(i+1)%3];
	}

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			basis[i][j] = (double)basis[i][j]/delta;
		}
	}
}