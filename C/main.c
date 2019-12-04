#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "functions.h"


void exercises(int length_array_N, int array_N[length_array_N], int exercise){

	for (int i = 0; i < length_array_N; i++){	
		printf("JACOBI");
		call_methods(array_N[i], 1, exercise);
		printf("\n");
	}
	
	for (int i = 0; i < length_array_N; i++){	
		printf("GAUSS-SEIDEL");
		call_methods(array_N[i], 2, exercise);
		printf("\n");
	}
	
	for (int i = 0; i < length_array_N; i++){	
		printf("SOR");
		call_methods(array_N[i], 3, exercise);
		printf("\n");
	}

}




int main(int argc, char const *argv[])
{

	int array_N[4] = {8, 16, 32, 64};
	/*Para pegar o tamanho do array, dividimos pelo tamanho de cada inteiro(4 bytes)*/
	int length_array_N = sizeof(array_N)/4;	

	//terceiro argumento é em relação ao exercício
	//1A -> 1; 1B -> 2; 2A -> 3; 2B -> 4;
	exercises(length_array_N, array_N, 1);
	exercises(length_array_N, array_N, 2);
	exercises(length_array_N, array_N, 3);
	exercises(length_array_N, array_N, 4);
	/*
	void (*test_ptr)(double, int) = &test_pointer;

	double val = 10.23232435442343; int val2 = 2;
	(*test_ptr)(val, val2);
	*/
	
	return 0;
}
