#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define MAXITER 100000
#define alpha 3.0
// 8, 16, 32, 64, 128, 256, 512


void print_matrix(int N, double m[N+1][N+1]){
	for(int i=0; i <= N; i++){	
		for(int j=0; j <= N; j++)
			printf("%lf ", m[i][j]);
		printf("\n");
	}
}



void initialize_matrix(int N, double m[N+1][N+1], double value){

	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m[i][j] = value;
}			


void initialize_U_answer_1A(int N, double m[N+1][N+1], double h){
	
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m[i][j] = (alpha * exp(h*i) * sin(h*j));
}	


void copy_m2_to_m1(int N, double m_1[N+1][N+1], double m_2[N+1][N+1]){
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m_1[i][j] = m_2[i][j];
}


void update_border(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){

	for(int i=0; i <= N; i++){
		U_old[i][0] = alpha * exp(h*i) * sin(h*0.0);
		U_old[i][N] = alpha * exp(h*i) * sin(h*N);
	}

	for(int j=0; j <= N; j++){
		U_old[0][j] = alpha * exp(h*0.0) * sin(h*j);
		U_old[N][j] = alpha * exp(h*N) * sin(h*j);
	}	

	//memcpy(U_old, U_new, (N+1)*(N+1));
	copy_m2_to_m1(N, U_new, U_old);

}


double calculate_diff(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){
	double sum = 0.0;
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			sum += pow(fabs(U_new[i][j] - U_old[i][j]), 2);

	return (h*h * sqrt(sum));
}


double calculate_max_abs_diff(int N, double U_answer[N+1][N+1], double U_new[N+1][N+1]){
	
	double max_value = 0.0, diff = 0.0;

	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++){
			diff = fabs(U_answer[i][j] - U_new[i][j]);

			if (diff > max_value)
				max_value = diff;
		}

	return max_value;
}


void jacobi_method(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]){
	for(int i=1; i <= N-1 ; i++)
		for(int j=1; j <= N-1 ; j++)
			U_new[i][j] =  0.25*(U_old[i-1][j] + U_old[i+1][j] + U_old[i][j-1] + U_old[i][j+1] + h*h *F[i][j]);
}

void g_s_method(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]){
	for(int i=1; i <= N-1 ; i++)
		for(int j=1; j <= N-1 ; j++)
			U_new[i][j] =  0.25*(U_new[i-1][j] + U_old[i+1][j] + U_new[i][j-1] + U_old[i][j+1] + h*h *F[i][j]);
}

void sor_method(int N, double h, double w, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]){
	for(int i=1; i <= N-1 ; i++)
		for(int j=1; j <= N-1 ; j++)
			U_new[i][j] =  (w/4.0) * (U_new[i-1][j] + U_old[i+1][j] + U_new[i][j-1] + U_old[i][j+1] + h*h *F[i][j]) + (1-w)*U_old[i][j];  
}


void iterative_method(int N, int option, double h, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1], double U_answer[N+1][N+1]){
	
	double TOL = h*0.00001;
	/* M_PI, constante definida em math.h	 */
	double w = 2.0/(1 + sin(M_PI*h));
	clock_t start = clock();

	for (int iter=0; iter < MAXITER; iter++){
		if (option == 1)
			jacobi_method(N, h, U_new, U_old, F);
		if (option == 2)
			g_s_method(N, h, U_new, U_old, F);
		if (option == 3)
			sor_method(N, h, w, U_new, U_old, F);
		
		if ( calculate_diff(N, h, U_new, U_old) <= TOL) {
			printf("\nConvergiu, N=%d, iteracoes=%d \n", N, iter);
			break;
		}
		
		//memcpy(U_old, U_new, (N+1)*(N+1));

		copy_m2_to_m1(N, U_old, U_new);
	}

	clock_t stop = clock();

	double max_value = calculate_max_abs_diff(N, U_answer, U_new);

	double elapsed = (double)(stop - start)  / CLOCKS_PER_SEC;
    printf("Tempo em segundos: %lf\nMax_value = %.9g\n", elapsed, max_value);

}


void call_methods(int N, int option){
	double h = 1.0/N;
	
	double (*U_answer)[N+1] = malloc((N+1) * sizeof(*U_answer));
	double (*U_new)[N+1] = malloc((N+1) * sizeof(*U_new));
	double (*U_old)[N+1] = malloc((N+1) * sizeof(*U_old));
	double (*F)[N+1] = malloc((N+1) * sizeof(*F));

	initialize_matrix(N, U_new, 0.0);
	initialize_matrix(N, U_old, 0.0);
	/* Parte A a matriz F é toda nula */
	initialize_matrix(N, F, 0.0);

	initialize_U_answer_1A(N, U_answer, h);
	
	update_border(N, h, U_new, U_old);
	
	iterative_method(N, option, h, U_new, U_old, F, U_answer);

	free(U_answer); free(U_new); free(U_old); free(F);
}



int main(int argc, char const *argv[])
{

	int N[4] = {8, 16, 32, 64};

	/*
	for (int i = 0; i < 4; ++i)
	{	
		printf("JACOBI");
		call_methods(N[i], 1);
		printf("\n");
	}
	
	for (int i = 0; i < 4; ++i)
	{	
		printf("GAUSS-SEIDEL");
		call_methods(N[i], 2);
		printf("\n");
	}
	*/

	for (int i = 0; i < 4; ++i)
	{	
		printf("SOR");
		call_methods(N[i], 3);
		printf("\n");
	}

	
	/*
	printf("\n");
	print_matrix(N, U_answer);
	printf("\n");
	print_matrix(N, U_old);
	printf("\n");
	print_matrix(N, U_new);	
	*/

	

	return 0;
}
