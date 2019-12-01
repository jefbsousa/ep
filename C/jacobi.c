#include <stdio.h>
#include <math.h>

#define N 4
#define MAXITER 8000
#define h 1.0/N 
#define alpha 3.0


void print_matrix(double m[N+1][N+1]){
	for(int i=0; i <= N; i++){
		for(int j=0; j <= N; j++)
			printf("%lf ", m[i][j]);
		printf("\n");
	}
}

void initialize_matrix( double m[N+1][N+1], double value){
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m[i][j] = value;
}			

void initialize_U_answer( double m[N+1][N+1]){
	
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m[i][j] = (alpha * exp(h*i) * sin(h*j));
}	

void copy_m2_to_m1(double m_1[N+1][N+1], double m_2[N+1][N+1]){
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m_1[i][j] = m_2[i][j];
}


void update_border(double U_new[N+1][N+1], double U_old[N+1][N+1]){

	for(int i=0; i <= N; i++){
		U_old[i][0] = alpha * exp(h*i) * sin(h*0.0);
		U_old[i][N] = alpha * exp(h*i) * sin(h*N);
	}

	for(int j=0; j <= N; j++){
		U_old[0][j] = alpha * exp(h*0.0) * sin(h*j);
		U_old[N][j] = alpha * exp(h*N) * sin(h*j);
	}	

	copy_m2_to_m1(U_new, U_old);
}


void jacobi_method(double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]){
	for(int i=1; i <= N-1 ; i++)
		for(int j=1; j <= N-1 ; j++)
			U_new[i][j] =  0.25*(U_old[i-1][j] + U_old[i+1][j] + U_old[i][j-1] + U_old[i][j+1] + h*h *F[i][j]);
}



int main(int argc, char const *argv[])
{
	
	double U_new[N+1][N+1], U_old[N+1][N+1], U_answer[N+1][N+1], F[N+1][N+1];

	initialize_matrix(U_new, 0.0); initialize_matrix(U_old, 0.0); initialize_matrix(F, 0.0);
	initialize_U_answer(U_answer);
	update_border(U_new, U_old);


	for (int iter=0; iter <= 10; iter++){
		jacobi_method(U_new, U_old, F);
		/*
		if ( calculate_diff(U_new, U_old, N) <= TOL) {
			printf("Convergiu")
			break
		}
		*/
		copy_m2_to_m1(U_old, U_new);
	}

	printf("\n");
	print_matrix(U_answer);
	printf("\n");
	print_matrix(U_old);


	return 0;
}
