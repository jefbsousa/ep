
#include "functions.h"


void initialize_U_answer_1A(int N, double m[N+1][N+1], double h){
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m[i][j] = (alpha * exp(h*i) * sin(h*j));
}

void initialize_U_answer_1B(int N, double m[N+1][N+1], double h){
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m[i][j] = alpha * cos(M_PI* h*i) * sin(M_PI* h*j);
}

void update_F_1A(int N, double F[N+1][N+1]){

	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			F[i][j] = 0.0;
}

void update_F_1B(int N, double F[N+1][N+1]){
	double h = 1.0/N, constant = 02 * alpha * (M_PI*M_PI);

	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			F[i][j] = constant * cos(M_PI* h*i) * sin(M_PI* h*j);

}

void update_F_2A(int N, double F[N+1][N+1]){

	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			F[i][j] = -100./8.84 ;
}

void update_F_2B(int N, double F[N+1][N+1]){
	double h = 1.0/N;
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			F[i][j] = -10.0/75.0 * sin(M_PI * (i*h + j*h)) * pow(10, 4);
}

void update_border_1A(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){

	for(int i=0; i <= N; i++){
		U_old[i][0] = alpha * exp(h*i) * sin(h*0.0);
		U_old[i][N] = alpha * exp(h*i) * sin(h*N);
	}

	for(int j=0; j <= N; j++){
		U_old[0][j] = alpha * exp(h*0.0) * sin(h*j);
		U_old[N][j] = alpha * exp(h*N) * sin(h*j);
	}	

	copy_m2_to_m1(N, U_new, U_old);

}

void update_border_1B(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){

	for(int i=0; i <= N; i++){
		U_old[i][0] = alpha * cos(M_PI* h*i) * sin(M_PI* h*0.0);
		U_old[i][N] = alpha * cos(M_PI* h*i) * sin(M_PI* h*N);
	}

	for(int j=0; j <= N; j++){
		U_old[0][j] = alpha * cos(M_PI* h*0.0) * sin(M_PI* h*j);
		U_old[N][j] = alpha * cos(M_PI* h*N) * sin(M_PI* h*j);
	}	

	copy_m2_to_m1(N, U_new, U_old);

}

/* Ex. 2) A) e B) têm os mesmos valores de bordo*/
void update_border_2AB(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){
	/*Lados A e B no espaçamento entre as placas(Figura 1(b) do enunciado do EP)*/
	for(int i=0; i <= N; i++){
		U_old[i][0] = 0.0;
		U_old[i][N] = 110.0;
	}

	for(int j=0; j <= N; j++){
		U_old[0][j] = 110. * sin((M_PI/2.0)*(h*j));
		U_old[N][j] = 110. * sin((M_PI/2.0)*(h*j));
	}	
	copy_m2_to_m1(N, U_new, U_old);

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

void write_results_disk_1(int option, int N, int n_iter, double elapsed, double max_err_value){
	FILE *f = fopen("results_1.csv", "a");
	
	if (option == 1)
		fprintf(f, "Jacobi;%d;%d;%lf;%lf\n",N,n_iter,elapsed,max_err_value);
	if (option == 2)
		fprintf(f, "Gauss-Seidel;%d;%d;%lf;%lf\n",N,n_iter,elapsed,max_err_value);
	if (option == 3)
		fprintf(f, "SOR;%d;%d;%lf;%lf\n",N,n_iter,elapsed,max_err_value);

	fclose(f); 
}

void write_results_disk_2(int option, int N, int n_iter, double elapsed){
	FILE *f = fopen("results_2.csv", "a");
	
	if (option == 1)
		fprintf(f, "Jacobi;%d;%d;%lf\n",N,n_iter,elapsed);
	if (option == 2)
		fprintf(f, "Gauss-Seidel;%d;%d;%lf\n",N,n_iter,elapsed);
	if (option == 3)
		fprintf(f, "SOR;%d;%d;%lf\n",N,n_iter,elapsed);

	fclose(f); 
}
