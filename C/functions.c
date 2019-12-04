
#include "functions.h"


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


void copy_m2_to_m1(int N, double m_1[N+1][N+1], double m_2[N+1][N+1]){
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m_1[i][j] = m_2[i][j];
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




double calculate_diff(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){
	double sum = 0.0;
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			sum += pow(fabs(U_new[i][j] - U_old[i][j]), 2);

	return (h*h * sqrt(sum));
}

void write_results_disk_1A(int option, int N, int n_iter, double elapsed, double max_err_value){
	FILE *f = fopen("/home/jeferson/usp/edp/ep/C/results_1A.csv", "a");
	
	if (option == 1)
		fprintf(f, "Jacobi;%d;%d;%lf;%lf\n",N,n_iter,elapsed,max_err_value);
	if (option == 2)
		fprintf(f, "Gauss-Seidel;%d;%d;%lf;%lf\n",N,n_iter,elapsed,max_err_value);
	if (option == 3)
		fprintf(f, "SOR;%d;%d;%lf;%lf\n",N,n_iter,elapsed,max_err_value);

	fclose(f); 
}


void iterative_method(int N, int option, double h, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1], double U_answer[N+1][N+1]){
	
	/* M_PI, constante definida em math.h	 */
	double TOL = h*0.00001, w = 2.0/(1 + sin(M_PI*h));
	int n_iter = 0;

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
			n_iter = iter;
			break;
		}

		copy_m2_to_m1(N, U_old, U_new);
		n_iter = iter;
	}

	clock_t stop = clock();

	double max_value = calculate_max_abs_diff(N, U_answer, U_new);

	double elapsed = (double)(stop - start)  / CLOCKS_PER_SEC;
    printf("Tempo em segundos: %lf\nMax_value = %.9g\n", elapsed, max_value);

    //write_results_disk_1A(option,  N,  n_iter,  elapsed,  max_value);
}


void call_methods(int N, int option, int exercise_opt){
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
	
	void (*update_border_ptr)(int, double, double[N+1][N+1], double[N+1][N+1]);
	if (exercise_opt == 1)
		update_border_ptr = &update_border_1A;
	


	(*update_border_ptr)(N, h, U_new, U_old);
	

	iterative_method(N, option, h, U_new, U_old, F, U_answer);

	free(U_answer); free(U_new); free(U_old); free(F);
}


void test_pointer(double arg, int arg2){
	printf("%.9g, %d\n", arg, arg2);
}

