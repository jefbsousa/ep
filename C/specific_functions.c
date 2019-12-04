
#include "functions.h"


void initialize_U_answer_1A(int N, double m[N+1][N+1], double h){
	
	for(int i=0; i <= N; i++)
		for(int j=0; j <= N; j++)
			m[i][j] = (alpha * exp(h*i) * sin(h*j));
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
		U_old[i][0] = alpha * exp(h*i) * sin(h*0.0);
		U_old[i][N] = alpha * exp(h*i) * sin(h*N);
	}

	for(int j=0; j <= N; j++){
		U_old[0][j] = alpha * exp(h*0.0) * sin(h*j);
		U_old[N][j] = alpha * exp(h*N) * sin(h*j);
	}	

	copy_m2_to_m1(N, U_new, U_old);

}

void update_border_2A(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){

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

void update_border_2B(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]){

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
