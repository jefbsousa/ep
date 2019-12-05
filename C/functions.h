#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAXITER 150000
#define alpha 3.0


void print_matrix(int N, double m[N+1][N+1]);

void initialize_matrix(int N, double m[N+1][N+1], double value);

void initialize_U_answer_1A(int N, double m[N+1][N+1], double h);
void initialize_U_answer_1B(int N, double m[N+1][N+1], double h);

void copy_m2_to_m1(int N, double m_1[N+1][N+1], double m_2[N+1][N+1]);

double calculate_max_abs_diff(int N, double U_answer[N+1][N+1], double U_new[N+1][N+1]);

void jacobi_method(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]);
void g_s_method(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]);
void sor_method(int N, double h, double w, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]);

void update_border_1A(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]);
void update_border_1B(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]);
void update_border_2AB(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]);

double calculate_diff(int N, double h, double U_new[N+1][N+1], double U_old[N+1][N+1]);

void write_results_disk_1(int option, int N, int n_iter, double elapsed, double max_err_value);
void write_results_disk_2(int option, int N, int n_iter, double elapsed);

void iterative_method(int N, int option, int exercise_opt, double U_new[N+1][N+1], double U_old[N+1][N+1], double F[N+1][N+1]);
void call_methods(int N, int option, int exercise_opt);

void update_F_1A(int N, double F[N+1][N+1]);
void update_F_1B(int N, double F[N+1][N+1]);
void update_F_2A(int N, double F[N+1][N+1]);
void update_F_2B(int N, double F[N+1][N+1]);