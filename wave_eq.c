#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void ini_cond_sine(float *u_space_0, int space_steps, int n_of_zeroes);

void ini_cond_sharp_pulse(float *u_space_0, int space_steps, int intensity);

void ini_cond_gauss_pulse(float *u_space_0, int space_steps, float center, float width);

void seq_step(float *u_space_0,float *u_space_1, float *u_space_2, float C,int time_steps, int space_steps);

void seq_replace(float **u_space_0,float **u_space_1, float **u_space_2);

int main(){

	// Initialization of variables
	
	// Space steps are in mm, time steps are in ms.
	
	int length_space = 1; // In meters
	int length_time = 5; // In seconds
	
	int space_steps = 1024*length_space;
	int time_steps = 1000*length_time;
	int res = 100; // Number of frames we want to export
	
	float c = 2.5; // Speed of the wave
	float dx = 1.0 / (space_steps - 1);
	float dt = 0.8 * dx / c; // This is to ensure the CFT condition.
	float C = (c*dt/dx)*(c*dt/dx);
	
	float* u_space_0 = malloc(space_steps*sizeof(double));
	float* u_space_1 = malloc(space_steps*sizeof(double));
	float* u_space_2 = malloc(space_steps*sizeof(double));
	
	// Initial Conditions
	
	clock_t start_time = clock();

	ini_cond_gauss_pulse(u_space_0,space_steps,0.25,0.04);
	
	// First time iteration
	
	for (int i = 1; i<space_steps;i++){
		u_space_1[i] = u_space_0[i] + (C/2.0)*(u_space_0[i-1]-2*u_space_0[i]+u_space_0[i+1]);
	}
	u_space_1[0] = 0;
	u_space_1[space_steps-1] = 0;

	clock_t end_time = clock();
	float time_initializing = ((float)(end_time - start_time) / CLOCKS_PER_SEC);

	// Output of initial wave conditions
	
	FILE *output_file = fopen("output_wave.csv","w");
	
	for (int i = 0; i<space_steps;i++){
		fprintf(output_file,"%f,",u_space_0[i]);
	}
	fprintf(output_file,"\n");
	
	// Output of the first moment simulated.

	for (int i = 0; i<space_steps;i++){
		fprintf(output_file,"%f,",u_space_1[i]);
	}
	fprintf(output_file,"\n");
	
	// Time iterations onwards
	
	start_time = clock();
	
	for (int n = 2; n < time_steps; n++){
		seq_step(u_space_0, u_space_1, u_space_2, C, time_steps, space_steps);

		if (n % ((time_steps)/res) == 0){
			for (int m = 0; m < space_steps;m++){
			fprintf(output_file,"%f,",u_space_2[m]);
			};
			fprintf(output_file,"\n");
		}
		
		seq_replace(&u_space_0,&u_space_1,&u_space_2);
	}
	
	
	fclose(output_file);
	
	end_time = clock();
	float time_seq = ((float)(end_time - start_time) / CLOCKS_PER_SEC);
	
	FILE *timing_file = fopen("timing_results.txt","w");
	
	if (timing_file != NULL){
	fprintf(timing_file, "Wave Simulation Timing Results\n");
	fprintf(timing_file, "===============================\n");
	fprintf(timing_file, "Points: %d\n", space_steps);
	fprintf(timing_file, "Time steps: %d\n", time_steps);
	fprintf(timing_file, "Initialization time: %.6f seconds\n", time_initializing);
	fprintf(timing_file, "Sequential time: %.6f seconds\n", time_seq);
	fclose(timing_file);
	}
	
	free(u_space_0);
	free(u_space_1);
	free(u_space_2);
}

void ini_cond_sine(float *u_space_0, int space_steps, int n_of_zeroes){
	for (int i = 0; i<=space_steps;i++){
		u_space_0[i] = sin(n_of_zeroes * M_PI * (float)i/(space_steps-1));
	}
	u_space_0[0] = 0;
	u_space_0[space_steps-1] = 0;
}

void ini_cond_sharp_pulse(float *u_space_0, int space_steps, int intensity){
	for (int i = 0; i<=space_steps;i++){
		u_space_0[i] = 0;
	}
	u_space_0[space_steps/2] = intensity;
}

void ini_cond_gauss_pulse(float *u_space_0, int space_steps, float center, float width){
	for (int i = 0; i < space_steps; i++) {
		float x = (float)i / (space_steps - 1);
		u_space_0[i] = exp(-powf(x - center, 2) / (2 * powf(width, 2)));
	}
}

void seq_step(float *u_space_0,float *u_space_1, float *u_space_2, float C,int time_steps, int space_steps){
		
	u_space_2[0] = 0;
	u_space_2[space_steps-1] = 0;
	
	for (int i = 1; i < space_steps-1;i++){
		u_space_2[i] = -u_space_0[i] + 2.0*u_space_1[i] + (C)*(u_space_1[i-1]+u_space_1[i+1]-2.0*u_space_1[i]);
	}
	
	u_space_2[0] = 0;
	u_space_2[space_steps-1] = 0;
}

void seq_replace(float **u_space_0,float **u_space_1, float **u_space_2){
	
	float* temp = *u_space_0;
	*u_space_0 = *u_space_1;
	*u_space_1 = *u_space_2;
	*u_space_2 = temp;
}

// Source: https://hplgit.github.io/num-methods-for-PDEs/doc/pub/wave/pdf/wave-4print-A4-2up.pdf

