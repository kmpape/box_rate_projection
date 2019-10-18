#include <assert.h>
#include <stdio.h>
#include "box_rate_projection.h"


/* Some inputs/results for 2 dimensional projection */
float x0_2dim_1[2] = {3.0000, 0.0000};
float res_2dim_1[2] = {1.0000, 0.0000};
float x0_2dim_2[2] = {2.1213, 2.1213};
float res_2dim_2[2] = {1.0000, 1.0000};
float x0_2dim_3[2] = {0.0000, 3.0000};
float res_2dim_3[2] = {0.0000, 1.0000};
float x0_2dim_4[2] = {-2.1213, 2.1213};
float res_2dim_4[2] = {-0.5000, 0.5000};
float x0_2dim_5[2] = {-3.0000, 0.0000};
float res_2dim_5[2] = {-1.0000, 0.0000};
float x0_2dim_6[2] = {-3.0000, 0.0000};
float res_2dim_6[2] = {-1.0000, 0.0000};
float x0_2dim_7[2] = {-2.1213, -2.1213};
float res_2dim_7[2] = {-1.0000, -1.0000};
float x0_2dim_8[2] = {-0.0000, -3.0000};
float res_2dim_8[2] = {-0.0000, -1.0000};
float x0_2dim_9[2] = {2.1213, -2.1213};
float res_2dim_9[2] = {0.5000, -0.5000};
float x0_2dim_10[2] = {0.1, 0.1};
float res_2dim_10[2] = {0.1, 0.1};


/* Some inputs/results for DYKSTRA unit tests */
float x0_1[10] = {1.0166, 27.8734, -11.6667, -18.5430, -11.4068, -10.9334, -4.3361, -1.6847, -2.1853, 5.4133};
float res_1[10] = {1.0000, 1.0000, -0.0000, -1.0000, -1.0000, -1.0001, -1.0000, -1.0000, -0.0000, 1.0000};

float x0_2[10] = {8.8840, -11.4707, -10.6887, -8.0950, -29.4428, 14.3838, 3.2519, -7.5493, 13.7030, -17.1152};
float res_2[10] = {0.0001, -0.9999, -1.0001, -0.9999, -0.9998, 0.0003, -0.0003, -1.0001, -0.0001, -1.0000};

float x0_3[10] = {-1.0224, -2.4145, 3.1921, 3.1286, -8.6488, -0.3005, -1.6488, 6.2771, 10.9327, 11.0927};
float res_3[10] = {-1.0000, -0.1112, 0.8888, 0.0000, -1.0000, -0.3005, -0.0000, 1.0000, 0.9999, 1.0001};

float x0_4[10] = {-8.6365, 0.7736, -12.1412, -11.1350, -0.0685, 15.3263, -7.6967, 3.7138, -2.2558, 11.1736};
float res_4[10] = {-1.0000, 0.0000, -1.0000, -1.0000, -0.0000, 1.0000, -0.0000, 1.0000, -0.0000, 1.0000};

float x0_5[10] = {-10.8906, 0.3256, 5.5253, 11.0061, 15.4421, 0.8593, -14.9159, -7.4230, -10.6158, 23.5046};
float res_5[10] = {-1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 0.0000, -1.0000, -1.0000, 0.0000, 1.0000};

float x0_99[10] = {1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000};
float res_99[10] = {1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000};

float out[DIM];

float calc_error(const float * in1, const float * in2, const int len) {
	int i;
	float error = 0.0;
	for (i = 0; i < len; i++)
		error += (in1[i] - in2[i]) * (in1[i] - in2[i]);
	return error;
}

void debug_symmetric_dykstra(void) {
	float error;

	assert(DIM == 10);
	assert(A_MAX == 1.0);
	assert(R_MAX == 1.0);

	// 99
	//symmetric_box_rate_projection(x0_99, out);
	//error = calc_error(out, res_99, DIM);
	symmetric_box_rate_projection(x0_2, out);
	error = calc_error(out, res_2, DIM);
	printf("Trial 99 error = %.6f\n", error);
}

void test_symmetric_dykstra(void) {
	float error;

	assert(DIM == 10);
	assert(A_MAX == 1.0);
	assert(R_MAX == 1.0);

	// 1
	symmetric_box_rate_projection(x0_1, out);
	error = calc_error(out, res_1, DIM);
	printf("Trial 1 error = %.6f\n", error);

	// 2
	symmetric_box_rate_projection(x0_2, out);
	error = calc_error(out, res_2, DIM);
	printf("Trial 2 error = %.6f\n", error);

	// 3
	symmetric_box_rate_projection(x0_3, out);
	error = calc_error(out, res_3, DIM);
	printf("Trial 3 error = %.6f\n", error);

	// 4
	symmetric_box_rate_projection(x0_4, out);
	error = calc_error(out, res_4, DIM);
	printf("Trial 4 error = %.6f\n", error);

	// 5
	symmetric_box_rate_projection(x0_5, out);
	error = calc_error(out, res_5, DIM);
	printf("Trial 5 error = %.6f\n", error);
}

void test_2dim_projection(void) {
	float error;
	float out_2dim[2];
	assert(A_MAX == 1.0);
	assert(R_MAX == 1.0);

	// 1
	box_rate_proj_2dim(x0_2dim_1, out_2dim);
	error = calc_error(out_2dim, res_2dim_1, 2);
	printf("Trial 2DIM 1 error = %.6f\n", error);

	// 2
	box_rate_proj_2dim(x0_2dim_2, out_2dim);
	error = calc_error(out_2dim, res_2dim_2, 2);
	printf("Trial 2DIM 2 error = %.6f\n", error);

	// 3
	box_rate_proj_2dim(x0_2dim_3, out_2dim);
	error = calc_error(out_2dim, res_2dim_3, 2);
	printf("Trial 2DIM 3 error = %.6f\n", error);

	// 4
	box_rate_proj_2dim(x0_2dim_4, out_2dim);
	error = calc_error(out_2dim, res_2dim_4, 2);
	printf("Trial 2DIM 4 error = %.6f\n", error);

	// 5
	box_rate_proj_2dim(x0_2dim_5, out_2dim);
	error = calc_error(out_2dim, res_2dim_5, 2);
	printf("Trial 2DIM 5 error = %.6f\n", error);

	// 6
	box_rate_proj_2dim(x0_2dim_6, out_2dim);
	error = calc_error(out_2dim, res_2dim_6, 2);
	printf("Trial 2DIM 6 error = %.6f\n", error);

	// 7
	box_rate_proj_2dim(x0_2dim_7, out_2dim);
	error = calc_error(out_2dim, res_2dim_7, 2);
	printf("Trial 2DIM 7 error = %.6f\n", error);

	// 8
	box_rate_proj_2dim(x0_2dim_8, out_2dim);
	error = calc_error(out_2dim, res_2dim_8, 2);
	printf("Trial 2DIM 8 error = %.6f\n", error);

	// 9
	box_rate_proj_2dim(x0_2dim_9, out_2dim);
	error = calc_error(out_2dim, res_2dim_9, 2);
	printf("Trial 2DIM 9 error = %.6f\n", error);

	// 10
	box_rate_proj_2dim(x0_2dim_10, out_2dim);
	error = calc_error(out_2dim, res_2dim_10, 2);
	printf("Trial 2DIM 10 error = %.6f\n", error);
}

int main() {
	//debug_symmetric_dykstra();

	//test_2dim_projection();

	test_symmetric_dykstra();

	return 0;
}

