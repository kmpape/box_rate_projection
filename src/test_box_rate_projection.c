#include <assert.h>
#include <stdio.h>
#include "box_rate_projection.h"
#include "brp_mpc_test_data.h"

#define ABS_ERROR_TOL (1e-4f)

/* Some inputs/results for 2 dimensional projection */
brp_float x0_2dim_1[2] = {3.0000, 0.0000};
brp_float res_2dim_1[2] = {1.0000, 0.0000};
brp_float x0_2dim_2[2] = {2.1213, 2.1213};
brp_float res_2dim_2[2] = {1.0000, 1.0000};
brp_float x0_2dim_3[2] = {0.0000, 3.0000};
brp_float res_2dim_3[2] = {0.0000, 1.0000};
brp_float x0_2dim_4[2] = {-2.1213, 2.1213};
brp_float res_2dim_4[2] = {-0.5000, 0.5000};
brp_float x0_2dim_5[2] = {-3.0000, 0.0000};
brp_float res_2dim_5[2] = {-1.0000, 0.0000};
brp_float x0_2dim_6[2] = {-3.0000, 0.0000};
brp_float res_2dim_6[2] = {-1.0000, 0.0000};
brp_float x0_2dim_7[2] = {-2.1213, -2.1213};
brp_float res_2dim_7[2] = {-1.0000, -1.0000};
brp_float x0_2dim_8[2] = {-0.0000, -3.0000};
brp_float res_2dim_8[2] = {-0.0000, -1.0000};
brp_float x0_2dim_9[2] = {2.1213, -2.1213};
brp_float res_2dim_9[2] = {0.5000, -0.5000};
brp_float x0_2dim_10[2] = {0.1, 0.1};
brp_float res_2dim_10[2] = {0.1, 0.1};


/* Some inputs/results for DYKSTRA unit tests with dimension = 10 */
const int brp_dim = 10;

brp_float x0_1[10] = {1.0166, 27.8734, -11.6667, -18.5430, -11.4068, -10.9334, -4.3361, -1.6847, -2.1853, 5.4133};
brp_float res_1[10] = {1.0000, 1.0000, -0.0000, -1.0000, -1.0000, -1.0001, -1.0000, -1.0000, -0.0000, 1.0000};

brp_float x0_2[10] = {8.8840, -11.4707, -10.6887, -8.0950, -29.4428, 14.3838, 3.2519, -7.5493, 13.7030, -17.1152};
brp_float res_2[10] = {0.0001, -0.9999, -1.0001, -0.9999, -0.9998, 0.0003, -0.0003, -1.0001, -0.0001, -1.0000};

brp_float x0_3[10] = {-1.0224, -2.4145, 3.1921, 3.1286, -8.6488, -0.3005, -1.6488, 6.2771, 10.9327, 11.0927};
brp_float res_3[10] = {-1.0000, -0.1112, 0.8888, 0.0000, -1.0000, -0.3005, -0.0000, 1.0000, 0.9999, 1.0001};

brp_float x0_4[10] = {-8.6365, 0.7736, -12.1412, -11.1350, -0.0685, 15.3263, -7.6967, 3.7138, -2.2558, 11.1736};
brp_float res_4[10] = {-1.0000, 0.0000, -1.0000, -1.0000, -0.0000, 1.0000, -0.0000, 1.0000, -0.0000, 1.0000};

brp_float x0_5[10] = {-10.8906, 0.3256, 5.5253, 11.0061, 15.4421, 0.8593, -14.9159, -7.4230, -10.6158, 23.5046};
brp_float res_5[10] = {-1.0000, 0.0000, 1.0000, 1.0000, 1.0000, 0.0000, -1.0000, -1.0000, 0.0000, 1.0000};

brp_float x0_99[10] = {1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000};
brp_float res_99[10] = {1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000};

brp_float out[10];

brp_float calc_error(const brp_float * in1, const brp_float * in2, const int len) {
	int i;
	brp_float error = 0.0;
	for (i = 0; i < len; i++)
		error += (in1[i] - in2[i]) * (in1[i] - in2[i]);
	return error;
}

void test_symmetric_dykstra(void) {
	printf("Testing test_symmetric_dykstra()...\n");
	brp_float error;
#ifdef BRP_EMBEDDED
	assert(BRP_DIM == 10);
#else
	BRP_initialize(brp_dim);
#endif

	assert(BRP_A_MAX == 1.0);
	assert(BRP_R_MAX == 1.0);
	// 1
	symmetric_box_rate_projection(x0_1, out);
	error = calc_error(out, res_1, brp_dim);
	assert(error < ABS_ERROR_TOL);
	// 2
	symmetric_box_rate_projection(x0_2, out);
	error = calc_error(out, res_2, brp_dim);
	assert(error < ABS_ERROR_TOL);
	// 3
	symmetric_box_rate_projection(x0_3, out);
	error = calc_error(out, res_3, brp_dim);
	assert(error < ABS_ERROR_TOL);
	// 4
	symmetric_box_rate_projection(x0_4, out);
	error = calc_error(out, res_4, brp_dim);
	assert(error < ABS_ERROR_TOL);
	// 5
	symmetric_box_rate_projection(x0_5, out);
	error = calc_error(out, res_5, brp_dim);
	assert(error < ABS_ERROR_TOL);
#ifndef BRP_EMBEDDED
	BRP_finalize();
#endif
	printf("Test test_symmetric_dykstra() passed.\n");
}

void test_2dim_projection(void) {
	printf("Testing test_2dim_projection()...\n");
	brp_float error;
	brp_float out_2dim[2];
	assert(BRP_A_MAX == 1.0);
	assert(BRP_R_MAX == 1.0);
	// 1
	box_rate_proj_2dim(x0_2dim_1, out_2dim);
	error = calc_error(out_2dim, res_2dim_1, 2);
	assert(error < ABS_ERROR_TOL);
	// 2
	box_rate_proj_2dim(x0_2dim_2, out_2dim);
	error = calc_error(out_2dim, res_2dim_2, 2);
	assert(error < ABS_ERROR_TOL);
	// 3
	box_rate_proj_2dim(x0_2dim_3, out_2dim);
	error = calc_error(out_2dim, res_2dim_3, 2);
	assert(error < ABS_ERROR_TOL);
	// 4
	box_rate_proj_2dim(x0_2dim_4, out_2dim);
	error = calc_error(out_2dim, res_2dim_4, 2);
	assert(error < ABS_ERROR_TOL);
	// 5
	box_rate_proj_2dim(x0_2dim_5, out_2dim);
	error = calc_error(out_2dim, res_2dim_5, 2);
	assert(error < ABS_ERROR_TOL);
	// 6
	box_rate_proj_2dim(x0_2dim_6, out_2dim);
	error = calc_error(out_2dim, res_2dim_6, 2);
	assert(error < ABS_ERROR_TOL);
	// 7
	box_rate_proj_2dim(x0_2dim_7, out_2dim);
	error = calc_error(out_2dim, res_2dim_7, 2);
	assert(error < ABS_ERROR_TOL);
	// 8
	box_rate_proj_2dim(x0_2dim_8, out_2dim);
	error = calc_error(out_2dim, res_2dim_8, 2);
	assert(error < ABS_ERROR_TOL);
	// 9
	box_rate_proj_2dim(x0_2dim_9, out_2dim);
	error = calc_error(out_2dim, res_2dim_9, 2);
	assert(error < ABS_ERROR_TOL);
	// 10
	box_rate_proj_2dim(x0_2dim_10, out_2dim);
	error = calc_error(out_2dim, res_2dim_10, 2);
	assert(error < ABS_ERROR_TOL);
	printf("Test test_2dim_projection() passed.\n");
}

void test_brp_mpc_projection_1(void) {
	printf("Testing test_brp_mpc_projection_1()...\n");
	brp_float error;
	brp_float brp_mpc_out[brp_mpc_horiz_1 * brp_mpc_num_inputs_1];

	// Initialize Projection
	BRP_MPC_initialize(brp_mpc_horiz_1, brp_mpc_num_inputs_1,
			brp_mpc_rate_limits_1, brp_mpc_ampl_limits_1);
	BRP_MPC_initialize_projection(brp_old_input_1);

	// Project
	BRP_MPC_box_rate_projection(brp_mpc_x0_1, brp_mpc_out);

	// Comperate with Matlab result
	error = calc_error(brp_mpc_out, brp_mpc_xp_1,
			brp_mpc_horiz_1 * brp_mpc_num_inputs_1);

	assert(error < ABS_ERROR_TOL);

	BRP_MPC_finalize();

	printf("Test test_brp_mpc_projection_1() passed.\n");
}

void test_brp_mpc_projection_2(void) {
	printf("Testing test_brp_mpc_projection_2()...\n");
	brp_float error;
	brp_float brp_mpc_out[brp_mpc_horiz_2 * brp_mpc_num_inputs_2];

	// Initialize Projection
	BRP_MPC_initialize(brp_mpc_horiz_2, brp_mpc_num_inputs_2,
			brp_mpc_rate_limits_2, brp_mpc_ampl_limits_2);
	BRP_MPC_initialize_projection(brp_old_input_2);

	// Project
	BRP_MPC_box_rate_projection(brp_mpc_x0_2, brp_mpc_out);

	// Comperate with Matlab result
	error = calc_error(brp_mpc_out, brp_mpc_xp_2,
			brp_mpc_horiz_2 * brp_mpc_num_inputs_2);

	assert(error < ABS_ERROR_TOL);

	BRP_MPC_finalize();

	printf("Test test_brp_mpc_projection_2() passed.\n");
}


int main() {

	test_2dim_projection();
	test_symmetric_dykstra();
	test_brp_mpc_projection_1();
	test_brp_mpc_projection_2();

	return 0;
}

