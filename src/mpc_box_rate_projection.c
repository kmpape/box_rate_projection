#include <stdio.h>
#include "box_rate_projection.h"

#ifndef BRP_EMBEDDED
#include <stdlib.h> // for malloc
#include <assert.h>
brp_float * BRP_MPC_proj_set_1;
brp_float * BRP_MPC_corr_set_1;
brp_float * BRP_MPC_corr_set_2;
brp_float * BRP_MPC_arr_to_proj;
int BRP_MPC_HORIZ;
int BRP_MPC_DIM;
int BRP_MPC_NUM_INPUTS;
int BRP_MPC_is_initialized = 0;

// Limits for symmetric projection
brp_float *BRP_MPC_R_MAX, *BRP_MPC_A_MAX;
brp_float *BRP_MPC_RDIAG_MIN, *BRP_MPC_RDIAG_MAX;
brp_float *BRP_MPC_A_MAX_MINUS_R_MAX, *BRP_MPC_R_MAX_MINUS_A_MAX;

// Limits and corner points for non-symmetric projection
//brp_float *BRP_MPC_SIGN_U0;
brp_float *BRP_MPC_A0_MAX, *BRP_MPC_A0_MIN, *BRP_MPC_A1_MAX, *BRP_MPC_A1_MIN; // max(-a, -r + u_old)
brp_float *BRP_MPC_C1, *BRP_MPC_C2, *BRP_MPC_C3, *BRP_MPC_C4;
brp_float *BRP_MPC_DIAG_LIM1, *BRP_MPC_DIAG_LIM2;

void BRP_MPC_initialize(const int horiz, const int num_inputs, const brp_float * rate_limits,
				    const brp_float * ampl_limits) {
	int i;

	BRP_MPC_HORIZ = horiz;
	BRP_MPC_NUM_INPUTS = num_inputs;
	BRP_MPC_DIM = BRP_MPC_HORIZ * BRP_MPC_NUM_INPUTS;

	// Allocate Arrays
	BRP_MPC_proj_set_1 = (brp_float *)malloc(BRP_MPC_HORIZ * sizeof(brp_float));
	assert(BRP_MPC_proj_set_1);
	BRP_MPC_corr_set_1 = (brp_float *)malloc(BRP_MPC_HORIZ * sizeof(brp_float));
	if (!BRP_MPC_corr_set_1) {
		free(BRP_MPC_proj_set_1);
		assert(BRP_MPC_corr_set_1);
	}
	BRP_MPC_corr_set_2 = (brp_float *)malloc(BRP_MPC_HORIZ * sizeof(brp_float));
	if (!BRP_MPC_corr_set_2) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1);
		assert(BRP_MPC_corr_set_2);
	}
	BRP_MPC_arr_to_proj = (brp_float *)malloc(BRP_MPC_HORIZ * sizeof(brp_float));
	if (!BRP_MPC_arr_to_proj) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		assert(BRP_MPC_arr_to_proj);
	}
	BRP_MPC_RDIAG_MIN = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_RDIAG_MIN) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj);
		assert(BRP_MPC_RDIAG_MIN);
	}
	BRP_MPC_RDIAG_MAX = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_RDIAG_MAX) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN);
		assert(BRP_MPC_RDIAG_MAX);
	}
	BRP_MPC_A_MAX_MINUS_R_MAX = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_A_MAX_MINUS_R_MAX) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		assert(BRP_MPC_A_MAX_MINUS_R_MAX);
	}
	BRP_MPC_R_MAX_MINUS_A_MAX = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_R_MAX_MINUS_A_MAX) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX);
		assert(BRP_MPC_R_MAX_MINUS_A_MAX);
	}
	BRP_MPC_A0_MAX = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_A0_MAX) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX);
		assert(BRP_MPC_A0_MAX);
	}
	BRP_MPC_A0_MIN = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_A0_MIN) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		assert(BRP_MPC_A0_MIN);
	}
	BRP_MPC_C1 = (brp_float *)malloc(2 * BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_C1) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN);
		assert(BRP_MPC_C1);
	}
	BRP_MPC_C2 = (brp_float *)malloc(2 * BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_C2) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1);
		assert(BRP_MPC_C2);
	}
	BRP_MPC_C3 = (brp_float *)malloc(2 * BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_C3) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2);
		assert(BRP_MPC_C3);
	}
	BRP_MPC_C4 = (brp_float *)malloc(2 * BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_C4) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
		assert(BRP_MPC_C4);
	}
	BRP_MPC_DIAG_LIM1 = (brp_float *)malloc(2 * BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_DIAG_LIM1) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
		free(BRP_MPC_C4);
		assert(BRP_MPC_DIAG_LIM1);
	}
	BRP_MPC_DIAG_LIM2 = (brp_float *)malloc(2 * BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_DIAG_LIM2) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
		free(BRP_MPC_C4); free(BRP_MPC_DIAG_LIM1);
		assert(BRP_MPC_DIAG_LIM2);
	}
	BRP_MPC_A_MAX = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_A_MAX) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
		free(BRP_MPC_C4); free(BRP_MPC_DIAG_LIM1); free(BRP_MPC_DIAG_LIM2);
		assert(BRP_MPC_A_MAX);
	}
	BRP_MPC_R_MAX = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_R_MAX) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
		free(BRP_MPC_C4); free(BRP_MPC_DIAG_LIM1); free(BRP_MPC_DIAG_LIM2); free(BRP_MPC_A_MAX);
		assert(BRP_MPC_R_MAX);
	}
	BRP_MPC_A1_MAX = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_A1_MAX) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
		free(BRP_MPC_C4); free(BRP_MPC_DIAG_LIM1); free(BRP_MPC_DIAG_LIM2); free(BRP_MPC_A_MAX);
		free(BRP_MPC_R_MAX);
		assert(BRP_MPC_A1_MAX);
	}
	BRP_MPC_A1_MIN = (brp_float *)malloc(BRP_MPC_NUM_INPUTS * sizeof(brp_float));
	if (!BRP_MPC_A1_MIN) {
		free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
		free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
		free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
		free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
		free(BRP_MPC_C4); free(BRP_MPC_DIAG_LIM1); free(BRP_MPC_DIAG_LIM2); free(BRP_MPC_A_MAX);
		free(BRP_MPC_R_MAX); free(BRP_MPC_A1_MAX);
		assert(BRP_MPC_A1_MIN);
	}

	// Assign provided limits and constant parameters
	for (i = 0; i < BRP_MPC_NUM_INPUTS; i++) {
		BRP_MPC_A_MAX[i] = ampl_limits[i];
		BRP_MPC_R_MAX[i] = rate_limits[i];
		BRP_MPC_RDIAG_MIN[i] = rate_limits[i] - 2 * ampl_limits[i];
		BRP_MPC_RDIAG_MAX[i] = -BRP_MPC_RDIAG_MIN[i];
		BRP_MPC_A_MAX_MINUS_R_MAX[i] = ampl_limits[i] - rate_limits[i];
		BRP_MPC_R_MAX_MINUS_A_MAX[i] = -BRP_MPC_A_MAX_MINUS_R_MAX[i];
	}

	BRP_MPC_is_initialized = 1;
}

void BRP_MPC_finalize(void) {
	BRP_MPC_is_initialized = 0;
	free(BRP_MPC_proj_set_1); free(BRP_MPC_corr_set_1); free(BRP_MPC_corr_set_2);
	free(BRP_MPC_arr_to_proj); free(BRP_MPC_RDIAG_MIN); free(BRP_MPC_RDIAG_MAX);
	free(BRP_MPC_A_MAX_MINUS_R_MAX); free(BRP_MPC_R_MAX_MINUS_A_MAX); free(BRP_MPC_A0_MAX);
	free(BRP_MPC_A0_MIN); free(BRP_MPC_C1); free(BRP_MPC_C2); free(BRP_MPC_C3);
	free(BRP_MPC_C4); free(BRP_MPC_DIAG_LIM1); free(BRP_MPC_DIAG_LIM2);
	free(BRP_MPC_A_MAX); free(BRP_MPC_R_MAX); free(BRP_MPC_A1_MIN); free(BRP_MPC_A1_MAX);
}

void BRP_MPC_initialize_projection(const brp_float * old_input) {
	int i;
	assert(BRP_MPC_is_initialized);

	for (i = 0; i < BRP_MPC_NUM_INPUTS; i++) {
		// Compute new u0 limits
		// a0_min = max(-a, -r + u_old);
		BRP_MPC_A0_MIN[i] = BRP_max(-BRP_MPC_A_MAX[i], old_input[i] - BRP_MPC_R_MAX[i]);
		// a0_max = min(a, r + u_old)
		BRP_MPC_A0_MAX[i] = BRP_min(BRP_MPC_A_MAX[i], old_input[i] + BRP_MPC_R_MAX[i]);
		// a1_min = max(-a, -r + a0_min)
		BRP_MPC_A1_MIN[i] = BRP_max(-BRP_MPC_A_MAX[i], BRP_MPC_A0_MIN[i] - BRP_MPC_R_MAX[i]);
		// a1_max = min(a, r + obj.a0_max)
		BRP_MPC_A1_MAX[i] = BRP_min(BRP_MPC_A_MAX[i], BRP_MPC_A0_MAX[i] + BRP_MPC_R_MAX[i]);
		// Compute corner points
		// C1 = [a0_min; min(r + a0_min, a)]
		BRP_MPC_C1[2 * i] = BRP_MPC_A0_MIN[i];
		BRP_MPC_C1[2 * i + 1] = BRP_min(BRP_MPC_R_MAX[i] + BRP_MPC_A0_MIN[i], BRP_MPC_A_MAX[i]);
		// C2 = [min(a0_max, a - r); min(a, r + a0_max)]
		BRP_MPC_C2[2 * i] = BRP_min(BRP_MPC_A0_MAX[i], BRP_MPC_A_MAX_MINUS_R_MAX[i]);
		BRP_MPC_C2[2 * i + 1] = BRP_min(BRP_MPC_A_MAX[i], BRP_MPC_R_MAX[i] + BRP_MPC_A0_MAX[i]);
		// C3 = [a0_max; -r + a0_max]
		BRP_MPC_C3[2 * i] = BRP_MPC_A0_MAX[i];
		BRP_MPC_C3[2 * i + 1] = BRP_MPC_A0_MAX[i] - BRP_MPC_R_MAX[i];
		// C4 = [max(r - a, a0_min); max(-a, -r + a0_min)]
		BRP_MPC_C4[2 * i] = BRP_max(BRP_MPC_R_MAX_MINUS_A_MAX[i], BRP_MPC_A0_MIN[i]);
		BRP_MPC_C4[2 * i + 1] = BRP_max(-BRP_MPC_A_MAX[i], BRP_MPC_A0_MIN[i]-BRP_MPC_R_MAX[i]);

		// Compute limits of diagonal for regions A1 and A2
		BRP_MPC_DIAG_LIM1[2 * i] = BRP_MPC_C1[2 * i] + BRP_MPC_C1[2 * i + 1];
		BRP_MPC_DIAG_LIM1[2 * i + 1] = BRP_MPC_C2[2 * i] + BRP_MPC_C2[2 * i + 1];
		BRP_MPC_DIAG_LIM2[2 * i] = BRP_MPC_C4[2 * i] + BRP_MPC_C4[2 * i + 1];
		BRP_MPC_DIAG_LIM2[2 * i + 1] = BRP_MPC_C3[2 * i] + BRP_MPC_C3[2 * i + 1];
	}
}

/*
 * Projects a 2-dimensional vector the non-symmetric first input set
 */
void BRP_MPC_box_rate_proj_2dim(const brp_float * restrict in, brp_float * restrict out, const int i_input) {
	const brp_float rdiag = in[0] + in[1];

	if (rdiag <= BRP_MPC_RDIAG_MIN[i_input]) {
		if (in[1] >= BRP_MPC_R_MAX_MINUS_A_MAX[i_input]) { // map to corner C1
			out[0] = -BRP_MPC_A_MAX[i_input];
			out[1] = BRP_MPC_R_MAX_MINUS_A_MAX[i_input];
		} else if (in[0] >= BRP_MPC_R_MAX_MINUS_A_MAX[i_input]) { // map to corner C4
			out[0] = BRP_MPC_R_MAX_MINUS_A_MAX[i_input];
			out[1] = -BRP_MPC_A_MAX[i_input];
		} else { // project onto lower left corner of the box
			out[0] = BRP_max(-BRP_MPC_A_MAX[i_input], in[0]);
			out[1] = BRP_max(-BRP_MPC_A_MAX[i_input], in[1]);
		}
	} else if (rdiag >= BRP_MPC_RDIAG_MAX[i_input]) {
		if (in[1] <= BRP_MPC_A_MAX_MINUS_R_MAX[i_input]) { // map to corner C3
			out[0] = BRP_MPC_A_MAX[i_input];
			out[1] = BRP_MPC_A_MAX_MINUS_R_MAX[i_input];
		} else if (in[0] <= BRP_MPC_A_MAX_MINUS_R_MAX[i_input]) { // map to corner C2
			out[0] = BRP_MPC_A_MAX_MINUS_R_MAX[i_input];
			out[1] = BRP_MPC_A_MAX[i_input];
		} else { // project onto upper right corner of the box
			out[0] = BRP_min(BRP_MPC_A_MAX[i_input], in[0]);
			out[1] = BRP_min(BRP_MPC_A_MAX[i_input], in[1]);
		}
	} else { // project onto diagonal lines defined by rate constraints
		const brp_float ldiag = BRP_max(-BRP_MPC_R_MAX[i_input], BRP_min(BRP_MPC_R_MAX[i_input], in[0] - in[1])); // project onto [1,-1]-diagonal and truncate
		out[0] = 0.5 * ldiag + 0.5 * rdiag; // reconstruct first coordinate in original basis
		out[1] = -0.5 * ldiag + 0.5 * rdiag; // reconstruct second coordinate in original basis
	}
}

/*
 * Projects a 2-dimensional vector onto the symmetric rate-amplitude constraint set
 */
void BRP_MPC_box_rate_proj_2dim_nonsym(const brp_float * restrict in, brp_float * restrict out, const int i_input) {
	// [1; 1]' * x0
	const brp_float rdiag = in[0] + in[1];
	// [-1; 1] * x0
	const brp_float ldiag = in[1] - in[0];

	const int ix = 2 * i_input;
	const int iy = 2 * i_input + 1;

	if ((rdiag >= BRP_MPC_DIAG_LIM1[ix]) && (rdiag <= BRP_MPC_DIAG_LIM1[iy]) && (ldiag >= BRP_MPC_R_MAX[i_input])) { // A1
		// pl = min(r, pl);
        // xp = 0.5 * [-1; 1] * pl + 0.5 * [1; 1] * pr;
		//if (i_input == 0)
		//	printf("A1\n");
		const brp_float tmp = BRP_min(BRP_MPC_R_MAX[i_input], ldiag);
		out[0] = 0.5 * (rdiag - tmp);
		out[1] = 0.5 * (tmp + rdiag);
	} else if ((rdiag >= BRP_MPC_DIAG_LIM2[ix]) && (rdiag <= BRP_MPC_DIAG_LIM2[iy]) && (ldiag <= -BRP_MPC_R_MAX[i_input])) { // A2
        // pl = max(-r, pl);
        // xp = 0.5 * [-1; 1] * pl + 0.5 * [1; 1] * pr;
		//if (i_input == 0)
		//	printf("A2\n");
		const brp_float tmp = BRP_max(-BRP_MPC_R_MAX[i_input], ldiag);
		out[0] = 0.5 * (-tmp + rdiag);
		out[1] = 0.5 * (tmp + rdiag);
	} else if ((rdiag <= BRP_MPC_DIAG_LIM1[ix]) && (in[0] <= BRP_MPC_C1[ix]) && (in[1] >= BRP_MPC_C1[iy])) { // C1
		//if (i_input == 0)
		//	printf("C1\n");
		// ((pr <= not_norm_p1(1)) && (x0(1) <= C1(1)) && (x0(2) >= C1(2))) % C1
		out[0] = BRP_MPC_C1[ix];
		out[1] = BRP_MPC_C1[iy];
	} else if ((rdiag <= BRP_MPC_DIAG_LIM1[ix]) && (in[0] <= BRP_MPC_C4[ix]) && (in[1] <= BRP_MPC_C4[iy])) { // C4
		//if (i_input == 0)
		//	printf("C4\n");
		// ((pr <= not_norm_p2(1)) && (x0(1) >= C4(1)) && (x0(2) <= C4(2))) % C4
		out[0] = BRP_MPC_C4[ix];
		out[1] = BRP_MPC_C4[iy];
	} else if ((rdiag >= BRP_MPC_DIAG_LIM1[iy]) && (in[0] <= BRP_MPC_C2[ix]) && (in[1] >= BRP_MPC_C2[iy])) { // C2
		//if (i_input == 0)
		//	printf("C2\n");
		// ((pr >= not_norm_p1(2)) && (x0(1) <= C2(1)) && (x0(2) >= C2(2))) % C2
		out[0] = BRP_MPC_C2[ix];
		out[1] = BRP_MPC_C2[iy];
	} else if ((rdiag >= BRP_MPC_DIAG_LIM2[iy]) && (in[0] >= BRP_MPC_C3[ix]) && (in[1] <= BRP_MPC_C3[iy])) { // C3
		//if (i_input == 0)
		//	printf("C3\n");
		// ((pr >= not_norm_p2(2)) && (x0(1) >= C3(1)) && (x0(2) <= C3(2))) % C3
		out[0] = BRP_MPC_C3[ix];
		out[1] = BRP_MPC_C3[iy];
	} else {
		//if (i_input == 0)
		//	printf("B\n");
		// [max(a0_min, min(x0(1), a0_max))
		//  max(a1_min, min(x0(2), a1_max))]
		out[0] = BRP_max(BRP_MPC_A0_MIN[i_input], BRP_min(in[0], BRP_MPC_A0_MAX[i_input]));
		out[1] = BRP_max(BRP_MPC_A1_MIN[i_input], BRP_min(in[1], BRP_MPC_A1_MAX[i_input]));
	}
}

/*
 * This function is used when S = S_0 n S_1 n ... n S_N-2 is split into
 * S_even = S_0 n S_2 n ... AND S_odd = S_1 n S_3 n ...
 * where the projections for S_even and S_odd are known. To project onto S_even,
 * pass start_index = 0. To project onto S_odd, pass start_index = 1.
 */
void BRP_MPC_box_rate_proj_odd_even(const brp_float * restrict in, brp_float * restrict out, const int start_index,
									const int dim, const int i_input) {
	int i_set;
	for (i_set = start_index; i_set < (dim - 1); i_set+=2) {
		BRP_MPC_box_rate_proj_2dim(in + i_set, out + i_set, i_input); // modifies out[i_set] and out[i_set + 1]
	}
	// Need to copy elements not modified inside the loop (these are left unchanged)
	if (start_index == 1) {
		out[0] = in[0];
		if (dim % 2 == 0) { // BRP_MPC_DIM is odd
			out[dim - 1] = in[dim - 1];
		}
	} else { // start_index = 0
		if (BRP_MPC_DIM % 2 != 0) { // BRP_MPC_DIM is even
			out[dim - 1] = in[dim - 1];
		}
	}
}

void BRP_MPC_box_rate_proj_start1(const brp_float * restrict in, brp_float * restrict out, const int i_input) {
	int i_set;
	for (i_set = 1; i_set < (BRP_MPC_HORIZ - 1); i_set+=2) {
		BRP_MPC_box_rate_proj_2dim(in + i_set, out + i_set, i_input); // modifies out[i_set] and out[i_set + 1]
	}
	// Need to copy elements not modified inside the loop (these are left unchanged)
	out[0] = in[0];
	if (BRP_MPC_HORIZ % 2 == 0) { // BRP_MPC_DIM is odd
		out[BRP_MPC_HORIZ - 1] = in[BRP_MPC_HORIZ - 1];
	}
}

void BRP_MPC_box_rate_proj_start0(const brp_float * restrict in, brp_float * restrict out, const int i_input) {
	int i_set;
	// First set is non-symmetric due to old_input
	BRP_MPC_box_rate_proj_2dim_nonsym(in, out, i_input);
	for (i_set = 2; i_set < (BRP_MPC_HORIZ - 1); i_set += 2) {
		BRP_MPC_box_rate_proj_2dim(in + i_set, out + i_set, i_input); // modifies out[i_set] and out[i_set + 1]
	}
	// Need to copy elements not modified inside the loop (these are left unchanged)
	if (BRP_MPC_HORIZ % 2 != 0) { // BRP_MPC_DIM is even
		out[BRP_MPC_HORIZ - 1] = in[BRP_MPC_HORIZ - 1];
	}
}

/* Dykstra's algorithm for two sets */
void BRP_MPC_box_rate_projection(const brp_float * restrict in, brp_float * restrict out) {
	int i_iter, i_input;
	brp_float abs_error, last_iter_inf_norm;
	// NOTE: we use out for proj_set_2

	assert(BRP_MPC_is_initialized == 1);
	for (i_input = 0; i_input < BRP_MPC_NUM_INPUTS; i_input++,
		 in+=BRP_MPC_HORIZ, out+=BRP_MPC_HORIZ) {
		/* First iteration */
		BRP_MPC_box_rate_proj_start1(in, BRP_MPC_proj_set_1, i_input);
		BRP_vec_sub(in, BRP_MPC_proj_set_1, BRP_MPC_corr_set_1, BRP_MPC_HORIZ);

		BRP_MPC_box_rate_proj_start0(BRP_MPC_proj_set_1, out, i_input);
		BRP_vec_sub(BRP_MPC_proj_set_1, out, BRP_MPC_corr_set_2, BRP_MPC_HORIZ);

		for (i_iter = 1; i_iter < BRP_MAX_ITER; i_iter++) {
			/* Set 1 */
			BRP_vec_add(out, BRP_MPC_corr_set_1, BRP_MPC_arr_to_proj, BRP_MPC_HORIZ);
			BRP_MPC_box_rate_proj_start1(BRP_MPC_arr_to_proj, BRP_MPC_proj_set_1, i_input);
			BRP_vec_sub(BRP_MPC_arr_to_proj, BRP_MPC_proj_set_1, BRP_MPC_corr_set_1, BRP_MPC_HORIZ);

			/* Set 2 */
			BRP_vec_add(BRP_MPC_proj_set_1, BRP_MPC_corr_set_2, BRP_MPC_arr_to_proj, BRP_MPC_HORIZ);
			BRP_MPC_box_rate_proj_start0(BRP_MPC_arr_to_proj, out, i_input);
			BRP_vec_sub(BRP_MPC_arr_to_proj, out, BRP_MPC_corr_set_2, BRP_MPC_HORIZ);

			/* Check for termination: always check first iterate */
			if ((i_iter == 1) || ((BRP_CHECK_TERMINATION) && (i_iter % BRP_CHECK_TERMINATION == 0))) {
				//abs_error = BRP_inf_norm_error(out, BRP_MPC_proj_set_1, BRP_MPC_DIM);
				abs_error = BRP_inf_norm_error(out, BRP_MPC_proj_set_1, BRP_MPC_HORIZ);
				if (abs_error == 0) {
					break;
				} else {
					last_iter_inf_norm = BRP_inf_norm(out, BRP_MPC_HORIZ);
					if (last_iter_inf_norm == 0) {
						break;
					} else {
						if ((abs_error < BRP_EPS_ABS) && (abs_error < BRP_EPS_REL * last_iter_inf_norm))
							break;
					}
				}
			}
		}
	}
}
#endif
