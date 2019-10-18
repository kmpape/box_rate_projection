#include <stdio.h>
#include "box_rate_projection.h"

/* CONSTANTS FOR DYKSTRA'S ALGORITHM */
#define N_SETS (DIM - 1)
#define FIRST_SET_NUM (0)
#define LAST_SET_NUM (N_SETS - 1)

/* ARRAYS FOR DYKSTRA'S ALGORITHM */
float proj_set_1[DIM];
//float proj_set_2[DIM]; // we use <out> for this
float corr_set_1[DIM];
float corr_set_2[DIM];
float arr_to_proj[DIM];

/* CONSTANTS FOR 2-DIMENSIONAL PROJECTION */
const float RDIAG_MIN = R_MAX - 2 * A_MAX;     // used for point location
const float RDIAG_MAX = 2 * A_MAX - R_MAX;     // used for point location
const float A_MAX_MINUS_R_MAX = A_MAX - R_MAX; // used for fixed-point projection
const float R_MAX_MINUS_A_MAX = R_MAX - A_MAX; // used for fixed-point projection

float min(float in1, float in2) {
	return (in1 > in2) ? in2 : in1;
}

float max(float in1, float in2) {
	return (in1 > in2) ? in1 : in2;
}

float abs_float(float in) {
	return (in > 0) ? in : -in;
}

float inf_norm(const float * in) {
	int i;
	float max_val = abs_float(in[0]);
	for (i = 1; i < DIM; i++)
		max_val = max(abs_float(in[i]), max_val);
	return max_val;
}

float inf_norm_error(const float * in1, const float * in2) {
	int i;
	float max_val = abs_float(in1[0] - in2[0]);
	for (i = 1; i < DIM; i++)
		max_val = max(abs_float(in1[i] - in2[i]), max_val);
	return max_val;
}

void vec_sub(const float * restrict in1, const float * restrict in2, float * restrict out) {
	int i;
	for (i = 0; i < DIM; i++) {
		out[i] = in1[i] - in2[i];
	}
}

void vec_add(const float * restrict in1, const float * restrict in2, float * restrict out) {
	int i;
	for (i = 0; i < DIM; i++) {
		out[i] = in1[i] + in2[i];
	}
}

void vec_copy(const float * restrict in1, float * restrict out) {
	int i;
	for (i = 0; i < DIM; i++) {
		out[i] = in1[i];
	}
}

float * get_set_ptr(float * const in, const int set_num) { // set_num goes from 0 to N_SETS - 1
	return (in + set_num * DIM);
}

/*
 * Projects a 2-dimension vector onto the symmetric rate-amplitude constraint set
 */
void box_rate_proj_2dim(const float * restrict in, float * restrict out) {
	const float rdiag = in[0] + in[1];

	if (rdiag <= RDIAG_MIN) {
		if (in[1] >= R_MAX_MINUS_A_MAX) { // map to corner C1
			out[0] = -A_MAX;
			out[1] = R_MAX_MINUS_A_MAX;
		} else if (in[0] >= R_MAX_MINUS_A_MAX) { // map to corner C4
			out[0] = R_MAX_MINUS_A_MAX;
			out[1] = -A_MAX;
		} else { // project onto lower left corner of the box
			out[0] = max(-A_MAX, in[0]);
			out[1] = max(-A_MAX, in[1]);
		}
	} else if (rdiag >= RDIAG_MAX) {
		if (in[1] <= A_MAX_MINUS_R_MAX) { // map to corner C3
			out[0] = A_MAX;
			out[1] = A_MAX_MINUS_R_MAX;
		} else if (in[0] <= A_MAX_MINUS_R_MAX) { // map to corner C2
			out[0] = A_MAX_MINUS_R_MAX;
			out[1] = A_MAX;
		} else { // project onto upper right corner of the box
			out[0] = min(A_MAX, in[0]);
			out[1] = min(A_MAX, in[1]);
		}
	} else { // project onto diagonal lines defined by rate constraints
		const float ldiag = max(-R_MAX, min(R_MAX, in[0] - in[1])); // project onto [1,-1]-diagonal and truncate
		out[0] = 0.5 * ldiag + 0.5 * rdiag; // reconstruct first coordinate in original basis
		out[1] = -0.5 * ldiag + 0.5 * rdiag; // reconstruct second coordinate in original basis
	}
}

/*
 * This function is used when S = S_0 n S_1 n ... n S_N-2 is split into
 * S_even = S_0 n S_2 n ... AND S_odd = S_1 n S_3 n ...
 * where the projections for S_even and S_odd are known. To project onto S_even,
 * pass start_index = 0. To project onto S_odd, pass start_index = 1.
 */
void box_rate_proj_odd_even(const float * restrict in, float * restrict out, const int start_index) {
	int i_set;
	for (i_set = start_index; i_set < N_SETS; i_set+=2) {
		box_rate_proj_2dim(in + i_set, out + i_set); // modifies out[i_set] and out[i_set + 1]
	}
	// Need to copy elements not modified inside the loop (these are left unchanged)
	if (start_index == 1) {
		out[0] = in[0];
		if (DIM % 2 == 0) { // DIM is odd
			out[DIM - 1] = in[DIM - 1];
		}
	} else { // start_index = 0
		if (DIM % 2 != 0) { // DIM is even
			out[DIM - 1] = in[DIM - 1];
		}
	}
}

/* Dykstra's algorithm for two sets */
void symmetric_box_rate_projection(const float * restrict in, float * restrict out) {
	int i_iter;
	float abs_error, last_iter_inf_norm;
	// NOTE: we use out for proj_set_2

	/* First iteration */
	box_rate_proj_odd_even(in, proj_set_1, 1);
	vec_sub(in, proj_set_1, corr_set_1);

	box_rate_proj_odd_even(proj_set_1, out, 0);
	vec_sub(proj_set_1, out, corr_set_2);

	for (i_iter = 1; i_iter < MAX_ITER; i_iter++) {
		/* Set 1 */
		vec_add(out, corr_set_1, arr_to_proj);
		box_rate_proj_odd_even(arr_to_proj, proj_set_1, 1);
		vec_sub(arr_to_proj, proj_set_1, corr_set_1);

		/* Set 2 */
		vec_add(proj_set_1, corr_set_2, arr_to_proj);
		box_rate_proj_odd_even(arr_to_proj, out, 0);
		vec_sub(arr_to_proj, out, corr_set_2);

		/* Check for termination */
		if ((CHECK_TERMINATION) && (i_iter % CHECK_TERMINATION == 0)) {
			abs_error = inf_norm_error(out, proj_set_1);
			if (abs_error == 0) {
				break;
			} else {
				last_iter_inf_norm = inf_norm(out);
				if (last_iter_inf_norm == 0) {
					break;
				} else {
					if ((abs_error < EPS_ABS) && (abs_error < EPS_REL * last_iter_inf_norm))
						break;
				}
			}
		}
	}
}


#if 0 // old version
float proj_per_set[DIM * N_SETS];
float diff_per_set[DIM * N_SETS];

void box_rate_proj_2dim_inplace(float * in_out) {
	const float rdiag = in_out[0] + in_out[1];

	if (rdiag <= RDIAG_MIN) {
		if (in_out[1] >= R_MAX_MINUS_A_MAX) { // map to corner C1
			in_out[0] = -A_MAX;
			in_out[1] = R_MAX_MINUS_A_MAX;
		} else if (in_out[0] >= R_MAX_MINUS_A_MAX) { // map to corner C4
			in_out[0] = R_MAX_MINUS_A_MAX;
			in_out[1] = -A_MAX;
		} else { // project onto lower left corner of the box
			in_out[0] = max(-A_MAX, in_out[0]);
			in_out[1] = max(-A_MAX, in_out[1]);
		}
	} else if (rdiag >= RDIAG_MAX) {
		if (in_out[1] <= A_MAX_MINUS_R_MAX) { // map to corner C3
			in_out[0] = A_MAX;
			in_out[1] = A_MAX_MINUS_R_MAX;
		} else if (in_out[0] <= A_MAX_MINUS_R_MAX) { // map to corner C2
			in_out[0] = A_MAX_MINUS_R_MAX;
			in_out[1] = A_MAX;
		} else { // project onto upper right corner of the box
			in_out[0] = min(A_MAX, in_out[0]);
			in_out[1] = min(A_MAX, in_out[1]);
		}
	} else { // project onto diagonal lines defined by rate constraints
		const float ldiag = max(-R_MAX, min(R_MAX, in_out[0] - in_out[1])); // project onto [1,-1]-diagonal and truncate
		in_out[0] = 0.5 * ldiag + 0.5 * rdiag; // reconstruct first coordinate in original basis
		in_out[1] = -0.5 * ldiag + 0.5 * rdiag; // reconstruct second coordinate in original basis
	}
}

/* Dykstra's algorithm */
void symmetric_box_rate_projection(const float * restrict in, float * restrict out) {
	int i_iter, i_set;

	/* Initialize algorithm */
	vec_copy(in, get_set_ptr(proj_per_set, LAST_SET_NUM));

	//NOTE: First iteration done outside of the loop in order to avoid resetting diff_per_set to zero

	/*** First Iteration ***/
	/* First Set */
	vec_copy(in, get_set_ptr(proj_per_set, 0));
	box_rate_proj_2dim(in, proj_per_set);
	vec_sub(proj_per_set, in, diff_per_set);
	/* Other Sets */
	for (i_set = 1; i_set < N_SETS; i_set++) {
		vec_copy(get_set_ptr(proj_per_set, i_set - 1), get_set_ptr(proj_per_set, i_set));
		box_rate_proj_2dim(get_set_ptr(proj_per_set, i_set - 1) + i_set,
						   get_set_ptr(proj_per_set, i_set) + i_set);
		vec_sub(get_set_ptr(proj_per_set, i_set),
				get_set_ptr(proj_per_set, i_set - 1),
				get_set_ptr(diff_per_set, i_set));
	}

	/*** Other Iterations ***/
	for (i_iter = 1; i_iter < MAX_ITER; i_iter++) {
		/* First Set */
		vec_copy(get_set_ptr(proj_per_set, LAST_SET_NUM), get_set_ptr(proj_per_set, 0));
		vec_sub(get_set_ptr(proj_per_set, LAST_SET_NUM), diff_per_set, arr_to_proj);
		box_rate_proj_2dim(arr_to_proj, proj_per_set);
		vec_sub(proj_per_set, arr_to_proj, diff_per_set);

		/* Other Sets */
		for (i_set = 1; i_set < N_SETS; i_set++) {
			vec_copy(get_set_ptr(proj_per_set, i_set - 1), get_set_ptr(proj_per_set, i_set));
			vec_sub(get_set_ptr(proj_per_set, i_set - 1), get_set_ptr(diff_per_set, i_set), arr_to_proj);
			box_rate_proj_2dim(arr_to_proj + i_set, get_set_ptr(proj_per_set, i_set) + i_set);
			vec_sub(get_set_ptr(proj_per_set, i_set), arr_to_proj, get_set_ptr(diff_per_set, i_set));
		}

		/* Check for termination */
		if ((CHECK_TERMINATION) && (i_iter % CHECK_TERMINATION == 0)) {
			break; // TODO
		}
	}

	/*** Copy Result ***/
	vec_copy(get_set_ptr(proj_per_set, LAST_SET_NUM), out);
}
#endif
