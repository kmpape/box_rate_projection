#ifndef SYMMETRIC_BOX_RATE_PROJECTION_H_
#define SYMMETRIC_BOX_RATE_PROJECTION_H_

/* Define BRP_EMBEDDED to disable MALLOC and use a constant dimension BRP_DIM */
//#define BRP_EMBEDDED

//#define DOUBLE_PRECISION
#define DOUBLE_PRECISION
#ifdef DOUBLE_PRECISION
typedef double brp_float;
#else
typedef float brp_float;
#endif

#ifdef BRP_EMBEDDED
#define BRP_DIM (10) // use hard-coded problem dimension
#endif

#define BRP_MAX_ITER (100)        // max number of iterations
#define BRP_CHECK_TERMINATION (0) // 0 for no check, k>0 to check every k'th iteration
#define BRP_EPS_ABS (1e-3)
#define BRP_EPS_REL (1e-3)

/* CONSTANTS FOR SYMMETRIC BOX PROJECTION */
#ifdef DOUBLE_PRECISION
#define BRP_A_MAX (1.0) // max amplitude
#define BRP_R_MAX (1.0) // max rate
#else
#define BRP_A_MAX (1.0f) // max amplitude
#define BRP_R_MAX (1.0f) // max rate
#endif


#ifndef BRP_EMBEDDED
void BRP_initialize(const int dim);
void BRP_finalize(void);

void BRP_MPC_initialize(const int horiz, const int num_inputs, const brp_float * rate_limits, const brp_float * ampl_limits);
void BRP_MPC_finalize(void);
void BRP_MPC_initialize_projection(const brp_float * old_input);
void BRP_MPC_box_rate_projection(const brp_float * restrict in, brp_float * restrict out);
#endif

#ifdef __cplusplus
extern "C"
{
#endif
void box_rate_proj_2dim(const brp_float * restrict in, brp_float * restrict out); // uses hard-coded limits
void box_rate_proj_2dim_w_limits(const brp_float * restrict in, brp_float * restrict out, const brp_float r_max, const brp_float a_max); // uses provided limits
void symmetric_box_rate_projection(const brp_float * restrict in, brp_float * restrict out);
#ifdef __cplusplus
}
#endif

/* Utils */
brp_float BRP_min(brp_float in1, brp_float in2);
brp_float BRP_max(brp_float in1, brp_float in2);
brp_float BRP_abs_float(brp_float in);
brp_float BRP_inf_norm(const brp_float * in, const int len);
brp_float BRP_inf_norm_error(const brp_float * in1, const brp_float * in2, const int len);
void BRP_vec_sub(const brp_float * restrict in1, const brp_float * restrict in2, brp_float * restrict out, const int len);
void BRP_vec_add(const brp_float * restrict in1, const brp_float * restrict in2, brp_float * restrict out, const int len);
void BRP_vec_copy(const brp_float * restrict in1, brp_float * restrict out, const int len);

#endif /* SYMMETRIC_BOX_RATE_PROJECTION_H_ */
