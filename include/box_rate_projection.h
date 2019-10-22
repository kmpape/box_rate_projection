#ifndef SYMMETRIC_BOX_RATE_PROJECTION_H_
#define SYMMETRIC_BOX_RATE_PROJECTION_H_

/* Define BRP_EMBEDDED to disable MALLOC and use a constant dimension BRP_DIM */
//#define BRP_EMBEDDED

#ifdef BRP_EMBEDDED
#define BRP_DIM (10) // use hard-coded problem dimension
#endif

#define BRP_MAX_ITER (100)        // max number of iterations
#define BRP_CHECK_TERMINATION (10) // 0 for no check, k>0 to check every k'th iteration
#define BRP_EPS_ABS (1e-3)
#define BRP_EPS_REL (1e-3)

/* CONSTANTS FOR SYMMETRIC BOX PROJECTION */
#define BRP_A_MAX (1.0f) // max amplitude
#define BRP_R_MAX (1.0f) // max rate


#ifndef BRP_EMBEDDED
void BRP_initialize(const int dim);
void BRP_finalize(void);
#endif
void box_rate_proj_2dim(const float * restrict in, float * restrict out);
void symmetric_box_rate_projection(const float * restrict in, float * restrict out);

#endif /* SYMMETRIC_BOX_RATE_PROJECTION_H_ */
