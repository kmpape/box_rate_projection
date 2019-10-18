/*
 * symmetric_box_rate_projection.h
 *
 *  Created on: Oct 17, 2019
 *      Author: idris
 */

#ifndef SYMMETRIC_BOX_RATE_PROJECTION_H_
#define SYMMETRIC_BOX_RATE_PROJECTION_H_

/* CONSTANTS FOR DYKSTRA'S ALGORITHM */
#define MAX_ITER (100)        // max number of iterations
#define CHECK_TERMINATION (10) // 0 for no check, k>0 to check every k'th iteration
#define DIM (10)              // hard-coded length of vector to project to avoid use of malloc
#define EPS_ABS (1e-3)
#define EPS_REL (1e-3)

/* CONSTANTS FOR SYMMETRIC BOX PROJECTION */
#define A_MAX (1.0f) // max amplitude
#define R_MAX (1.0f) // max rate

void box_rate_proj_2dim(const float * restrict in, float * restrict out);
void symmetric_box_rate_projection(const float * restrict in, float * restrict out);

#endif /* SYMMETRIC_BOX_RATE_PROJECTION_H_ */
