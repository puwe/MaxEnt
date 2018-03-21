/* Local Maximum Entropy Approximation
 * p(x0, y0) = exp(-(y0-x0)'*beta'*beta*(y0-x0) + lambda*(y0-x0);
 * y0: material point postion
 * x0: nodal position
 * lambda: paramter 
 * beta: lattice basis
 * Z: sum(p(x, y0)) over x0 
 * Using Newton-Raphson to min(lnZ) s.t. y0 = x0*p(x0, y0);
 * Hui Zhou (hzhou@student.ethz.ch);
*/

#ifndef MaxEnt_h
#define MaxEnt_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <assert.h>


void MaxEnt_Y(double *y1, double *p0, const double *x1, const double *y0, const double *lambda, const double *beta0, const int D, const int N, const double *x0);

void MaxEnt_F(double *F, const double *x1, const double *y0, const double *lambda, const double *beta0, const int D, const int N, const double *x0);

void MaxEnt_NR(double *lambda0, const double *y0, const int maxsteps, const double eps, const double *beta0, const int D, const int N, const double *x0);


double MaxEnt_Z(const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0);

void MaxEnt_dlnZ_dL(double *dlnZ_dL, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0);

void MaxEnt_hlnZ_hL(double *hlnZ_hL, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0);

void MaxEnt_dlnZ_dY(double *dlnZ_dY, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0);

void MaxEnt_hlnZ_hY(double *hlnZ_hY, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0);

void MaxEnt_dP_dY(double *dP_dY, const double *y0, const double *lambda, const double *beta0, const int D, const int N, const double *x0);

void mat_inv(double *A, const int dim);

#endif // MaxEnt_h
