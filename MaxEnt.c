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
#include "MaxEnt.h"

void MaxEnt_Y(double *y1, double *p0, const double *x1, const double *y0, const double *lambda, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double p = 0;
	double z = 0;
	for(int i = 0; i < N; ++i){
		double ds2 = 0;
		double ds1 = 0;
		for(short j = 0; j < D; ++j){
			double hj = 0;
			for(short k = 0; k < D; ++k){
				hj += beta0[i*D2+j*D+k]*(y0[k] - x0[D*i+k]);
			}
			//hj = fabs(hj-nearbyint(hj));
			ds1 += lambda[j]*(y0[j] - x0[D*i+j]);
			ds2 += hj*hj;
		}
		p = exp(-ds2 + ds1);
		z += p; 
		p0[i] = p;
	}
	for(int i = 0; i < N; ++i){
		p0[i] /= z;	
	}

	for(short j = 0; j < D; ++j){
		y1[j] = 0;
	}
	for(int i = 0; i < N; ++i){
		for(short j = 0; j < D; ++j){
			y1[j] += p0[i]*x1[D*i+j];		
		}
	}
}

double MaxEnt_Z(const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double p = 0;
	double z = 0;
	for(int i = 0; i < N; ++i){
		double ds2 = 0;
		double ds1 = 0;
		for(short j = 0; j < D; ++j){
			double hj = 0;
			for(short k = 0; k < D; ++k){
				hj += beta0[i*D2+j*D+k]*(y0[k] - x0[D*i+k]);
			}
			//hj = fabs(hj-nearbyint(hj));
			ds1 += lambda[j]*(y0[j] - x0[D*i+j]);
			ds2 += hj*hj;
		}
		p = exp(-ds2 + ds1);
		z += p; 
	}
	return z;
}

void MaxEnt_dlnZ_dL(double *dlnZ_dL, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double p = 0;
	double z = MaxEnt_Z(lambda, y0, beta0, D, N, x0);

	for(int i = 0; i < D; ++i){
		dlnZ_dL[i] = 0;
	}

	for(int i = 0; i < N; ++i){
		double ds2 = 0;
		double ds1 = 0;
		for(short j = 0; j < D; ++j){
			double hj = 0;
			for(short k = 0; k < D; ++k){
				hj += beta0[i*D2+j*D+k]*(y0[k] - x0[D*i+k]);
			}
			//hj = fabs(hj-nearbyint(hj));
			ds1 += lambda[j]*(y0[j] - x0[D*i+j]);
			ds2 += hj*hj;
		}
		p = exp(-ds2 + ds1);
		for(short j = 0; j < D; ++j){
			dlnZ_dL[j] += (y0[j] - x0[D*i+j])*p/z;
		}
	}

}

void MaxEnt_hlnZ_hL(double *hlnZ_hL, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double p = 0;
	double *dlnZ_dL = (double *) malloc(D*sizeof(double)); assert(dlnZ_dL!=NULL);

	double z = MaxEnt_Z(lambda, y0, beta0, D, N, x0);
	MaxEnt_dlnZ_dL(dlnZ_dL, lambda, y0, beta0, D, N, x0);

	for(int i = 0; i < D2; ++i){
		hlnZ_hL[i] = 0;
	}

	for(int i = 0; i < N; ++i){
		double ds2 = 0;
		double ds1 = 0;
		for(short j = 0; j < D; ++j){
			double hj = 0;
			for(short k = 0; k < D; ++k){
				hj += beta0[i*D2+j*D+k]*(y0[k] - x0[D*i+k]);
			}
			//hj = fabs(hj-nearbyint(hj));
			ds1 += lambda[j]*(y0[j] - x0[D*i+j]);
			ds2 += hj*hj;
		}
		p = exp(-ds2 + ds1);
		for(short j = 0; j < D; ++j){
			for(short k = 0; k < D; ++k){
				hlnZ_hL[D*j+k] += (y0[j] - x0[D*i+j])*(y0[k] - x0[D*i+k])*p/z - (y0[j] - x0[D*i+j])*p/z*dlnZ_dL[k];
			}
		}
	}

	free(dlnZ_dL);
}

void MaxEnt_dlnZ_dY(double *dlnZ_dY, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double p = 0;
	double z = MaxEnt_Z(lambda, y0, beta0, D, N, x0);
	double *A = (double *) malloc(D2*sizeof(double)); assert(A!=NULL);
	double *AA = (double *) malloc(D2*sizeof(double)); assert(AA!=NULL);
	double *bh = (double *) malloc(D*sizeof(double)); assert(bh!=NULL);
	for(int i = 0; i < D; ++i){
		dlnZ_dY[i] = 0;
		for(int j = 0; j < D; ++j){
			AA[i*D+j] = A[D*i+j] + A[D*j+i];	
		}
	}

	for(int i = 0; i < N; ++i){
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, D, D, D, 1, beta0+i*D2, D, beta0+i*D2, D, 0, A, D);
		for(int j = 0; j < D; ++j){
			for(int k = 0; k < D; ++k){
				AA[j*D+k] = A[D*j+k] + A[D*k+j];	
			}
		}

		double ds2 = 0;
		double ds1 = 0;
		for(short j = 0; j < D; ++j){
			double hj = 0;
			double aj = 0;
			for(short k = 0; k < D; ++k){
				hj += beta0[i*D2+j*D+k]*(y0[k] - x0[D*i+k]);
				aj += AA[j*D+k]*(y0[k] - x0[D*i+k]);
			}
			//hj = fabs(hj-nearbyint(hj));
			ds1 += lambda[j]*(y0[j] - x0[D*i+j]);
			ds2 += hj*hj;
			bh[j] = aj;
		}
		p = exp(-ds2 + ds1);
		for(short j = 0; j < D; ++j){
			dlnZ_dY[j] += (-bh[j] + lambda[j])*p/z;
		}
	}
	free(A);
	free(AA);
	free(bh);
}

void MaxEnt_hlnZ_hY(double *hlnZ_hY, const double *lambda, const double *y0, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double p = 0;
	double *dlnZ_dY = (double *) malloc(D*sizeof(double)); assert(dlnZ_dY!=NULL);

	double z = MaxEnt_Z(lambda, y0, beta0, D, N, x0);
	MaxEnt_dlnZ_dY(dlnZ_dY, lambda, y0, beta0, D, N, x0);

	double *A = (double *) malloc(D2*sizeof(double)); assert(A!=NULL);
	double *AA = (double *) malloc(D2*sizeof(double)); assert(AA!=NULL);
	double *bh = (double *) malloc(D*sizeof(double)); assert(bh!=NULL);

	for(short i = 0; i < D; ++i){
		for(short j = 0; j < D; ++j){
			hlnZ_hY[D*i+j] = 0;
		}
	}

	for(int i = 0; i < N; ++i){
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, D, D, D, 1, beta0+i*D2, D, beta0+i*D2, D, 0, A, D);
		for(short j = 0; j < D; ++j){
			for(short k = 0; k < D; ++k){
				AA[D*j+k] = A[D*j+k] + A[D*k+j];	
			}
		}

		double ds2 = 0;
		double ds1 = 0;
		for(short j = 0; j < D; ++j){
			double hj = 0;
			double aj = 0;
			for(short k = 0; k < D; ++k){
				hj += beta0[i*D2+j*D+k]*(y0[k] - x0[D*i+k]);
				aj += AA[j*D+k]*(y0[k] - x0[D*i+k]);
			}
			//hj = fabs(hj-nearbyint(hj));
			ds1 += lambda[j]*(y0[j] - x0[D*i+j]);
			ds2 += hj*hj;
			bh[j] = aj;
		}	
		p = exp(-ds2 + ds1);
		for(short j = 0; j < D; ++j){
			for(short k = 0; k < D; ++k){
				hlnZ_hY[D*j+k] += (y0[j] - x0[D*i+j])*(-bh[k] + lambda[k])*p/z - (y0[j] - x0[D*i+j])*p/z*dlnZ_dY[k];
			}
		}
	}

	for(short j = 0; j < D; ++j){
		for(short k = 0; k < D; ++k){
			double delta = (j==k) ? 1:0;
			hlnZ_hY[D*j+k] += delta;
		}
	}
	free(dlnZ_dY);
	free(A);
	free(AA);
	free(bh);
}
/*
void mat_inv2x2(double *mat){
	double a = mat[0];
	double b = mat[1];
	double c = mat[2];
	double d = mat[3];
	
	double det = a*d - b*c;
	mat[0] = d/det;
	mat[1] = -b/det;
	mat[2] = -c/det;
	mat[3] = a/det;
}

void mat_inv3x3(double *A){
	double a11 = A[0];
	double a12 = A[1];
	double a13 = A[2];

	double a21 = A[3];
	double a22 = A[4];
	double a23 = A[5];

	double a31 = A[6];
	double a32 = A[7];
	double a33 = A[8];

	double D1 = a11*a22*a33;
	double D2 = a21*a32*a13;
	double D3 = a31*a12*a23;
	double D4 = -a31*a22*a13;
	double D5 = -a12*a21*a33;
	double D6 = -a11*a23*a32;

	double D = D1 + D2 + D3 + D4 + D5 + D6;

	A[0] = (a22*a33 - a32*a23)/D;
	A[1] = (a13*a32 - a33*a12)/D;
	A[2] = (a12*a23 - a22*a13)/D;
	
	A[3] = (a23*a31 - a33*a21)/D;
	A[4] = (a11*a33 - a13*a31)/D;
	A[5] = (a13*a21 - a11*a23)/D;

	A[6] = (a21*a32 - a31*a22)/D;
	A[7] = (a12*a31 - a32*a11)/D;
	A[8] = (a11*a22 - a21*a12)/D;
}*/

void mat_inv(double *A, const int dim){
	assert(dim>=1&&dim<=3);
	switch (dim){
		case 1:{
			double a = A[0];
			//assert(a!=0);
			A[0] = 1.0/a;
			break;
		}
		case 2:{
			double a = A[0];
			double b = A[1]; 
			double c = A[2];
			double d = A[3];
			
			double det = a*d - b*c;
			//assert(det!=0);
			A[0] = d/det;
			A[1] = -b/det;
			A[2] = -c/det;
			A[3] = a/det;
			break;
		}
		case 3:{
			double a11 = A[0];
			double a12 = A[1];
			double a13 = A[2];

			double a21 = A[3];
			double a22 = A[4];
			double a23 = A[5];

			double a31 = A[6];
			double a32 = A[7];
			double a33 = A[8];

			double D1 = a11*a22*a33;
			double D2 = a21*a32*a13;
			double D3 = a31*a12*a23;
			double D4 = -a31*a22*a13;
			double D5 = -a12*a21*a33;
			double D6 = -a11*a23*a32;

			double D = D1 + D2 + D3 + D4 + D5 + D6;
			//assert(D!=0);
			A[0] = (a22*a33 - a32*a23)/D;
			A[1] = (a13*a32 - a33*a12)/D;
			A[2] = (a12*a23 - a22*a13)/D;
			
			A[3] = (a23*a31 - a33*a21)/D;
			A[4] = (a11*a33 - a13*a31)/D;
			A[5] = (a13*a21 - a11*a23)/D;

			A[6] = (a21*a32 - a31*a22)/D;
			A[7] = (a12*a31 - a32*a11)/D;
			A[8] = (a11*a22 - a21*a12)/D;
			break;
		}
	}
}

void MaxEnt_dP_dY(double *dP_dY, const double *y0, const double *lambda, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double p = 0;
	double z = MaxEnt_Z(lambda, y0, beta0, D, N, x0);
	double *hlnZ_hL = (double *) malloc(D2*sizeof(double)); assert(hlnZ_hL!=NULL);
	double *hlnZ_hY = (double *) malloc(D2*sizeof(double)); assert(hlnZ_hY!=NULL);
	double *A = (double *) malloc(D2*sizeof(double)); assert(A!=NULL);
	double *AA = (double *) malloc(D2*sizeof(double)); assert(AA!=NULL);

	for(int i = 0; i < N; ++i){
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, D, D, D, 1, beta0+i*D2, D, beta0+i*D2, D, 0, A, D);

		for(int j = 0; j < D; ++j){
			for(int k = 0; k < D; ++k){
				AA[j*D+k] = A[D*j+k] + A[D*k+j];	
			}
		}

		MaxEnt_hlnZ_hL(hlnZ_hL, lambda, y0, beta0, D, N, x0);
		mat_inv(hlnZ_hL, D);
		/*lapack_int *ipiv = (lapack_int*) malloc(D*sizeof(lapack_int)); assert(ipiv!=NULL);
		lapack_int info;
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, D, D, hlnZ_hL, D, ipiv);
		info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, D, hlnZ_hL, D, ipiv);*/

		MaxEnt_hlnZ_hY(hlnZ_hY, lambda, y0, beta0, D, N, x0);

		for(int j = 0; j < D; ++j){
			for(int k = 0; k < D; ++k){
				A[D*j+k] = AA[D*j+k];	
			}
		}

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, D, D, D, 1, hlnZ_hL, D, hlnZ_hY, D, 1, A, D);


		double ds2 = 0;
		double ds1 = 0;
		for(short j = 0; j < D; ++j){
			double hj = 0;
			for(short k = 0; k < D; ++k){
				hj += beta0[i*D2+j*D+k]*(y0[k] - x0[D*i+k]);
			}
			//hj = fabs(hj-nearbyint(hj));
			ds1 += lambda[j]*(y0[j] - x0[D*i+j]);
			ds2 += hj*hj;
		}
		p = exp(-ds2 + ds1);
		for(short j = 0; j < D; ++j){
			double dp = 0;
			for(short k = 0; k < D; ++k){
				dp += A[D*j+k]*(y0[k] - x0[D*i+k]);
			}
			dP_dY[D*i+j] = -p*dp/z;
		}
	}

	free(hlnZ_hL);
	free(hlnZ_hY);
	free(A);
	free(AA);
	//free(ipiv);
	//free(work);
}

void MaxEnt_F(double *F, const double *x1, const double *y0, const double *lambda, const double *beta0, const int D, const int N, const double *x0){
	double *dP_dY = (double *) malloc(N*D*sizeof(double)); assert(dP_dY!=NULL);
	MaxEnt_dP_dY(dP_dY, y0, lambda, beta0, D, N, x0);
	for(short j = 0; j < D; ++j){
		for(short k = 0; k < D; ++k){
			F[j*D+k] = 0;
		}
	}
	for(int i = 0; i < N; ++i){
		for(short j = 0; j < D; ++j){
			for(short k = 0; k < D; ++k){
				F[j*D+k] += dP_dY[D*i+k]*x1[D*i+j];
			}
		}	
	}
	free(dP_dY);
}

void MaxEnt_NR(double *lambda0, const double *y0, const int maxsteps, const double eps, const double *beta0, const int D, const int N, const double *x0){
	const int D2 = D*D;
	double *dL= (double *) malloc(D*sizeof(double)); assert(dL!=NULL);
	double *hL= (double *) malloc(D2*sizeof(double)); assert(hL!=NULL);
	for(int n = 0; n < maxsteps; ++n){
		MaxEnt_dlnZ_dL(dL, lambda0, y0, beta0, D, N, x0);
		MaxEnt_hlnZ_hL(hL, lambda0, y0, beta0, D, N, x0);

		double d2p = 0;
		for(short i = 0; i < D; ++i){
			d2p += dL[i]*dL[i];
		}
		if(d2p < eps/D){
			break;
		}

		mat_inv(hL, D);

		for(short i = 0; i < D; ++i){
			double dpi = 0;
			for(short j = 0; j < D; ++j){
				dpi += hL[D*i+j]*dL[j];
			}
			lambda0[i] -= dpi;
		}	
	}
	free(dL);
	free(hL);
}
