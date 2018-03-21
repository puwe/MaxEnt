#include <time.h>
#include "MaxEnt.h"

// Y <- A*X + b
void Affine(double *Y, double *A, int m, double *X, int n, const double *b){
	//double *Af = (double *) malloc(m*n*sizeof(double)); assert(Af!=NULL);
	double alpha = 1.0;
	double beta = 1.0;
	for(int i = 0; i < m; ++i){
		Y[i] = b[i];
		//for(int j = 0; j < n; ++j){
		//	Af[i+j*m] = A[n*i+j];
		//}
	}
	//char trans = 'n';
	int lda = n;
	int incx = 1;
	int incy = 1;
	//dgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, A, lda, X, incx, beta, Y, incy);
	//free(Af);
}


void test_affinetrans(){
	const int D = 2;

	double L1 = 1;
	double L2 = 1;

	int Na = 9;
	int Nb = 100;

	double ha = (L2+L1)*2/(Na-1);
	double hb = (L2+L1)/(Nb-1);

	int nax = (int) round(L1/ha+0.5);
	int nbx = (int) round(L1/4/hb+0.5);
	
	double *Xa = (double *) malloc(Na*D*sizeof(double));
	printf("Xa\n");
	for(int n = 0; n < Na; ++n){
		Xa[D*n] = n%nax*ha;
		Xa[D*n+1] = n/nax*ha;
	}
	for(int n = 0; n < Na; ++n){
		for(int i = 0; i < D; ++i){
			printf("%f ", Xa[D*n+i]);
		}
		printf("\n");
	}
	printf("Xb\n");
	double *Xb = (double *) malloc(Nb*D*sizeof(double));
	for(int n = 0; n < Nb; ++n){
		Xb[D*n] = L1/2 + n%nbx*hb;
		Xb[D*n+1] = L2/2 + n/nbx*hb;
	}
	for(int n = 0; n < Nb; ++n){
		for(int i = 0; i < D; ++i){
			printf("%f ", Xb[D*n+i]);
		}
		printf("\n");
	}
	// beta^2 = 1/(2*sig^2) sig = 2pi/a
	double a0 = 1;
	double *lambda = (double*) malloc((D)*sizeof(double));
	double *Lambda = (double*) malloc(D*Nb*sizeof(double));
	double *beta = (double*) malloc((D*D)*sizeof(double));
	double *Beta = (double*) malloc((D*D)*Nb*sizeof(double));
	for(int i = 0; i < D; ++i){	
		lambda[i] = 0; //(double) rand()/RAND_MAX;
		for(int j = 0; j < D; ++j){
			beta[D*i+j] = a0/(sqrt(2)*2*M_PI);  
		}
	}


	FILE *pf = fopen("MaxEntAtoms.txt", "w"); assert(pf!=NULL);
	for(int n = 0; n < Na; ++n){
		for(int i = 0; i < D; ++i){
			fprintf(pf, "%f ", Xa[D*n+i]);
		}
		fprintf(pf, "\n");
	}
	fclose(pf);

	pf = fopen("MaxEntNodes.txt", "w"); assert(pf!=NULL);
	for(int n = 0; n < Nb; ++n){
		for(int i = 0; i < D; ++i){
			fprintf(pf, "%f ", Xb[D*n+i]);
		}
		fprintf(pf, "\n");
	}
	fclose(pf);

	pf = fopen("MaxEntParams.txt", "w"); assert(pf!=NULL);
	for(int n = 0; n < Nb; ++n){		
	
		MaxEnt_NR(lambda, Xb+D*n, 2e3, 1e-6, beta, D, Na, Xa);
				
		//printf("\n");
		//for(short i = 0; i < D; ++i){
		//	printf("%f \n", lambda[i]);
		//}
		for(int i = 0; i < D; ++i){
			Lambda[D*n+i] = lambda[i];
			fprintf(pf, "%f ", lambda[i]);
		}
		for(int i = 0; i < D*D; ++i){
			Beta[D*D*n+i] = beta[i];
			fprintf(pf, "%f ", beta[i]);	
		}
		fprintf(pf, "\n");
	}
	fclose(pf);

	lapack_int info;
	double *A = (double*) malloc(D*D*sizeof(double));
	double *b = (double*) malloc(D*sizeof(double));
	//double *ba = (double*) malloc(D*Na*sizeof(double));
	double *Xa1 = (double *) malloc(Na*D*sizeof(double));
	double *Xb1 = (double *) malloc(Nb*D*sizeof(double));
	double *Beta1 = (double*) malloc((D*D)*Nb*sizeof(double));
	double *Lambda1 = (double*) malloc(D*Nb*sizeof(double));
	for(int n = 0; n < Nb; ++n){
		for(int i = 0; i < D; ++i){
			Lambda1[n*D+i] = 0;
			for(int j = 0; j < D; ++j){
				Beta1[n*D*D+i*D+j] = 0;
			}
		}
	}
	srand(time(NULL));
	for(int i = 0; i < D; ++i){
		double r = (double) rand()/RAND_MAX;
		b[i] = r;
	}
	lapack_int iseed[4] = {1, 3, 5, 7};

	//info = LAPACKE_dlagge(LAPACK_ROW_MAJOR, D, D, D-1, D-1, b, A, D, iseed);
	info = LAPACKE_dlagsy(LAPACK_ROW_MAJOR, D, D-1, b, A, D, iseed);
	for(int i = 0; i < D; ++i){
		for(int j = 0; j < D; ++j){
			//A[i*D+j] = (i==j) ? 1:0;
			printf("%6.2f", A[i*D+j]);
		}
		printf("\n");
	}
	//A[1] = 0.5;

	pf = fopen("MaxEntAtoms1.txt", "w"); assert(pf!=NULL);
	for(int n = 0; n < Na; ++n){
		for(int i = 0; i < D; ++i){
			double r = (double) rand()/RAND_MAX;
			b[i] = 0;//-r*(L1+L2);
		}

		Affine(Xa1+D*n, A, D, Xa+D*n, D, b);
		for(int i = 0; i < D; ++i){
			fprintf(pf, "%f ", Xa1[D*n+i]);
		}
		fprintf(pf, "\n");
	}
	fclose(pf);

	pf = fopen("MaxEntParams1.txt", "w"); assert(pf!=NULL);
	//lapack_int info;
	lapack_int* ipiv = (lapack_int*) malloc(D*sizeof(lapack_int));
	//double *beta1 = (double*) malloc(D*D*sizeof(double));
	double *lambda1 = (double*) malloc(D*sizeof(double));
	double *Pa = (double*) malloc(Na*sizeof(double));
	//double X[3];
	double *F = (double*) malloc(D*D*sizeof(double));
	for(int n = 0; n < Nb; ++n){		
	
		MaxEnt_Y(Xb1+n*D, Pa, Xa1, Xb+D*n, Lambda+n*D, Beta+n*D*D, D, Na, Xa);
		MaxEnt_F(F, Xa1, Xb+D*n, Lambda+n*D, Beta+n*D*D, D, Na, Xa);
		//printf("%d ", n);
		//for(int i = 0; i < D; ++i){
		//	double xbi = 0;
		//	for(int na = 0; na < Na; ++na){	
		//		xbi += Pa[na]*Xa1[na*D+i];
		//	}
		//	Xb1[n*D+i] = xbi;
		//	printf("%f ", xbi);
		//}
		//for(int i = 0; i < D; ++i){
		//	printf("%f ", X[i]);
		//}
		for(int i = 0; i < D*D; ++i){
			printf("%f ", F[i]);
		}
		printf("\n");

		MaxEnt_NR(lambda1, Xb1+D*n, 200, 1e-6, Beta+n*D*D, D, Na, Xa1);
		//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, D, D, D, 1.0, A, D, Beta+n*D*D, D, 1.0, Beta1+n*D*D, D);	
		for(int i = 0; i < D*D; ++i){
			Beta1[n*D*D+i] = Beta[n*D*D+i];
		}
		//info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, D, D, A, D, ipiv, Beta1+n*D*D, D);	
		//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, D, D, D, 1.0, A, D, Lambda+n*D*D, D, 1.0, Lambda1+n*D*D, D);	
		for(int i = 0; i < D; ++i){
			Lambda1[n*D+i] = lambda1[i]; //Lambda[n*D+i];
		}
		//info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, D, 1, A, D, ipiv, Lambda1+n*D, 1);	

		//printf("\n");
		//for(short i = 0; i < D; ++i){
		//	printf("%f \n", Lambda[D*n+i]);
		//}
		for(int i = 0; i < D; ++i){
			fprintf(pf, "%f ", Lambda1[D*n+i]);
		}
		for(int i = 0; i < D*D; ++i){
			fprintf(pf, "%f ", Beta1[D*D*n+i]);	
		}
		fprintf(pf, "\n");

	}
	fclose(pf);

	pf = fopen("MaxEntNodes1.txt", "w"); assert(pf!=NULL);
	for(int n = 0; n < Nb; ++n){
		//for(int i = 0; i < D; ++i){
			//double r = (double) rand()/RAND_MAX;
			//b[i] = 0;
		//}

		//Affine(Xb1+D*n, A, D, Xb+D*n, D, b);
		for(int i = 0; i < D; ++i){
			fprintf(pf, "%f ", Xb1[D*n+i]);
		}
		fprintf(pf, "\n");
	}
	fclose(pf);


}

void test_lattice(){
	const int D = 2;
	const int D2 = D*D;
	double L1 = 3;
	double L2 = 3;

	int N1 = 3;
	int N2 = 3;
	int N = N1*N2;

	double H1 = L1/(N1-1);
	double H2 = L2/(N2-1);

	int n1 = 10;
	int n2 = 10;
	int n = n1*n2;
	
	double eps = 1e-3;
	double h1 = (L1-2*eps)/(n1-1);
	double h2 = (L2-2*eps)/(n2-1);

	double *x0 = (double *) malloc(N1*N2*D*sizeof(double));
	for(int i = 0; i < N1; ++i){
		double x = i*H1;
		for(int j = 0; j < N2; ++j){
			double y = j*H2;
			x0[D*(i*N2+j)] = x;
			x0[D*(i*N2+j)+1] = y;
		
		}
	}

	double *y0 = (double *) malloc(n1*n2*D*sizeof(double));
	for(int i = 0; i < n1; ++i){
		double x = eps + i*h1;
		for(int j = 0; j < n2; ++j){
			double y = eps + j*h2;
			y0[D*(i*n2+j)] = x;
			y0[D*(i*n2+j)+1] = y;
		
		}
	}

	// beta^2 = 1/(2*sig^2) sig = 2pi/a
	double a0 = 1;
	double lambda[D];
	double *Lambda = (double*) malloc(D*n*sizeof(double));
	double *Beta = (double*) malloc(D2*N*sizeof(double));

	double *beta = (double*) malloc(D2*n*sizeof(double));
	for(short j = 0; j < D; ++j){
			lambda[j] = 0; //(double) rand()/RAND_MAX;
	}

	for(int i = 0; i < n; ++i){	
		for(short j = 0; j < D; ++j){
			Lambda[i*D+j] = 0; //(double) rand()/RAND_MAX;
		}
	}

	for(int i = 0; i < N; ++i){	
		for(int j = 0; j < D; ++j){
			for(short k = 0; k < D; ++k){
				Beta[D2*i+j*D+k] = (j==k) ? a0/(sqrt(2)*2*M_PI):0;  
			}
		}
	}
	
	for(int i = 0; i < n; ++i){	
		for(int j = 0; j < D; ++j){
			for(short k = 0; k < D; ++k){
				beta[D2*i+j*D+k] = (j==k) ? a0/(sqrt(2)*2*M_PI):0;  
			}
		}
	}

	FILE *pf = fopen("MaxEntAtoms.txt", "w"); assert(pf!=NULL);
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < D; ++j){
			fprintf(pf, "%f ", x0[D*i+j]);
		}
		for(int k = 0; k < D2; ++k){
			fprintf(pf, "%f ", Beta[i*D2+k]);
		}
		fprintf(pf, "\n");
	}
	fclose(pf);

	pf = fopen("MaxEntNodes.txt", "w"); assert(pf!=NULL);
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < D; ++j){
			fprintf(pf, "%f ", y0[D*i+j]);
		}
		for(int k = 0; k < D2; ++k){
			fprintf(pf, "%f ", beta[i*D2+k]);
		}

		fprintf(pf, "\n");
	}
	fclose(pf);

	pf = fopen("MaxEntParams.txt", "w"); assert(pf!=NULL);
	for(int i = 0; i < n; ++i){		
		for(int j = 0; j < D; ++j){
			lambda[j] = 0; 
		}
		MaxEnt_NR(lambda, y0+D*i, 2e3, 1e-6, Beta, D, N, x0);
		
		for(int j = 0; j < D; ++j){
			Lambda[i*D+j] = lambda[j];
			fprintf(pf, "%f ", lambda[j]);
			//printf("%f ", lambda[j]);
		}
	
		fprintf(pf, "\n");
		//printf("\n");
	}
	fclose(pf);

	double *x1 = (double *) malloc(N1*N2*D*sizeof(double));
	double *y1 = (double *) malloc(n1*n2*D*sizeof(double));
	pf = fopen("MaxEntAtoms1.txt", "w"); assert(pf!=NULL);
	double b[D];
	double A[D2];
	for(int j = 0; j < D; ++j){
		double r = (double) rand()/RAND_MAX;
		b[j] = r;//-r*(L1+L2);
	}
	lapack_int info;
	lapack_int iseed[4] = {1, 2, 3, 4};
	info = LAPACKE_dlagge(LAPACK_ROW_MAJOR, D, D, D-1, D-1, b, A, D, iseed);
	//info = LAPACKE_dlagsy(LAPACK_ROW_MAJOR, D, D-1, b, A, D, iseed);
	for(int i = 0; i < D; ++i){
		for(int j = 0; j < D; ++j){
			//A[i*D+j] = (i==j) ? 1:0;
			printf("%6.2f", A[i*D+j]);
		}
		printf("\n");
	}

	double *Beta1 = (double*) malloc(D*D*N*sizeof(double));
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < D; ++j){
			double r = (double) rand()/RAND_MAX;
			x1[D*i+j] = r;
		}

		cblas_dgemv(CblasRowMajor, CblasNoTrans, D, D, 1.0, A, D, x0+i*D, 1, 1.0, x1+i*D, 1 );
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, D, D, D, 1.0, A, D, Beta+i*D2, D, 0.0, Beta1+i*D2, D);
		
		for(int j = 0; j < D; ++j){
			fprintf(pf, "%f ", x1[D*i+j]);
		}
		for(int k = 0; k < D2; ++k){
			fprintf(pf, "%f ", Beta1[D2*i+k]);
		}

		fprintf(pf, "\n");
	}
	fclose(pf);
	
	pf = fopen("MaxEntParams1.txt", "w"); assert(pf!=NULL);
	double *lambda1 = (double*) malloc(D*sizeof(double));
	double *Lambda1 = (double*) malloc(D*n*sizeof(double));
	double *beta1 = (double*) malloc(D*D*n*sizeof(double));

	double *Pa = (double*) malloc(N*sizeof(double));
	double *F = (double*) malloc(D*D*sizeof(double));
	for(int i = 0; i < n; ++i){		
	
		MaxEnt_Y(y1+i*D, Pa, x1, y0+D*i, Lambda+i*D, Beta, D, N, x0);
		MaxEnt_F(F, x1, y0+D*i, Lambda+i*D, Beta, D, N, x0);
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, D, D, D, 1.0, F, D, beta+i*D2, D, 0.0, beta1+i*D2, D);
		for(int j = 0; j < D2; ++j){
			printf("%f ", F[j]);
		}
		printf("\n");
		for(int j = 0; j < D; ++j){
			lambda1[j] = 0; 
		}
		MaxEnt_NR(lambda1, y1+D*i, 200, 1e-6, Beta1, D, N, x1);
		
		for(int j = 0; j < D; ++j){
			Lambda1[i*D+j] = lambda1[j]; //Lambda[n*D+i];
		}

		for(int j = 0; j < D; ++j){
			fprintf(pf, "%f ", Lambda1[D*i+j]);
		}
		fprintf(pf, "\n");

	}
	fclose(pf);

	pf = fopen("MaxEntNodes1.txt", "w"); assert(pf!=NULL);
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < D; ++j){
			fprintf(pf, "%f ", y1[D*i+j]);
		}
		for(int j = 0; j < D2; ++j){
			fprintf(pf, "%f ", beta1[D2*i+j]);	
		}
		fprintf(pf, "\n");
	}
	fclose(pf);



}

int main(){
	srand(time(NULL));
	test_lattice();

}
