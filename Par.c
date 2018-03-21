#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "vtkCellType.h"

typedef struct {
	int num;
	double point[3];
	double director[9];
	int lattice[3];
} Cell;


int main(){
	double L1 = 4;
	double L2 = 4;
	double L3 = 4;

	double origin[3] = {0, 0, 0};
	double a1[3] = {1 ,0, 0};
	double a2[3] = {0, 1, 0};
	double a3[3] = {0, 0, 1};

	double a1_n = 0;
	double a2_n = 0;
	double a3_n = 0;

	for(int i = 0; i < 3; ++i){
		a1_n += a1[i]*a1[i];
		a2_n += a2[i]*a2[i];
		a3_n += a3[i]*a3[i];
	}
	a1_n = sqrt(a1_n);
	a2_n = sqrt(a2_n);
	a3_n = sqrt(a3_n);


	double *director = (double *) malloc(9*sizeof(double));
	for(short j = 0; j < 3; ++j){
		director[3*j] = a1[j]/2;
		director[3*j+1] = a2[j]/2;
		director[3*j+2] = a3[j]/2;
	}

	int N1 = (int) round(L1/a1_n);
	int N2 = (int) round(L2/a2_n);
	int N3 = (int) round(L3/a3_n);
	int N = N1*N2*N3;

	Cell *Cells = (Cell *) malloc(N*sizeof(Cell));

	for(int i = 0; i < N1; ++i){
		for(int j = 0; j < N2; ++j){
			for(int k = 0; k < N3; ++k){
				int n = i*N2*N3 + j*N3 + k;
				Cells[n].num = n;
				for(short a = 0; a < 3; ++a){
					Cells[n].point[a] = origin[a] + i*a1[a] + j*a2[a] + k*a3[a];
				}
				for(short b = 0; b < 3; ++b){
					for(short c = 0; c < 3; ++c){
						Cells[n].director[3*b+c] = director[3*b+c];
					}
				}
				for(short d = 0; d < 3; ++d){
					Cells[n].lattice[d] = 2;
				}
			}
		}	
	}
	
	int np = 0;
	int nc = 0;
	int ncsize = 0;
	int ntypes = 13;
	short celltypes[ntypes];
	for(short i = 0; i < ntypes; ++i){
		celltypes[i] = 0;
	}
	for(int n = 0; n < N; ++n){
		int l[3] = {0, 0, 0};
		for(short d = 0; d < 3; ++d){
			l[d] = Cells[n].lattice[d];
		}
		int vtk_celltype = 0;
		int dim = 0;
		double ds2[3] = {0, 0, 0};
		for(short i = 0; i < 3; ++i){
			double di_n = 0;
			double li = l[i];
			for(short j = 0; j < 3; ++j){
				double dij = Cells[n].director[3*j+i];
				di_n += li*li*dij*dij;
			}
			ds2[i] = di_n;
			if(di_n > 0){
				dim++;
			}
		}
		int vertices = 0;
		int subcells = 0;
		switch (dim){
			case 0:
				//vtk_celltype = VTK_VERTEX;
				vtk_celltype = 1;
				vertices = 1;
				subcells = 1;
				break;
			case 1:
				//vtk_celltype = VTK_LINE;
				vtk_celltype = 3;
				vertices = 2;
				for(short d = 0; d < 3; ++d){
					if(ds2[d]!=0){
						subcells = l[d]; 
					}
				}
				break;
			case 2:
				//vtk_celltype = VTK_QUAD;
				vtk_celltype = 9;
				vertices = 4;
				for(short d = 0; d < 3; ++d){
					if(ds2[d]==0){
						subcells = l[(d+1)%3]*l[(d+2)%3]; 
					}
				}

				break;
			case 3:
				//vtk_celltype = VTK_HEXAHEDRON;
				vtk_celltype = 12;
				vertices = 8;
				subcells = l[0]*l[1]*l[2];
				break;
		}
		celltypes[vtk_celltype] += subcells;

		np += (1+l[0])*(1+l[1])*(1+l[2]);
		nc += subcells;
		ncsize += subcells*vertices;
	}
	for(int i = 1; i < ntypes; ++i){
		celltypes[0] += celltypes[i];
	}
	char header[256] = "# vtk DataFile Version 2.0";
	char title[256] = "Test";
	char datatype[128] = "ASCII";
	char dataset[128] = "UNSTRUCTURED_GRID";
	int point_data = np;
	int cell_data = nc;
	int cell_size = ncsize;
	FILE *pf = fopen("Cells.vtk", "w");
	fprintf(pf, "%s\n%s\n%s\n", header, title, datatype);
	fprintf(pf, "DATASET %s\n", dataset);

	fprintf(pf, "POINTS %d %s\n", point_data, "float");
	for(int n = 0; n < N; ++n){
		double p[3] = {0, 0, 0};
		for(short k = 0; k < 3; ++k){
			p[k] = Cells[n].point[k];
			//fprintf(pf, "%f ", p[k]);
		}
		//fprintf(pf, "\t");
		int L[3] = {0, 0, 0};
		for(short l = 0; l < 3; ++l){
			L[l] = Cells[n].lattice[l];
		}
		double A[9] = {0, 0, 0, 
					   0, 0, 0, 
					   0, 0, 0};
		for(short a = 0; a < 3; ++a){
			for(short j = 0; j < 3; ++j){
				A[3*a+j] = Cells[n].lattice[a]*Cells[n].director[3*a+j];
			}
		}
		for(int l1 = 0; l1 <= L[0]; ++l1){
			for(int l2 = 0; l2 <= L[1]; ++l2){
				for(int l3 = 0; l3 <= L[2]; ++l3){
					for(short i = 0; i < 3; ++i){
						double pl = p[i] + l1*A[3*i] + l2*A[3*i+1] + l3*A[3*i+2];
						fprintf(pf, "%e ", pl);
					}
					fprintf(pf, "\n");
				}
			}
		}
	}
	fprintf(pf, "\n");

	fprintf(pf, "CELLS %d %d\n", cell_data, cell_data + cell_size);
	np = 0;
	nc = 0;
	for(int n = 0; n < N; ++n){
		int L[3] = {0, 0, 0};
		for(short d = 0; d < 3; ++d){
			L[d] = Cells[n].lattice[d];
		}
		int dim = 0;
		double ds2[3] = {0, 0, 0};
		for(short i = 0; i < 3; ++i){
			double di_n = 0;	
			double li = L[i];
			for(short j = 0; j < 3; ++j){
				double dij = Cells[n].director[3*j+i];
				di_n += li*li*dij*dij;
			}
			ds2[i] = di_n;
			if(di_n > 0){
				dim++;
			}
		}
		int vertices = 0;
		switch (dim){
			case 0:
				vertices = 1;
				fprintf(pf, "%d ", vertices);
				fprintf(pf, "%d \n", np);
				break;
			case 1:
				vertices = 2;
				for(int l = 0; l < 3; ++l){
					if(ds2[l]!=0){
						for(int l1 = 0; l1 < L[l]; ++l1){	
							int pi = np + l1;	
							fprintf(pf, "%d ", vertices);
							fprintf(pf, "%d %d\n", pi, pi+1);
						}
						break;
					}
				}
				break;
			case 2:
				vertices = 4;
				for(int l = 0; l < 3; ++l){
					if(ds2[l]==0){
						int L1 = 0;
						int L2 = 0;
						switch (l){
							case 0:
								L1 = 1;
								L2 = 2;
								break;
							case 1:
								L1 = 0;
								L2 = 2;
								break;
							case 2:
								L1 = 0;
								L2 = 1;
								break;
						}
						int pi[4] = {0, 0, 0, 0};
						for(int l1 = 0; l1 < L[L1]; ++l1){
							for(int l2 = 0; l2 < L[L2]; ++l2){
								pi[0] = np + l1*(L[L2]+1) + l2;	
								pi[1] = np + l1*(L[L2]+1) + l2 + 1;	
								pi[2] = np + (l1+1)*(L[L2]+1) + l2;	
								pi[3] = np + (l1+1)*(L[L1]+1) + l2 + 1;	
								fprintf(pf, "%d ", vertices);
								for(short j = 0; j < 4; ++j){	
									fprintf(pf, "%d ", pi[j]);
								}
								fprintf(pf, "\n");
							}
						}
					}
				}
				break;
			case 3:
				vertices = 8;
				int pi[8] = {0, 0, 0, 0, 0, 0, 0, 0};
				for(int l1 = 0; l1 < L[0]; ++l1){
					for(int l2 = 0; l2 < L[1]; ++l2){
						for(int l3 = 0; l3 < L[2]; ++l3){
							pi[0] = np + l1*(L[1]+1)*(L[2]+1) + l2*(L[2]+1) + l3;	
							pi[1] = np + l1*(L[1]+1)*(L[2]+1) + l2*(L[2]+1) + l3 + 1;	
							pi[2] = np + l1*(L[1]+1)*(L[2]+1) + (l2+1)*(L[2]+1) + l3 + 1;	
							pi[3] = np + l1*(L[1]+1)*(L[2]+1) + (l2+1)*(L[2]+1) + l3;	
							pi[4] = np + (l1+1)*(L[1]+1)*(L[2]+1) + l2*(L[2]+1) + l3;	
							pi[5] = np + (l1+1)*(L[1]+1)*(L[2]+1) + l2*(L[2]+1) + l3 + 1;	
							pi[6] = np + (l1+1)*(L[1]+1)*(L[2]+1) + (l2+1)*(L[2]+1) + l3 + 1;	
							pi[7] = np + (l1+1)*(L[1]+1)*(L[2]+1) + (l2+1)*(L[2]+1) + l3;	

							fprintf(pf, "%d ", vertices);
							for(short j = 0; j < 8; ++j){	
								fprintf(pf, "%d ", pi[j]);
							}
							fprintf(pf, "\n");
						}
					}	
				}
				break;
		}
		np += (1+L[0])*(1+L[1])*(1+L[2]);
	}

	fprintf(pf, "CELL_TYPES %d\n", celltypes[0]);
	for(short t = 1; t < ntypes; ++t){
		for(int i = 0; i < celltypes[t]; ++i){
			if(celltypes[t]!=0){
				fprintf(pf, "%d\n", t);
			}
		}
	
	}

	//fprintf(pf, "POINT_DATA %d", point_data);
	//fprintf(pf, "")
	
	return 0;
}
