/*
	Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright (C) 2011  Xiang Zhou

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <cmath>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"

using namespace std;

extern "C" void sgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA, float *A, int *LDA, float *B, int *LDB, float *BETA, float *C, int *LDC);
extern "C" void spotrf_(char *UPLO, int *N, float *A, int *LDA, int *INFO);
extern "C" void spotrs_(char *UPLO, int *N, int *NRHS, float *A, int *LDA, float *B, int *LDB, int *INFO);
extern "C" void ssyev_(char* JOBZ, char* UPLO, int *N, float *A, int *LDA, float *W, float *WORK, int *LWORK, int *INFO);
extern "C" void ssyevr_(char* JOBZ, char *RANGE, char* UPLO, int *N, float *A, int *LDA, float *VL, float *VU, int *IL, int *IU, float *ABSTOL, int *M, float *W, float *Z, int *LDZ, int *ISUPPZ, float *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);
extern "C" void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);
extern "C" void dpotrs_(char *UPLO, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, int *INFO);
extern "C" void dsyev_(char* JOBZ, char* UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO);
extern "C" void dsyevr_(char* JOBZ, char *RANGE, char* UPLO, int *N, double *A, int *LDA, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);


//cholesky decomposition, A is distroyed
void lapack_float_cholesky_decomp (gsl_matrix_float *A)
{
	int N=A->size1, LDA=A->size1, INFO;
	char UPLO='L';
	
	if (N!=(int)A->size2) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_decomp."<<endl; return;}
	
	spotrf_(&UPLO, &N, A->data, &LDA, &INFO);
	if (INFO!=0) {cout<<"Cholesky decomposition unsuccessful in lapack_cholesky_decomp."<<endl; return;}	
	
	return;
}

//cholesky decomposition, A is distroyed
void lapack_cholesky_decomp (gsl_matrix *A)
{
	int N=A->size1, LDA=A->size1, INFO;
	char UPLO='L';
	
	if (N!=(int)A->size2) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_decomp."<<endl; return;}
	
	dpotrf_(&UPLO, &N, A->data, &LDA, &INFO);
	if (INFO!=0) {cout<<"Cholesky decomposition unsuccessful in lapack_cholesky_decomp."<<endl; return;}	
	
	return;
}

//cholesky solve, A is decomposed, 
void lapack_float_cholesky_solve (gsl_matrix_float *A, const gsl_vector_float *b, gsl_vector_float *x)
{
	int N=A->size1, NRHS=1, LDA=A->size1, LDB=b->size, INFO;
	char UPLO='L';
	
	if (N!=(int)A->size2 || N!=LDB) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_solve."<<endl; return;}
	
	gsl_vector_float_memcpy (x, b);
	spotrs_(&UPLO, &N, &NRHS, A->data, &LDA, x->data, &LDB, &INFO);
	if (INFO!=0) {cout<<"Cholesky solve unsuccessful in lapack_cholesky_solve."<<endl; return;}	
	
	return;
}

//cholesky solve, A is decomposed, 
void lapack_cholesky_solve (gsl_matrix *A, const gsl_vector *b, gsl_vector *x)
{
	int N=A->size1, NRHS=1, LDA=A->size1, LDB=b->size, INFO;
	char UPLO='L';
	
	if (N!=(int)A->size2 || N!=LDB) {cout<<"Matrix needs to be symmetric and same dimension in lapack_cholesky_solve."<<endl; return;}
	
	gsl_vector_memcpy (x, b);
	dpotrs_(&UPLO, &N, &NRHS, A->data, &LDA, x->data, &LDB, &INFO);
	if (INFO!=0) {cout<<"Cholesky solve unsuccessful in lapack_cholesky_solve."<<endl; return;}	
	
	return;
}


void lapack_sgemm (char *TransA, char *TransB, float alpha, const gsl_matrix_float *A, const gsl_matrix_float *B, float beta, gsl_matrix_float *C)
{
	int M, N, K1, K2, LDA=A->size1, LDB=B->size1, LDC=C->size2;
	
	if (*TransA=='N' || *TransA=='n') {M=A->size1; K1=A->size2;}
	else if (*TransA=='T' || *TransA=='t') {M=A->size2; K1=A->size1;}
	else {cout<<"need 'N' or 'T' in lapack_dgemm"<<endl; return;}
	
	if (*TransB=='N' || *TransB=='n') {N=B->size2; K2=B->size1;}
	else if (*TransB=='T' || *TransB=='t')  {N=B->size1; K2=B->size2;}
	else {cout<<"need 'N' or 'T' in lapack_dgemm"<<endl;  return;}
	
	if (K1!=K2) {cout<<"A and B not compatible in lapack_dgemm"<<endl; return;}
	if (C->size1!=(size_t)M || C->size2!=(size_t)N) {cout<<"C not compatible in lapack_dgemm"<<endl; return;}
	
	gsl_matrix_float *A_t=gsl_matrix_float_alloc (A->size2, A->size1);
	gsl_matrix_float_transpose_memcpy (A_t, A);
	gsl_matrix_float *B_t=gsl_matrix_float_alloc (B->size2, B->size1);
	gsl_matrix_float_transpose_memcpy (B_t, B);
	gsl_matrix_float *C_t=gsl_matrix_float_alloc (C->size2, C->size1);
	gsl_matrix_float_transpose_memcpy (C_t, C);
	
	sgemm_(TransA, TransB, &M, &N, &K1, &alpha, A_t->data, &LDA, B_t->data, &LDB, &beta, C_t->data, &LDC);
	gsl_matrix_float_transpose_memcpy (C, C_t);
	
	gsl_matrix_float_free (A_t);
	gsl_matrix_float_free (B_t);
	gsl_matrix_float_free (C_t);
	return;
}



void lapack_dgemm (char *TransA, char *TransB, double alpha, const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C)
{
	int M, N, K1, K2, LDA=A->size1, LDB=B->size1, LDC=C->size2;
	
	if (*TransA=='N' || *TransA=='n') {M=A->size1; K1=A->size2;}
	else if (*TransA=='T' || *TransA=='t') {M=A->size2; K1=A->size1;}
	else {cout<<"need 'N' or 'T' in lapack_dgemm"<<endl; return;}
	
	if (*TransB=='N' || *TransB=='n') {N=B->size2; K2=B->size1;}
	else if (*TransB=='T' || *TransB=='t')  {N=B->size1; K2=B->size2;}
	else {cout<<"need 'N' or 'T' in lapack_dgemm"<<endl;  return;}
	
	if (K1!=K2) {cout<<"A and B not compatible in lapack_dgemm"<<endl; return;}
	if (C->size1!=(size_t)M || C->size2!=(size_t)N) {cout<<"C not compatible in lapack_dgemm"<<endl; return;}
	
	gsl_matrix *A_t=gsl_matrix_alloc (A->size2, A->size1);
	gsl_matrix_transpose_memcpy (A_t, A);
	gsl_matrix *B_t=gsl_matrix_alloc (B->size2, B->size1);
	gsl_matrix_transpose_memcpy (B_t, B);
	gsl_matrix *C_t=gsl_matrix_alloc (C->size2, C->size1);
	gsl_matrix_transpose_memcpy (C_t, C);
	
	dgemm_(TransA, TransB, &M, &N, &K1, &alpha, A_t->data, &LDA, B_t->data, &LDB, &beta, C_t->data, &LDC);
	gsl_matrix_transpose_memcpy (C, C_t);
	
	gsl_matrix_free (A_t);
	gsl_matrix_free (B_t);
	gsl_matrix_free (C_t);
	return;
}



//eigen value decomposition, matrix A is destroyed, float seems to have problem with large matrices (in mac)
void lapack_float_eigen_symmv (gsl_matrix_float *A, gsl_vector_float *eval, gsl_matrix_float *evec, const size_t flag_largematrix)
{
	if (flag_largematrix==1) {
		int N=A->size1, LDA=A->size1, INFO, LWORK=-1;
		char JOBZ='V', UPLO='L';
				
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_eigen_symmv."<<endl; return;}
		
		//	float temp[1];
		//	ssyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, temp, &LWORK, &INFO);
		//	if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_eigen_symmv."<<endl; return;}
		//	LWORK=(int)temp[0];
		
		LWORK=3*N;
		float *WORK=new float [LWORK];	
		ssyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, WORK, &LWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_eigen_symmv."<<endl; return;}
		
		gsl_matrix_float_view A_sub=gsl_matrix_float_submatrix(A, 0, 0, N, N);
		gsl_matrix_float_memcpy (evec, &A_sub.matrix);
		gsl_matrix_float_transpose (evec);
		
		delete [] WORK;
	} else {	
		int N=A->size1, LDA=A->size1, LDZ=A->size1, INFO, LWORK=-1, LIWORK=-1;
		char JOBZ='V', UPLO='L', RANGE='A';
		float ABSTOL=1.0E-7;
		
		//VL, VU, IL, IU are not referenced; M equals N if RANGE='A'
		float VL=0.0, VU=0.0;
		int IL=0, IU=0, M;
		
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_float_eigen_symmv."<<endl; return;}
		
		int *ISUPPZ=new int [2*N];
				
		float WORK_temp[1];
		int IWORK_temp[1];
		ssyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK_temp, &LWORK, IWORK_temp, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_float_eigen_symmv."<<endl; return;}	
		LWORK=(int)WORK_temp[0]; LIWORK=(int)IWORK_temp[0];	
		 
		//LWORK=26*N;
		//LIWORK=10*N;
		float *WORK=new float [LWORK];
		int *IWORK=new int [LIWORK];
		
		ssyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_float_eigen_symmv."<<endl; return;}
		
		gsl_matrix_float_transpose (evec);
		
		delete [] ISUPPZ;
		delete [] WORK;
		delete [] IWORK;
	}
	
	
	return;
}



//eigen value decomposition, matrix A is destroyed
void lapack_eigen_symmv (gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec, const size_t flag_largematrix)
{
	if (flag_largematrix==1) {
		int N=A->size1, LDA=A->size1, INFO, LWORK=-1;
		char JOBZ='V', UPLO='L';		
		
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_eigen_symmv."<<endl; return;}
		
		//	double temp[1];
		//	dsyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, temp, &LWORK, &INFO);
		//	if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_eigen_symmv."<<endl; return;}		
		//	LWORK=(int)temp[0];
		
		LWORK=3*N;
		double *WORK=new double [LWORK];	
		dsyev_(&JOBZ, &UPLO, &N, A->data, &LDA, eval->data, WORK, &LWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_eigen_symmv."<<endl; return;}
		
		gsl_matrix_view A_sub=gsl_matrix_submatrix(A, 0, 0, N, N);
		gsl_matrix_memcpy (evec, &A_sub.matrix);
		gsl_matrix_transpose (evec);
		
		delete [] WORK;
	} else {	
		int N=A->size1, LDA=A->size1, LDZ=A->size1, INFO, LWORK=-1, LIWORK=-1;
		char JOBZ='V', UPLO='L', RANGE='A';
		double ABSTOL=1.0E-7;
		
		//VL, VU, IL, IU are not referenced; M equals N if RANGE='A'
		double VL=0.0, VU=0.0;
		int IL=0, IU=0, M;
		
		if (N!=(int)A->size2 || N!=(int)eval->size) {cout<<"Matrix needs to be symmetric and same dimension in lapack_eigen_symmv."<<endl; return;}
		
		int *ISUPPZ=new int [2*N];
		
		double WORK_temp[1];
		int IWORK_temp[1];

		dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK_temp, &LWORK, IWORK_temp, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Work space estimate unsuccessful in lapack_eigen_symmv."<<endl; return;}	
		LWORK=(int)WORK_temp[0]; LIWORK=(int)IWORK_temp[0];	

		//LWORK=26*N;
		//LIWORK=10*N;
		double *WORK=new double [LWORK];
		int *IWORK=new int [LIWORK];
		
		dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A->data, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, eval->data, evec->data, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
		if (INFO!=0) {cout<<"Eigen decomposition unsuccessful in lapack_eigen_symmv."<<endl; return;}

		gsl_matrix_transpose (evec);
		
		delete [] ISUPPZ;
		delete [] WORK;
		delete [] IWORK;
	}
	
	return;
}

//DO NOT set eigen values to be positive
double EigenDecomp (gsl_matrix *G, gsl_matrix *U, gsl_vector *eval, const size_t flag_largematrix)
{
#ifdef WITH_LAPACK
	lapack_eigen_symmv (G, eval, U, flag_largematrix);
#else
	gsl_eigen_symmv_workspace *w=gsl_eigen_symmv_alloc (G->size1);
	gsl_eigen_symmv (G, eval, U, w);
	gsl_eigen_symmv_free (w);	
#endif	
	/*
	for (size_t i=0; i<eval->size; ++i) {
		if (gsl_vector_get (eval, i)<1e-10) {
//			cout<<gsl_vector_get (eval, i)<<endl;
			gsl_vector_set (eval, i, 0);			
		}
	}
	*/
	//calculate track_G=mean(diag(G))	
	double d=0.0;
	for (size_t i=0; i<eval->size; ++i) {
		d+=gsl_vector_get(eval, i);
	}
	d/=(double)eval->size;
	
	return d;
}


//DO NOT set eigen values to be positive
double EigenDecomp (gsl_matrix_float *G, gsl_matrix_float *U, gsl_vector_float *eval, const size_t flag_largematrix)
{
#ifdef WITH_LAPACK
	lapack_float_eigen_symmv (G, eval, U, flag_largematrix);
#else
	//gsl doesn't provide float precision eigen decomposition; plus, float precision eigen decomposition in lapack may not work on OS 10.4
	//first change to double precision
	gsl_matrix *G_double=gsl_matrix_alloc (G->size1, G->size2);
	gsl_matrix *U_double=gsl_matrix_alloc (U->size1, U->size2);
	gsl_vector *eval_double=gsl_vector_alloc (eval->size);
	for (size_t i=0; i<G->size1; i++) {
		for (size_t j=0; j<G->size2; j++) {
			gsl_matrix_set(G_double, i, j, gsl_matrix_float_get(G, i, j));
		}
	}	
	gsl_eigen_symmv_workspace *w_space=gsl_eigen_symmv_alloc (G->size1);
	gsl_eigen_symmv (G_double, eval_double, U_double, w_space);
	gsl_eigen_symmv_free (w_space);	
	
	//change back to float precision
	for (size_t i=0; i<G->size1; i++) {
		for (size_t j=0; j<G->size2; j++) {
			gsl_matrix_float_set(K, i, j, gsl_matrix_get(G_double, i, j));
		}
	}
	for (size_t i=0; i<U->size1; i++) {
		for (size_t j=0; j<U->size2; j++) {
			gsl_matrix_float_set(U, i, j, gsl_matrix_get(U_double, i, j));
		}
	}
	for (size_t i=0; i<eval->size; i++) {
		gsl_vector_float_set(eval, i, gsl_vector_get(eval_double, i));
	}	
	
	//delete double precision matrices
	gsl_matrix_free (G_double);
	gsl_matrix_free (U_double);
	gsl_vector_free (eval_double);
#endif
	/*
	for (size_t i=0; i<eval->size; ++i) {
		if (gsl_vector_float_get (eval, i)<1e-10) {
			gsl_vector_float_set (eval, i, 0);
		}
	}
	*/
	//calculate track_G=mean(diag(G))	
	double d=0.0;
	for (size_t i=0; i<eval->size; ++i) {
		d+=gsl_vector_float_get(eval, i);
	}
	d/=(double)eval->size;
	
	return d;
}


double CholeskySolve(gsl_matrix *Omega, gsl_vector *Xty, gsl_vector *OiXty)
{
	double logdet_O=0.0;
	
#ifdef WITH_LAPACK
	lapack_cholesky_decomp(Omega);
	for (size_t i=0; i<Omega->size1; ++i) {
		logdet_O+=log(gsl_matrix_get (Omega, i, i));
	}	
	logdet_O*=2.0;	
	lapack_cholesky_solve(Omega, Xty, OiXty);	
#else	
	int status = gsl_linalg_cholesky_decomp(Omega);
	if(status == GSL_EDOM) {
		cout << "## non-positive definite matrix" << endl; 
		//		exit(0); 
	}
	
	for (size_t i=0; i<Omega->size1; ++i) {
		logdet_O+=log(gsl_matrix_get (Omega, i, i));
	}
	logdet_O*=2.0;	
	
	gsl_vector_memcpy (OiXty, Xty);
	gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, Omega, OiXty); 
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, OiXty); 	
	//	gsl_linalg_cholesky_solve(XtX, Xty, iXty);
#endif
	
	return logdet_O;
}


double CholeskySolve(gsl_matrix_float *Omega, gsl_vector_float *Xty, gsl_vector_float *OiXty)
{
	double logdet_O=0.0;
	
#ifdef WITH_LAPACK
	lapack_float_cholesky_decomp(Omega);
	for (size_t i=0; i<Omega->size1; ++i) {
		logdet_O+=log(gsl_matrix_float_get (Omega, i, i));
	}	
	logdet_O*=2.0;	
	lapack_float_cholesky_solve(Omega, Xty, OiXty);	
#else
	gsl_matrix *Omega_double=gsl_matrix_alloc (Omega->size1, Omega->size2);
	double d;
	for (size_t i=0; i<Omega->size1; ++i) {
		for (size_t j=0; j<Omega->size2; ++j) {
			d=(double)gsl_matrix_float_get (Omega, i, j);
			gsl_matrix_set (Omega_double, i, j, d);
		}
	}
	
	int status = gsl_linalg_cholesky_decomp(Omega_double);
	if(status == GSL_EDOM) {
		cout << "## non-positive definite matrix" << endl; 
		//		exit(0); 
	}	
	
	for (size_t i=0; i<Omega->size1; ++i) {
		for (size_t j=0; j<Omega->size2; ++j) {
			d=gsl_matrix_get (Omega_double, i, j);
			if (j==i) {logdet_O+=log(d);}
			gsl_matrix_float_set (Omega, i, j, (float)d);
		}
	}
	logdet_O*=2.0;	
	
	gsl_vector_float_memcpy (OiXty, Xty);
	gsl_blas_strsv(CblasLower, CblasNoTrans, CblasNonUnit, Omega, OiXty); 
	gsl_blas_strsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, OiXty); 	
	//	gsl_linalg_cholesky_solve(XtX, Xty, iXty);
	
	gsl_matrix_free (Omega_double);
#endif
	
	return logdet_O;
}	


//LU decomposition
void LUDecomp (gsl_matrix *LU, gsl_permutation *p, int *signum)
{
	gsl_linalg_LU_decomp (LU, p, signum);
	return;
}

void LUDecomp (gsl_matrix_float *LU, gsl_permutation *p, int *signum)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	//LU decomposition
	gsl_linalg_LU_decomp (LU_double, p, signum);
	
	//copy float matrix to double
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_float_set (LU, i, j, gsl_matrix_get(LU_double, i, j));
		}
	}
	
	//free matrix
	gsl_matrix_free (LU_double);
	return;
}


//LU invert
void LUInvert (const gsl_matrix *LU, const gsl_permutation *p, gsl_matrix *inverse)
{
	gsl_linalg_LU_invert (LU, p, inverse);
	return;
}

void LUInvert (const gsl_matrix_float *LU, const gsl_permutation *p, gsl_matrix_float *inverse)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	gsl_matrix *inverse_double=gsl_matrix_alloc (inverse->size1, inverse->size2);
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	//LU decomposition
	gsl_linalg_LU_invert (LU_double, p, inverse_double);
	
	//copy float matrix to double
	for (size_t i=0; i<inverse->size1; i++) {
		for (size_t j=0; j<inverse->size2; j++) {
			gsl_matrix_float_set (inverse, i, j, gsl_matrix_get(inverse_double, i, j));
		}
	}
	
	//free matrix
	gsl_matrix_free (LU_double);
	gsl_matrix_free (inverse_double);
	return;
}

//LU lndet
double LULndet (gsl_matrix *LU)
{
	double d;
	d=gsl_linalg_LU_lndet (LU);
	return d;
}

double LULndet (gsl_matrix_float *LU)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	double d;
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	//LU decomposition
	d=gsl_linalg_LU_lndet (LU_double);
	
	//copy float matrix to double
	/*
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_float_set (LU, i, j, gsl_matrix_get(LU_double, i, j));
		}
	}
	*/
	//free matrix
	gsl_matrix_free (LU_double);
	return d;
}


//LU solve
void LUSolve (const gsl_matrix *LU, const gsl_permutation *p, const gsl_vector *b, gsl_vector *x)
{
	gsl_linalg_LU_solve (LU, p, b, x);
	return;
}

void LUSolve (const gsl_matrix_float *LU, const gsl_permutation *p, const gsl_vector_float *b, gsl_vector_float *x)
{
	gsl_matrix *LU_double=gsl_matrix_alloc (LU->size1, LU->size2);
	gsl_vector *b_double=gsl_vector_alloc (b->size);
	gsl_vector *x_double=gsl_vector_alloc (x->size);	
	
	//copy float matrix to double	
	for (size_t i=0; i<LU->size1; i++) {
		for (size_t j=0; j<LU->size2; j++) {
			gsl_matrix_set (LU_double, i, j, gsl_matrix_float_get(LU, i, j));
		}
	}
	
	for (size_t i=0; i<b->size; i++) {
		gsl_vector_set (b_double, i, gsl_vector_float_get(b, i));
	}
	
	for (size_t i=0; i<x->size; i++) {
		gsl_vector_set (x_double, i, gsl_vector_float_get(x, i));
	}
	
	//LU decomposition
	gsl_linalg_LU_solve (LU_double, p, b_double, x_double);
	
	//copy float matrix to double
	for (size_t i=0; i<x->size; i++) {
		gsl_vector_float_set (x, i, gsl_vector_get(x_double, i));
	}
	
	//free matrix
	gsl_matrix_free (LU_double);
	gsl_vector_free (b_double);
	gsl_vector_free (x_double);
	return;
}


