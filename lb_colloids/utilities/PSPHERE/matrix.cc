// matrix.c : implementation file
//
//#include "matrix.h"


#include <stdio.h>
#include <malloc.h>

#include "matrix.h"
#include "nrutil.h"


MATRIXM::MATRIXM(int nrows, int ncols)
{
	
	int i;

	rows=nrows;
	cols=ncols;
	matrix=new double*[(unsigned) (rows)];
	if (!matrix) nrerror("allocation failure 1 in matrix()");
	matrix -= 1;
 	for(i=1;i<=rows;i++) {
		matrix[i]=new double[(unsigned) (cols)];
		if (!matrix[i]) nrerror("allocation failure 2 in matrix()");
		matrix[i] -= 1;
	}
}


MATRIXM::~MATRIXM()
{
	int i;
	int nrl,ncl;	
    
	nrl=1;
	ncl=1;

	for(i=rows;i>=nrl;i--){
		delete [] (matrix[i]+ncl);
	}

	delete [] (matrix+nrl);

}
/*-----------------------------------------------------------------------*/

VECTORM::VECTORM(int nrows)
{
	rows=nrows;
	vector=new double[(unsigned) (rows)];
	if (!vector) nrerror("allocation failure in vector()");
	vector-=1;
}

VECTORM::~VECTORM()
{

	delete [] (vector+1);

}
/*-----------------------------------------------------------------------*/
MATRIXM *read_matrix(FILE *binfile)
{
	int i,j;
	MATRIXM *m;
	int rows,cols;


	fread(&rows,sizeof(int),1,binfile);
	fread(&cols,sizeof(int),1,binfile);

	m=new MATRIXM(rows,cols);
	m->rows=rows;
	m->cols=cols;

	for (j=1;j<=m->cols;j++){
		for (i=1;i<=m->rows;i++){
			fread(&m->matrix[i][j],sizeof(double),1,binfile);

		}
	}
	return(m);
}
/*-----------------------------------------------------------------------*/
VECTORM *read_vector(FILE *binfile)
{
	int i;
	VECTORM *v;
	int rows,cols;

	fread(&rows,sizeof(int),1,binfile);
	fread(&cols,sizeof(int),1,binfile);
/*	printf("%d %d\n",rows,cols);*/
	if (cols!=1) {
		fprintf(stderr,"%d %d\n",rows,cols);
		nrerror("cols not equal to 1 in read_vector");
	}

	v=new VECTORM(rows);

	for (i=1;i<=v->rows;i++){
		fread(&v->vector[i],sizeof(double),1,binfile);
	}
	return(v);
}

/*-----------------------------------------------------------------------*/
int matrix_mul(MATRIXM *m1, MATRIXM *m2, MATRIXM *m3)
{
	int i,j,k;

	if (m1->cols!=m2->rows) nrerror("matrix_mul 1, cols&rows");
	if (m1->rows!=m3->rows) nrerror("matrix_mul 2, cols&rows");
	if (m2->cols!=m3->cols) nrerror("matrix_mul 3, cols&rows");

	for (i=1;i<=m2->cols;i++){
		for (j=1;j<=m1->rows;j++){
			m3->matrix[j][i]=0.0;
			for (k=1;k<=m2->cols;k++){
				m3->matrix[j][i]+=m1->matrix[j][k]*m2->matrix[k][i];
			}
		}

	}
	return(0);
}
/*-----------------------------------------------------------------------*/
int matrix_add(MATRIXM *m1, MATRIXM *m2, MATRIXM *m3)
{
	int i,j;

	if ((m1->cols!=m2->cols) || (m1->cols!=m3->cols)) nrerror("matrix_add 1, cols");
	if ((m1->rows!=m2->rows) || (m1->rows!=m3->rows)) nrerror("matrix_add 2, rows");

	for (i=1;i<=m2->cols;i++){
		for (j=1;j<=m1->rows;j++){
			m3->matrix[j][i]=m1->matrix[j][i]+m2->matrix[j][i];
		}

	}
	return(0);
}
/*-----------------------------------------------------------------------*/
int matrix_vector_mul(MATRIXM *m1,VECTORM *v1,VECTORM *v2)
{
	int j,k;

	if (m1->cols!=v1->rows) nrerror("matrix_vector_mul 1, cols&rows");
	if (m1->rows!=v2->rows) nrerror("matrix_vector_mul 2, cols&rows");

	for (j=1;j<=m1->rows;j++){
		v2->vector[j]=0.0;
		for (k=1;k<=v1->rows;k++){
			v2->vector[j]+=m1->matrix[j][k]*v1->vector[k];
		}
	}
	return(0);
}
/*-----------------------------------------------------------------------*/

int vector_add(VECTORM *v1,VECTORM *v2,VECTORM *v3)
{
	int j;

	if ((v1->rows!=v2->rows)||(v1->rows!=v3->rows)) nrerror("vector_add 1, rows");

	for (j=1;j<=v1->rows;j++){
		v3->vector[j]=v1->vector[j]+v2->vector[j];
	}
	return(0);
}
