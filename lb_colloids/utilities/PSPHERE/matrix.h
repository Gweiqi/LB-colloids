#ifndef MATRIX_000
#define MATRIX_000




/*typedef struct
	{
	int rows;
	int cols;
	double **matrix;
	} MATRIXM;

typedef struct 
	{
	int rows;
	double *vector;
	} VECTORM;
*/
class MATRIXM
	{
	public:
		int rows;
		int cols;
		double **matrix;
	public:
		MATRIXM(int,int);
		~MATRIXM();
	};

class VECTORM 
	{
	public:
		int rows;
		double *vector;
	public:
		VECTORM(int);
		~VECTORM();
	};


//VECTORM 	*alloc_vector(int);
//void 	free_vector(VECTORM *);
//MATRIXM 	*alloc_matrix(int,int);
//void 	free_matrix(MATRIXM *);
int 	matrix_mul(MATRIXM *,MATRIXM *,MATRIXM *);
int 	matrix_add(MATRIXM *,MATRIXM *,MATRIXM *);
int 	matrix_vector_mul(MATRIXM *,VECTORM *,VECTORM *);
int 	vector_add(VECTORM *,VECTORM *,VECTORM *);
MATRIXM	*read_matrix(FILE *);
VECTORM	*read_vector(FILE *);


#endif
