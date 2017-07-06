/* NRUTIL.C */
    
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "nrutil.h"
    
void nrerror(char *error_text)
{
  //	void exit();
  
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");

	exit(0);
}
    
    
