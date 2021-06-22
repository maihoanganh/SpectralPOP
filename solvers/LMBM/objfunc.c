
#include "objfunc.h"

#include <string.h>
#include <stdio.h>

// static char *objfunc_name = NULL;


static double (*objfunptr)(int,double*,double*);
void objfunc_(int *n, double *x, double *f, double *g)
{
  *f=objfunptr(*n,x,g);
}

void set_obj_func(void* fptr)
{
  objfunptr = fptr;
}
