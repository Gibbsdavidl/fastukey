#include "stdlib.h"
#include "stdio.h"
//#include "R.h"	
//#include "Rmath.h"
#include "math.h"


// vector_plus_const //
// takes a vector and a constant 
// and returns new memory with result.
double* vpc (double* v, double c, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v[i] + c;
  }
  return(result);
}


// vector_minus_const //
// takes a vector and a constant 
// and returns new memory with result.
double* vmc (double* v, double c, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v[i] - c;
  }
  return(result);
}


// vector_plus_vector //
// takes a vector and a vector
// and returns new memory with result.
double* vpv (double* v1, double* v2, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v1[i] + v2[i];
  }
  return(result);
}


// vector_minus_vector //
// takes a vector and a vector
// and returns new memory with result.
double* vmv (double* v1, double* v2, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v1[i] - v2[i];
  }
  return(result);
}


// const_minus_vector //
// takes two vectors and returns a new one //
double* cmv (double c, double* v, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = c - v[i];
  }
  return(result);
}


// vector_times_vector //
// takes two vectors and returns a new one //
double* vtv (double* v1, double* v2, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v1[i] * v2[i];
  }
  return(result);
}


// vector_times_const //
// takes two vectors and returns a new one //
double* vtc (double* v, double c, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v[i] * c;
  }
  return(result);
}


// vector_divided_by_const //
// takes two vectors and returns a new one //
double* vdc (double* v, double c, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v[i] / c;
  }
  return(result);
}


// vector_sum //
// takes one vector and returns the sum of elements //
double vsum (double* v, int n) 
{
  int i = 0;
  double result = 0.0;
  for(i = 0; i < n; i++) {
    result += v[i];
  }
  return(result);
}


// square each element of a vector //
double* vpow (double* v, int e, int n) {
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = pow(v[i], e);
  }
  return(result);
}


double mean(double* v, int n) {
  int i = 0;
  double sum = 0;
  double nn = (double)n;
  for (i = 0; i < n; ++i) {
    sum += v[i];
  }
  return(sum/nn);
}


double* rhobwC (double* x, double c1, int n){    
  double* res0, *res1, *res2, *res3, *res4, *res5;
  double* res6, *res7, *res8, *res9, *res10, *res11;
  double* ivec = malloc(n * sizeof(double));
  int i = 0;
  for(i=0; i<n; ++i) {
    if (fabs(x[i]) > c1) {
      ivec[i] = 1.0;
    } else {
      ivec[i] = 0;
    }
  }
  // watch for memory leaks!  each call to a vec function returns a vector!!
  res0 = vtc(ivec, (pow(c1,2)/6.0), n);      // (c1^2/6)*ivec
  res1 = cmv(1, ivec, n);                    // (1-ivec)
  res6 = vpow(x,2,n);                        // x^2
  res2 = vdc(res6,2.0,n);                    // (x^2/2)
  res7 = vpow(x,4,n);                        // (x^4)
  res3 = vdc(res7, (2*pow(c1,2)), n);        // (x^4)/(2*c1^2)
  res8 = vpow(x,6,n);                        // x^6
  res4 = vdc(res8, (6*pow(c1,4)), n);        // (x^6)/(6*c1^4) 
  res9 = vpv(res3, res4, n);                 // (x^4)/(2*c1^2) + (x^6)/(6*c1^4)
  res10 = vmv(res2, res9, n);                // (x^2/2) - res9
  res11 = vtv(res1, res10, n);               // (1-ivec) * res10
  res5 = vpv(res0, res11, n);                // (c1^2/6)*ivec + res11
  free(res0); free(res1); free(res2); free(res3); free(res4); free(res6);
  free(res7); free(res8); free(res9); free(res10); free(res11); free(ivec);
  return(res5);
  //return((c1^2/6)*ivec + (1-ivec)*(x^2/2-x^4/(2*c1^2)+x^6/(6*c1^4)));
}


double* psibwC (double* x, double c1, double n){
  double *res0, *res1, *res2, *res3, *res4, *res5, *res6;
  double* ivec = malloc(n * sizeof(double));
  int i = 0;
  for(i=0; i<n; ++i) {
    if (fabs(x[i]) > c1) {
      ivec[i] = 1.0;
    } else {
      ivec[i] = 0.0;
    }
  }
  res0 = cmv(1, ivec, n);         // (1-ivec)
  res1 = vdc(x, c1, n);           // (x/c1)
  res2 = vpow(res1, 2, n);        // (x/c1)^2
  res3 = cmv(1, res2, n);         // 1 - (x/c1)^2
  res4 = vpow(res3, 2, n);        // (1-(x/c1)^2)^2
  res5 = vtv(x, res4, n);         // x * (1-...)^2
  res6 = vtv(res0, res5, n);      // (1-ivec)*(x*(1-...)^2)
  free(res0); free(res1); free(res2); free(ivec);
  free(res3); free(res4); free(res5);
  //return((1-ivec)*(x*(1-(x/c1)^2)^2));
  return(res6);
}

void printit(double* x, int n) {
  int i = 0;
  for (i = 0; i < n; ++i) {
    printf("%f  ", x[i]);
  }
  printf("\n");
}


double* fast_rhobw (double* d, double k, double c1, int n)
{    
  int i = 0;
  double* g1 = malloc(n*sizeof(double));
  double res1, res6, res7, res8, ivec, x;
  double c2 = pow(c1,2)/6.0;
  double c3 = 2*pow(c1,2);
  double c4 = 6*pow(c1,4);
  for (i=0; i<n; ++i) {
    x = d[i]/k;
    if (fabs(x) > c1) {ivec=1.0;} else {ivec=0.0;}
    res1 = (ivec * c2);
    res6 = pow(x,2);
    res7 = pow(x,4);
    res8 = pow(x,6);
    g1[i] = res1 + (1-ivec)*((res6/2) - (res7/c3) + (res8/c4)); // res10
  }
  return(g1);
  //return((c1^2/6)*ivec + (1-ivec)*(x^2/2-x^4/(2*c1^2)+x^6/(6*c1^4)));
}


double* fast_psibw (double* d, double k, double c1, double n){
  double res0, res1, res2, res3,
    res4, res5, res6, ivec, x, x2;
  int i = 0;
  double* res = malloc(n*sizeof(double));
  for (i=0; i<n; ++i) {
    x = d[i]/k;
    x2 = d[i]/pow(k,2); // memory
    if (fabs(x) > c1) {ivec=1.0;} else {ivec=0.0;}
    res0 = 1 - pow((x/c1), 2);
    res1 = x * pow(res0, 2);
    res[i] = (1-ivec) * x2 * res1; 
  }
  //return((1-ivec)*(x*(1-(x/c1)^2)^2));
  return(res);
}


void fastksolvec (double* d, double* p,
		  double* c1, double* b0,
		  int* n, double* output) {
  int i = 0;
  double  k = 1;
  double  kold = 1;
  int     iter = 1;
  double  crit = 100;
  double  eps = 0.0000000001;
  double* r = NULL;
  double* psi = NULL;
  double  fk = 0;
  double  fkp = 0;
  while ((crit > eps)&&(iter<100)){
    kold = k;
    //printf("\nk: %f\n", k);
    r = fast_rhobw(d,k,*c1,*n); // memory
    //printit(r, *n);
    fk = mean(r, *n) - *b0;
    psi = fast_psibw(d,k,*c1,*n);  // memory
    //printit(psi, *n);
    fkp = -mean(psi, *n);
    //printf("fk: %f   fkp: %f \n", fk, fkp);
    if (fkp==0) {
      output = NULL;
      printf("no values close enough!\n");
      return;
    }
    k = k - fk/fkp;
    if (k < 0) {
      k = kold/2;
    }
    crit = fabs(k-kold);
    iter += 1;
    free(r); free(psi); 
  }
  output[0] = k;
  return;
}


int main(void) {

  double c1 = 5.06883;
  double b0 = 0.8564345;
  double p = 2;
  double d[] = {0.2724219,0.5817435,0.4753930,1.7951429,
		0.7584367,0.6522371,1.2774909,2.2665531,
		1.2767456,1.0847021};
  int n = 10;
  double* o;
  o = malloc(sizeof(double));
  fastksolvec (d, &p,
	       &c1, &b0,
	       &n, o);  
  printf("%f\n", *o);
  return(0);
}


