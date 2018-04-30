#include "myheader.h"

/* fn.outer */


/* fn.inner.trust */
function(mypar,Y,X,Z,A,family.mcml,nbeta,cache, ntrials )
void fnInnerTrust (double *mypar, double *Y, double *X, double *Z, double *A, int nbeta, double *cache, int ntrials)
{
	int beta [nbeta];
	int i;
	for (i = 0; i < nbeta;  i++){
		beta[i] = i;
	}
	
	int s [nbeta]
	for (i = 0; i < nbeta; i++){
		s[i] = -i;
	}

	
	/*eta<-X%*%beta+Z%*%A%*%s*/
	double *x = Calloc(*X, double);
	double *y = Calloc(*Y, double);
	double *z = Calloc(*Z, double);
	double *a = Calloc(*A, double);
	
	matvecmult(x, beta, n, nbeta, xbeta);
	matvecmult(z, a, n, n, ZA);
	matvecmult(ZA, s, n, nbeta, ZAs);
	addvec(xbeta, ZAs, n, eta);

	/*family.mcml<-getFamily(family.mcml)
	if(family.mcml$family.glmm=="bernoulli.glmm") {mu<-family.mcml$cp(eta)}
	if(family.mcml$family.glmm=="poisson.glmm"){mu<-family.mcml$cp(eta)}
	if(family.mcml$family.glmm=="binomial.glmm"){mu<-family.mcml$cp(eta, ntrials)}*/
	
	if (family.mcml = "bernoulli.glmm"){
		mu = family.mcml

	/*if(family.mcml$family.glmm=="bernoulli.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta),as.integer(1),as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}
	if(family.mcml$family.glmm=="poisson.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta), as.integer(2), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}
	if(family.mcml$family.glmm=="binomial.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta), as.integer(3), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}*/



}

