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

	n_trials = ntrials;	

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
	
	/* Being with family = Bernouli */
	cp3(eta, n, 1, ntrials, mu); 


	/*if(family.mcml$family.glmm=="bernoulli.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta),as.integer(1),as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}
	if(family.mcml$family.glmm=="poisson.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta), as.integer(2), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}
	if(family.mcml$family.glmm=="binomial.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta), as.integer(3), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}*/


	double value = elc(y, x, n, nbeta, eta, 1,  n_trials, 1.0, nbeta, nbeta*nbeta) - 0.5*matvecmult(s, nbeta, s, nbeta, val);

	
	/*gradient calculation */
	
	matvecmult(x, y-mu, n, n, db);
	matvecmult(a, z, n, n, az);
	matvecmult(az, y-mu, n, n, azymu);
	double ds = azymu-s;
	double gradient [2] = (db, ds);

	/*
	#hessian calculation
	if(family.mcml$family.glmm=="bernoulli.glmm") {cdub<-family.mcml$cpp(eta)}
	if(family.mcml$family.glmm=="poisson.glmm"){cdub<-family.mcml$cpp(eta)}
	if(family.mcml$family.glmm=="binomial.glmm"){cdub<-family.mcml$cpp(eta, ntrials)}
	cdub<-as.vector(cdub)
	cdub<-diag(cdub)
	kyoo<-nrow(A)
	piece1<- (-t(X)%*%cdub%*%X)
	piece2<- (-t(X)%*%cdub%*%Z%*%A)
	piece3<- (-A%*%t(Z)%*%cdub%*%X)
	piece4<- (-A%*%t(Z)%*%cdub%*%Z%*%A -diag(kyoo))
	hessian<- rbind(cbind(piece1,piece2),cbind(piece3,piece4))

	cache$s.twid<-s
	cache$beta.twid<-beta
	
	list(value=value,gradient=gradient,hessian=hessian)
	*/

	/*Bernoulli */
	cpp3(eta, n, 1, ntrials, cdub);
	double cdub_mat [ 
}



