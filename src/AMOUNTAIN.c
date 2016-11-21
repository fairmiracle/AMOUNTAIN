#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_cblas.h>
#define epsilon 1e-6
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
//gcc -c AMOUNTAIN.c -fPIC -std=c99;gcc -shared -o AMOUNTAIN.so AMOUNTAIN.o -lgsl -lgslcblas
double Sum(double *list, int length) {
    double sum = 0;
    for (int i = 0; i < length; i++) {
        sum += list[i];
    }
    return sum;
}

void Vmul(double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        x[i] = x[i]*y[i];
    }
}
/*
void OutPutVector(double* A, int N)
{
	//printf("Vector: %d\n",N);
	int i;	
	for(i = 0; i < N; i++)
			printf("%10.4f ",*(A+i));
	printf("\n");
}
*/

//S = which(list > theta)
//list[S]...
int Which(double *list, int n, double theta, double *cplist) {

	int cpn = 0;
    for (int i = 0; i < n; i++) {
    	if (list[i] > theta){
    		cpn += 1;
    		cplist[i] = list[i];
    	}
    }
    return cpn;
}

//S = which(list > theta)
//sum(y[S]) and sum(y[S]*y[S])
int Which2(double *list, double *y, int n, double theta, double *sum1, double *sum2) {

	//int cpn = 0;
    //for (int i = 0; i < n; i++) {
    //	if (list[i] >= theta){
    //		cpn += 1;
    //		cplist[i] = y[i];
    //	}
    //}
    int cpn = 0;
    *sum1 = 0;
    *sum2 = 0;
    for (int i = 0; i < n; i++) {
    	if (list[i] >= theta){
    		cpn += 1;
    		*sum1 += y[i];
    		*sum2 += y[i]*y[i];
    	}
    }
    return cpn;
}

// argmin_x 0.5*||x-y||^2 s.t. alpha\sum x_i+(1-alpha)*x_2^2=t
// EN: alpha x_i + (1-alpha)*x_2^2=t
// x_i(theta) = max(0,y_i-theta)/(1+(2-2*alpha)*theta)
// linear time refer to Efficient Euclidean Projections via Piecewise Root Finding and Its Application in
// Gradient Projection, rewrite it to normal form
// Increasing alpha means increase sparsity
void epENNORM(double *y, double *m_t, int *m_n, double *x, double *m_alpha){
	double Epsilon = 1e-9;
	double t = *m_t;
	double alpha = *m_alpha;
	int n = *m_n;
	int maxiter = 100;
	double theta = 0;
	memcpy(x, y, n*sizeof(double));
	Vmul(x, x, n);
	
	if(alpha*Sum(y, n) + (1-alpha)*Sum(x, n) <= t){
    	memcpy(x,y,n*sizeof(double));
    	return;
	} 
	else {
		//tmpx = (y-theta*alpha)/(1+(2-2*alpha)*theta);
		double *tmpx = Malloc(double,n);
		for (int i = 0; i < n; ++i)
		{
			tmpx[i] = (y[i]-theta*alpha)/(1+(2-2*alpha)*theta);
		}

		double s1 = 0;
		double s2 = 0;

		int cpn = Which2(tmpx, y, n, 0, &s1, &s2);
    	double a = alpha*alpha*(alpha-1)*cpn - 4*t*(1-alpha)*(1-alpha);
    	double b = -alpha*alpha*cpn - 4*t*(1-alpha);
    	double c = alpha*s1+(1-alpha)*s2-t;
    	for (int i = 0;  i < maxiter; i++) {
        	theta = (-b-sqrt(b*b-4*a*c))/(2*a);
  			//printf("%f,%f,%f,%f,%f\n", a,b,c,theta,(a*theta*theta + b*theta + c));

        	for (int j = 0; j < n; ++j)
			{
				tmpx[j] = (y[j]-theta*alpha)/(1+(2-2*alpha)*theta);
			}

    		cpn = Which2(tmpx, y, n, 0, &s1, &s2);
    		
    		a = alpha*alpha*(alpha-1)*cpn - 4*t*(1-alpha)*(1-alpha);
    		b = -alpha*alpha*cpn - 4*t*(1-alpha);
    		c = alpha*s1+(1-alpha)*s2-t;

    		//printf("%f,%f,%f,%f,%f,%f\n", Sum(tmpx,n),a,b,c,theta,(a*theta*theta + b*theta + c));
        	
    		if ((a*theta*theta + b*theta + c) < Epsilon){
    			break;
    		}
    	}
    	free(tmpx);
	}
	for (int i = 0; i < n; i++) {
    	x[i] = MAX(0,(y[i] - theta*alpha)/(1+(2-2*alpha)*theta));
	}
}
//argmin_x 0.5*||x-y||^2 s.t. \sum x_i+0.5*gamma*x_2=t
// gamma = 2*(1-alpha)/alpha
// EN: alpha x_i + (1-alpha)*x_2=t
// x_i(theta) = max(0,y_i-theta)/(1+theta*gamma)
// linear time refer to Efficient Euclidean Projections via Piecewise Root Finding and Its Application in Gradient Projection
//increase alpha means increase sparsity
void epEN(double *y, double t, int n, double *x, double alpha){
	int maxiter = 100;
	double gamma = 2*(1-alpha)/alpha;
	double theta = 0;
	memcpy(x, y, n*sizeof(double));
	Vmul(x, x, n);
	if( Sum(y, n) + 0.5*gamma*Sum(x, n) <= t){
		memcpy(x, y, n*sizeof(double));
		return;
	} 
	else {
		double *cplist = Malloc(double,n);
		memset(cplist,0,n*sizeof(double));
		int cpn = Which(y, n, theta, cplist);
		double a = t*gamma*gamma + 0.5*gamma*cpn;
		double b = 2*a/gamma;
		double c = - Sum(cplist, cpn);
		Vmul(cplist, cplist, cpn);
		c += t - 0.5*gamma*Sum(cplist, cpn);
		for (int i = 0; i < maxiter; i++) {
			theta = (-b + sqrt(b*b-4*a*c))/(2*a);
			cpn = Which(y, n, theta, cplist);
			a = t*gamma*gamma + 0.5*gamma*cpn;
			b = 2*a/gamma;
			c = - Sum(cplist, cpn);
			Vmul(cplist, cplist, cpn);
			c += t - 0.5*gamma*Sum(cplist, cpn);
			if ((a*theta*theta + b*theta + c) < epsilon){
				break;
			}
		}
		free(cplist);
	}
	for (int i = 0; i < n; i++) {
		x[i] = MAX(0,y[i] - theta)/(1 + theta*gamma);
	}
}

// 10000 nodes single layer,maxiter=1000,20x faster than R, 427.181s vs 19.132s
void miGPFixSS(double *W, double *z, double *x0, int *m_n, double *x, 
			double *func, double *m_a, double *m_lambda, int *m_maxiter){
	int n = *m_n;
	double a = *m_a;
	double lambda = *m_lambda;
	int maxiter = *m_maxiter;
	double t = 1;

	double *mWx = Malloc(double, n);		//-W*x
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, -1, W, n, x, 1, 0, mWx, 1);
	double *grad = Malloc(double, n);
	memcpy(grad,mWx,n*sizeof(double));
	cblas_daxpy(n, -lambda, z, 1, grad, 1);
	double f_x = 0.5*cblas_ddot(n, mWx, 1, x, 1) - lambda*cblas_ddot(n, z, 1, x, 1);
	
	//double *func = Malloc(double,maxiter);
	//memset(func,0,maxiter*sizeof(double));
	double *y = Malloc(double,n);
	double *x_cand = Malloc(double,n);
	double *diffx = Malloc(double,n);
	int iteration = 0;
	for (iteration = 0; iteration < maxiter; iteration++){
		func[iteration] = f_x;
		
		memcpy(y, x, n*sizeof(double));
		cblas_daxpy(n, -1, grad, 1, y, 1);
		
		//epEN(y, 1, n, x_cand, a);
		epENNORM(y,&t,&n,x_cand,&a);
		memcpy(diffx, x_cand, n*sizeof(double));
		cblas_daxpy(n, -1, x, 1, diffx, 1);
		
		if(sqrt(cblas_ddot(n, diffx, 1, diffx, 1)) < epsilon){
			break;
		}

 		memcpy(x, x_cand, n*sizeof(double));
 		cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, -1, W, n, x, 1, 0, mWx, 1);
		memcpy(grad,mWx,n*sizeof(double));
		cblas_daxpy(n, -lambda, z, 1, grad, 1);
		f_x = 0.5*cblas_ddot(n, mWx, 1, x, 1) - lambda*cblas_ddot(n, z, 1, x, 1);
		//if(abs(f_x - func[iteration] ) < epsilon){
		//	break;
		//}
	}
	*m_maxiter = iteration+1;
	free(mWx);
	free(grad);
	free(y);
	free(x_cand);
	free(diffx);
}

// 10000 nodes two layer, maxiter=100, 10x faster than R, 211.511s vs 15.972s
void miGPFixSSTwolayer(double *W1, double *z1, double *x0, int *m_n1, double *x, 
					double *W2, double *z2, double *y0, int *m_n2, double *y,
					double *A, double *func, double *m_a1, double *m_a2, 
					double *m_lambda1,double *m_lambda2, double *m_lambda3, int *m_maxiter){
	int n1 = *m_n1;
	int n2 = *m_n2;
	double a1 = *m_a1;
	double a2 = *m_a2;
	double lambda1 = *m_lambda1;
	double lambda2 = *m_lambda2;
	double lambda3 = *m_lambda3;
	int maxiter = *m_maxiter;
	double t = 1;
	
	double *mWx = Malloc(double, n1);		//-W1*x
	cblas_dgemv(CblasColMajor, CblasNoTrans, n1, n1, -1, W1, n1, x, 1, 0, mWx, 1);

	double *mWy = Malloc(double, n2);		//-W2*y
	cblas_dgemv(CblasColMajor, CblasNoTrans, n2, n2, -1, W2, n2, y, 1, 0, mWy, 1);

	double *Ay = Malloc(double, n1);		//A*y
	cblas_dgemv(CblasColMajor, CblasNoTrans, n1, n2, 1, A, n1, y, 1, 0, Ay, 1);

	double *tAx = Malloc(double, n2);		//t(A)*y
	
	double *grad_x = Malloc(double, n1);
	memcpy(grad_x,mWx,n1*sizeof(double));
	cblas_daxpy(n1, -lambda1, z1, 1, grad_x, 1);
	cblas_daxpy(n1, -lambda3, Ay, 1, grad_x, 1);

	double *grad_y = Malloc(double, n2);
	//memcpy(grad_y,mWy,n2*sizeof(double));
	//cblas_daxpy(n2, -lambda2, z2, 1, grad_y, 1);
	//cblas_daxpy(n2, -1, tAx, 1, grad_y, 1);

	double f_x = 0.5*cblas_ddot(n1, mWx, 1, x, 1) - lambda1*cblas_ddot(n1, z1, 1, x, 1) + 
				 0.5*cblas_ddot(n2, mWy, 1, y, 1) - lambda2*cblas_ddot(n2, z2, 1, y, 1) - 
				 lambda3*cblas_ddot(n1, x, 1, Ay, 1);
	
	//double *func = Malloc(double,maxiter);
	//memset(func,0,maxiter*sizeof(double));
	double *pix = Malloc(double,n1);
	double *x_cand = Malloc(double,n1);
	double *diffx = Malloc(double,n1);

	double *piy = Malloc(double,n2);
	double *y_cand = Malloc(double,n2);
	double *diffy = Malloc(double,n2);

	int iteration = 0;
	
	for (iteration = 0; iteration < maxiter; iteration++){
		func[iteration] = f_x;
		//printf("%f,%d\n", f_x,iteration);

		memcpy(pix, x, n1*sizeof(double));
		cblas_daxpy(n1, -1, grad_x, 1, pix, 1);
		epENNORM(pix,&t,&n1,x_cand,&a1);
		
		cblas_dgemv(CblasColMajor, CblasTrans, n1, n2, 1, A, n1, x_cand, 1, 0, tAx, 1);
		memcpy(grad_y,mWy,n2*sizeof(double));
		cblas_daxpy(n2, -lambda2, z2, 1, grad_y, 1);
		cblas_daxpy(n2, -lambda3, tAx, 1, grad_y, 1);

		memcpy(piy, y, n2*sizeof(double));
		cblas_daxpy(n2, -1, grad_y, 1, piy, 1);
		epENNORM(piy,&t,&n2,y_cand,&a2);

		memcpy(diffx, x_cand, n1*sizeof(double));
		cblas_daxpy(n1, -1, x, 1, diffx, 1);

		memcpy(diffy, y_cand, n2*sizeof(double));
		cblas_daxpy(n2, -1, y, 1, diffy, 1);
		
		//printf("diff: %f,%f\n", sqrt(cblas_ddot(n1, diffx, 1, diffx, 1)), 
		//	sqrt(cblas_ddot(n2, diffy, 1, diffy, 1)));
		
		if(sqrt(cblas_ddot(n1, diffx, 1, diffx, 1)) < epsilon && 
			sqrt(cblas_ddot(n2, diffy, 1, diffy, 1)) < epsilon){
			break;
		}

 		memcpy(x, x_cand, n1*sizeof(double));
 		memcpy(y, y_cand, n2*sizeof(double));

 		cblas_dgemv(CblasColMajor, CblasNoTrans, n2, n2, -1, W2, n2, y, 1, 0, mWy, 1);
		
 		//cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, -1, W, n, x, 1, 0, mWx, 1);
		//memcpy(grad,mWx,n*sizeof(double));
		//cblas_daxpy(n, -lambda, z, 1, grad, 1);

		cblas_dgemv(CblasColMajor, CblasNoTrans, n1, n1, -1, W1, n1, x, 1, 0, mWx, 1);
		cblas_dgemv(CblasColMajor, CblasNoTrans, n1, n2, 1, A, n1, y, 1, 0, Ay, 1);
		memcpy(grad_x,mWx,n1*sizeof(double));
		cblas_daxpy(n1, -lambda1, z1, 1, grad_x, 1);
		cblas_daxpy(n1, -lambda3, Ay, 1, grad_x, 1);

		cblas_dgemv(CblasColMajor, CblasNoTrans, n2, n2, -1, W2, n2, y, 1, 0, mWy, 1);
		

		f_x = 0.5*cblas_ddot(n1, mWx, 1, x, 1) - lambda1*cblas_ddot(n1, z1, 1, x, 1) + 
				 0.5*cblas_ddot(n2, mWy, 1, y, 1) - lambda2*cblas_ddot(n2, z2, 1, y, 1) - 
				 lambda3*cblas_ddot(n1, x, 1, Ay, 1);

		//if(abs(f_x - func[iteration] ) < epsilon){
		//	break;
		//}
	}
	*m_maxiter = iteration+1;
	free(mWx);
	free(mWy);
	free(Ay);
	free(tAx);
	free(grad_x);
	free(grad_y);
	free(pix);
	free(piy);
	free(x_cand);
	free(y_cand);
	free(diffx);
	free(diffy);
}

// when L layer share the same set of nodes, each layer has n nodes but with different z
// 10000 nodes with 5 layers, 10x faster than R, maxiter=100, 36.063s vs 2.527s
void miGPFixSSMultilayer(double *W, double *listz, int *m_L, double *x0, int *m_n, double *x, 
			double *func, double *m_a, double *m_lambda, int *m_maxiter){

	int n = *m_n;
	int L = *m_L;
	double a = *m_a;
	double lambda = *m_lambda;
	int maxiter = *m_maxiter;
	double t = 1;

	double *mWx = Malloc(double, n);		//-W*x
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, -1, W, n, x, 1, 0, mWx, 1);
	double *grad = Malloc(double, n);
	memcpy(grad,mWx,n*sizeof(double));


	//cblas_daxpy(n, -lambda, z, 1, grad, 1);
	double f_x = 0.5*cblas_ddot(n, mWx, 1, x, 1);
	
	double *z = Malloc(double, n);
	for (int i = 0; i < L; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			z[j] = listz[i*n+j];
		}
		cblas_daxpy(n, -lambda, z, 1, grad, 1);
		f_x -= lambda*cblas_ddot(n, z, 1, x, 1);
	}

	//double *func = Malloc(double,maxiter);
	//memset(func,0,maxiter*sizeof(double));
	double *y = Malloc(double,n);
	double *x_cand = Malloc(double,n);
	double *diffx = Malloc(double,n);
	int iteration = 0;
	for (iteration = 0; iteration < maxiter; iteration++){
		func[iteration] = f_x;
		
		memcpy(y, x, n*sizeof(double));
		cblas_daxpy(n, -1, grad, 1, y, 1);
		
		//epEN(y, 1, n, x_cand, a);
		epENNORM(y,&t,&n,x_cand,&a);
		memcpy(diffx, x_cand, n*sizeof(double));
		cblas_daxpy(n, -1, x, 1, diffx, 1);
		
		if(sqrt(cblas_ddot(n, diffx, 1, diffx, 1)) < epsilon){
			break;
		}

 		memcpy(x, x_cand, n*sizeof(double));
 		cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, -1, W, n, x, 1, 0, mWx, 1);
		memcpy(grad,mWx,n*sizeof(double));

		f_x = 0.5*cblas_ddot(n, mWx, 1, x, 1);
		for (int i = 0; i < L; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				z[j] = listz[i*n+j];
			}
			cblas_daxpy(n, -lambda, z, 1, grad, 1);
			f_x -= lambda*cblas_ddot(n, z, 1, x, 1);
		}

		//cblas_daxpy(n, -lambda, z, 1, grad, 1);
		//f_x = 0.5*cblas_ddot(n, mWx, 1, x, 1) - lambda*cblas_ddot(n, z, 1, x, 1);

	}
	*m_maxiter = iteration+1;
	free(mWx);
	free(grad);
	free(y);
	free(z);
	free(x_cand);
	free(diffx);

}
/*
int main(int argc, char const *argv[])
{
	double W[9]={1,2,0,4,5,0,0,0,3};
	double z[3]={1,2,-3};
	double x0[3]={0.3,0.2,0.3};
	double x[3]={0.3,0.2,0.3};
	double func[5]={0,0,0,0,0};
	int n = 3;
	double alpha = 0.5;
	double lambda = 1;
	int maxi = 10;

	//miGPFixSS(W, z, x0, &n, x, func, &alpha, &lambda, &maxi);

	//OutPutVector(x, 3);
	//OutPutVector(func, 3);
	return 0;
}
*/