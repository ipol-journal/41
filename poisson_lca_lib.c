/*
 * Copyright 2009, 2010 IPOL Image Processing On Line http://www.ipol.im/
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file poisson_lca_lib.c
 * @brief   gradient, divergence  and Poisson routines using DFT. 
 *
 * @author Catalina Sbert <catalina.sbert@uib.es/>
 * @author Ana Belén Petro <anabelen.petro@uib.es/>
 */


#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#include "auxiliary_lib.h"

/*#include "poisson_lca_lib.h"*/

/* M_PI is a POSIX definition */
#ifndef M_PI
/** macro definition for Pi */
#define M_PI 3.14159265358979323846
#endif                          /* !M_PI */


	

/**
 * @brief compute the gradient of a 2D array using the discrete cosinus transform 
 *
 * We supose that the input image (N X M array) is extended simmetrically across its sides to a 
*  2N X 2M array. Then we obtain a function wich is even around (-0.5 , j) ,(N-0.5,j ), 
 * (i,-0.5) and (i, M-0.5). To compute the gradient we use the following definition 
*  of cosinus transform, (we write it in dimension 1 for dimension two is similar)
*    
*    @f$       C(X)_k= 2\sum_{j=0}^{N-1} X_j cos( \frac{\pi k (j+1/2)}{N}) \f$
* Note that with the symmetry C(X)_N=0.
* 
* With this definition the sinus transform of the derivative is defined by 
*  
* @f$       S(X')_k=\frac{-\pi(k+1)}{N } C(X)_{k+1} , k=0,..., N-1  \f$
* 
*  This new function is odd around -1 and even aroun N-1
* Then we compute X'_k using the invers sinus  transform whose definition is 
*
* @f$  X'_k=\frac{2\sum_{j=0}^{N-1} S(X')_j sin (\frac{\pi(j+1)(k+1/2)}{N})}{2N} \f$
*
* Note that the result is odd around -0.5 and around N-0.5.
*
* @param data_x  partial derivative respect to the first variable
* @param data_y  partial derivative respect to the second variable
* @param data  input array
* @param nx x-the dimensions of the array
* @param ny y-the dimensions of the array
*
* @return data_x, data_y  or NULL if an error occured
*/ 


static double *gradient(double *data_x, double *data_y, double *data_in, size_t nx, size_t ny)
{
	
	int x,y,dim,dim4;
	double *ftout, *ftux, *ftuy;
	fftw_plan p;
	double normx, normy;
	/*  pointers to the right and bottom neighbour values */
	double *ptr_end, *ptr_in_x1,*ptr_in_y1, *ptr_x, *ptr_y;
	
	/* sanity check */
    if (NULL == data_in || NULL == data_x || NULL== data_y)
    {
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
    }
	
	dim=(int)nx*(int)ny;
	dim4=4*dim;
	
	/* allocate the memory of  fourier transform*/
	ftout=(double*)fftw_malloc(sizeof(double)*dim);
	ftux=(double*)fftw_malloc(sizeof(double)*dim);
	ftuy=(double*)fftw_malloc(sizeof(double)*dim);
	
	/* compute the cosinus transform of the input array*/
	
	p=fftw_plan_r2r_2d((int) ny,(int) nx, data_in, ftout,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
	fftw_execute(p); 
	
	
	
	normx=M_PI/nx; normy=M_PI/ny;
	
	/* compute the Fourier transform of partial derivatives respect to x and y*/
	
	ptr_in_x1=ftout+1;
	ptr_in_y1=ftout+nx;
	ptr_x=ftux;
	ptr_y=ftuy; 
	
	for(y=0;y<(int)ny;y++)
		for(x=0;x< (int) nx;x++)
		{
			if(x < (int) nx-1) *ptr_x=-normx*(x+1)*(*ptr_in_x1);
			else *ptr_x=0.;
			if(y < (int) ny-1) *ptr_y=-normy*(y+1)*(*ptr_in_y1);
			else *ptr_y=0.;
			ptr_x++;
			ptr_y++;
			ptr_in_x1++;
			ptr_in_y1++;
		}
	
	/* compute the inverse Fourier transform to obtain the gradient of the input array*/
	
	p=fftw_plan_r2r_2d((int) ny,(int) nx, ftux, data_x, FFTW_REDFT01,FFTW_RODFT01,FFTW_ESTIMATE);
	fftw_execute(p);
	
	p=fftw_plan_r2r_2d((int) ny, (int) nx, ftuy, data_y, FFTW_RODFT01,FFTW_REDFT01,FFTW_ESTIMATE);
	fftw_execute(p);
	
	/* normalize the inverse fourier transform*/         
	ptr_end=data_x+dim;
	ptr_x=data_x;
	ptr_y=data_y;
	while(ptr_x < ptr_end){
		*ptr_x=*ptr_x/dim4;
		*ptr_y=*ptr_y/dim4;
		ptr_x++; ptr_y++;
	}
	/*deallocate the memory*/  
	fftw_destroy_plan(p);
	fftw_free(ftout);
	fftw_free(ftux); fftw_free(ftuy);
	
	return data_x; return data_y;
	
}

/**
 * @brief compute the Fourier transform of the divergence of a vector field using the 
*  @brief       discrete sinus/cosinus transform
* 
*  We have two  N X M arrays V_1 and V_2. V_1 is  odd around (-0.5, j) and (N-0.5, j) and even around 
* (i,-0.5) and (i,M-0.5). V_2 is even  around (-0.5, j) and (N-0.5, j) and odd around 
* (i,-0.5) and (i,M-0.5).
* We compute the Fourier transform of  divergence of the vector field V= (V_1, V_2), 
*
*     @f$  div(V)=\frac{\partial V_1}{\partial x}+\frac{\partial V_2}{\partial y}  \f$
*
*   Due to the kind of symmetry to compute the Fourier transform of
* @f$ \frac{\partial V_1}{\partial x}\f$ ( respectively @f$ frac{\partial V_2}{\partial y}\f$)
* we use the sinus transform, which in dimension 1 is defined by
*
* @f$ S(X)_k=2\sum_{j=0}^{N-1} X_j sin(\frac{\pi(j+1/2)(k+1)}{N}) \f$ 
*
* The cosinus transform of the derivative of X is defined by
*
* @f$ C(X')_k=\frac{\pi k}{N} S(X)_{k-1} ,  k=1,... N-1,  C(X')_0=0 \f$ 
* 
* @param div  the cosinus transform of the divergence of a vector field
* @param data_x the first component of the vector field
* @param data_y the second component of the vector field
* @param nx x-the dimensions of the array
* @param ny y-the dimensions of the array 
*
* @return div or NULL 
*/

static double *divergence(double *div, double *data_x, double *data_y, size_t nx, size_t ny)
{
	
	
	int x,y, dim;
	double *ftout,*ptr_div,*ptr_x1;
	fftw_plan p;
	double normx, normy;
	
    /* sanity check*/
	if (NULL == div || NULL == data_x || NULL== data_y)
    {
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
    }
	
	
	dim=(int)nx*(int)ny;
	normx=M_PI/nx; normy=M_PI/ny;
	
	/* allocate the Fourier transform*/
	
	ftout=(double*)fftw_malloc(sizeof(double)*dim);
	
	/* compute the Fourier transform of the first component of the vector field*/
	p=fftw_plan_r2r_2d((int) ny,(int) nx, data_x, ftout,FFTW_REDFT10,FFTW_RODFT10,FFTW_ESTIMATE);
	fftw_execute(p); 
    
    /*compute the fourier transform of the derivative respect to x of the first component*/    
	ptr_x1=ftout;
	ptr_div=div;
	for(y=0;y< (int) ny; y++)
		for(x=0; x< (int) nx; x++){
			if(x==0) *ptr_div=0.;
			else {
				*ptr_div= normx*x* *(ptr_x1-1);
			}
			ptr_div++;
			ptr_x1++;
		}
	
    /* compute the Fourier transform of the second component of the vector field*/
	p=fftw_plan_r2r_2d((int) ny, (int) nx,data_y, ftout,FFTW_RODFT10,FFTW_REDFT10,FFTW_ESTIMATE);
	fftw_execute(p); 
	
	/* compute the fourier transform of the derivative respect to y of the second component and the divergence*/   
	ptr_x1=ftout;
	ptr_div=div+nx;
	for(y=1;y< (int) ny; y++)
		for(x=0; x< (int) nx; x++){
			*ptr_div+=normy*y*(*ptr_x1);
			ptr_div++;
			ptr_x1++;
		}
	
    /* deallocate the memory*/
	fftw_destroy_plan(p);
	fftw_free(ftout);
	
	return div;
	
}
/**
* @brief Solve a Poisson equation using FFTW
*
* Given a vector filed V=(ux,uy), this routine solves the Poisson equation
* @f$ \Delta u= div(V)\f$ with homogeneous Neumann boundary conditions
*
*
* @li compute  the  Fourier transform of the divergence of the vector field V, C(di.
* @li compute the cosinus transform of the laplacian
*  @f$ C(\Delta u)_{ij}=-[(\frac{\pi i}{nx})²+(\frac{\pi j}{ny})²] C(u)_{ij}\f$
* @li The cosinus transform of the solution of the Poisson equation is
* @f$ C(u)_{ij}=-\frac{C(div)}{(\frac{\pi i}{nx})²+(\frac{\pi j}{ny})²}
* @li the solution is computed by the inverse cosinus transform 
*
* @param ux first component of the vector
* @param uy second component of the vector
* @param nx x-dimension of the array
* @param ny y-dimension of the array
* 
* @return data_out
 
 
*/
double *solve_poisson_equation(double *data_out, double *ux, double *uy, size_t nx, size_t ny)
{
	
	
	double *div;
	double normx, normy, A;
	int x,y,l, dim, dim4;
	fftw_plan p;
	
	
	dim=(int)nx*(int)ny;
	dim4=4*dim;
	div=(double*)fftw_malloc(sizeof(double)*dim);

	
	/* compute the fourier transform of the divergence of the vector (ux,uy)*/
	
	(void) divergence(div, ux,uy, nx,ny);
	
	normx=M_PI/nx; normy=M_PI/ny;
	normx*=normx; normy*=normy;
	
	/* compute the fourier transform of the solution of the Poisson equation*/     
	
	for(y=0;y< (int) ny;y++){
		l=y*(int)nx;
		for(x=0;x< (int) nx;x++){
			if( x==0 && y==0) div[0]=0.;
			else{
                A=(normx*x*x+normy*y*y);
                div[l+x]=-div[l+x]/A;
			}
		}
	}
	
    /* compute the inverse fourier transform*/
	p=fftw_plan_r2r_2d((int) ny, (int) nx, div, data_out,FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
	fftw_execute(p); 
	
	for(l=0; l< dim; l++) data_out[l]/=dim4;
	
	return  data_out;
	
}

/**
* @brief  Linear option for modifying the gradient, solving a Poisson equation using the Fourier Transform
*
*  The input array is processed as follow:
*
* @li compute the gradient of the input array
* @li Modify the gradient: amplifying by a factor the gradient in the interior of the dark zones of the image
* @li solve the Poisson equation aplying FFT
*
* @param data_in input array
* @param T threshold for the dark regions
* @param a amplification factor for the gradient
* @param nx x-dimension of the array
* @param ny y-dimension of the array
* 
* @return data_out
*/
 double  *poisson_local(double *data_out, double *data_in, double T, double a, size_t nx, size_t ny)
{


	double *ux,*uy;
	double *ptr_x,*ptr_y, *ptr_data,  *ptr_xr, *ptr_xl,*ptr_yt,*ptr_yb;
	int x,y, dim;
	

	if (NULL == data_out || NULL == data_in )
	{
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
    }

	dim=(int)nx*(int)ny;
	
	/* allocate the memory for the Fourier transform*/
      
	ux=(double*)fftw_malloc(sizeof(double)*dim);
	uy=(double*)fftw_malloc(sizeof(double)*dim);
	

	/* compute the gradient of the input array*/

	gradient(ux, uy, data_in, nx,ny);
          
	/* multiply the gradient by a factor a in the interior of the dark zones which are defined by a threshold T*/
   
	ptr_x=ux; ptr_y=uy;
	ptr_data=data_in;
	ptr_xr=data_in+1;
	ptr_xl=data_in-1;
	ptr_yt=data_in-nx;
	ptr_yb=data_in+nx;
	for(y=0; y < (int) ny; y++)
		for(x=0; x < (int) nx; x++){
			if( x >0  && x < (int) nx-1  && y>0 && y < (int) ny-1)
			if(*ptr_xl <= T && *ptr_data <= T && *ptr_xr <= T && *ptr_yt <= T && *ptr_yb <= T){
				*ptr_x=a*(*ptr_x); *ptr_y=a* (*ptr_y);
			}
			ptr_data++; ptr_x++; ptr_y++; ptr_xr++; ptr_xl++; ptr_yt++; ptr_yb++;
		}

	  
	
	solve_poisson_equation(data_out,ux,uy, nx,ny);
	
	fftw_free(ux); fftw_free(uy);
	
    return data_out;
}


/**
* @brief  Power option for modifying the gradient, solving a Poisson equation using the Fourier Transform
*
*  The input array is processed as follow:
*
* @li compute the gradient of the input array
* @li Modify the gradient: power by alpha the gradient
* @li solve the Poisson equation
*
* @param data_in input array
* @param alpha Power function for the gradient
* @param nx x-dimension of the array
* @param ny y-dimension of the array
*
* @return data_out
*/

double  *poisson_global(double *data_out, double *data_in, double a, size_t nx, size_t ny)
{
	
	double *ux,*uy;
	int l,dim;
	double A;
	
	
	if (NULL == data_out || NULL == data_in )
    {
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
    }
	
	dim=(int)nx*(int)ny;
	
	/* allocate the memory for the Fourier transform*/
	
	ux=(double*)fftw_malloc(sizeof(double)*dim);
	uy=(double*)fftw_malloc(sizeof(double)*dim);
	
	
	/* compute and modify the gradient of the input array*/
	
	gradient(ux, uy, data_in, nx,ny);
	
	for(l=0;l< dim;l++){
		A=fabs(ux[l])+fabs(uy[l])+0.1;
		A=pow(A, a-1.);
		ux[l]*=A; 
		uy[l]*=A;
	}
	
	
	solve_poisson_equation(data_out, ux, uy, nx,ny);
	
	fftw_free(ux); fftw_free(uy);
	
	return data_out;
}

/**
 * @brief  The algorithm applied to the gray intensity and the color channels are  computed 
 *  such that the ratios R/G/B are preserved.
 *
 *
 * @param data input array
 * @param T the threshold to compute the dark regions
 * @param a the gradient amplification factor 
 * @param alpha Power function for the gradient
 * @param nx x-dimension of the array
 * @param ny y-dimension of the array
 * @param channel the number of channels of the input image
 * @param s the percentage of pixels saturated to minimum/maximum value
 *
 * @return data_norm  simplest color balance of the input with s% of saturation
 * @return local_out the result of the local option of the algorithm
 * @return global_out the result of the global option of the algoritm
 */




double *poisson_gray(double *data_norm, double *local_out, double *global_out, 
					 double *data, double T, double a, double alpha, size_t nx, size_t ny, int channel, float s)
{
	
	int dim;
	double *data_gray,*data_gray1;
	
	dim=(int)nx*(int)ny;
	
	data_gray = (double*) malloc(dim* sizeof(double));
	data_gray1 = (double*) malloc(dim * sizeof(double));
	
	
	if(channel == 1){
		
		simplest_color_balance(data_norm, data, dim, s);
		
		poisson_local(local_out, data_norm, T, a, nx,ny);
		
		simplest_color_balance(local_out, local_out, dim,s);
		
		poisson_global(global_out, data_norm, alpha, nx,ny);
		
		simplest_color_balance(global_out, global_out, dim,s);
		
	}
	else{
		
		gray_intensity(data_gray, data, dim);
		
		simplest_color_balance(data_gray1, data_gray, dim,s);
		
		color(data_norm, data, data_gray, data_gray1, dim);
		
		poisson_local(data_gray, data_gray1, T, a, nx,ny);
		
		simplest_color_balance(data_gray, data_gray, dim,s);
		
		color(local_out, data_norm, data_gray1, data_gray, dim);
		
		poisson_global(data_gray, data_gray1, alpha, nx,ny);
		
		simplest_color_balance(data_gray, data_gray, dim,s);
		
		color(global_out, data_norm, data_gray1, data_gray, dim);
		
		
	}
	

	free(data_gray); free(data_gray1);
	
	return data_norm;
	return local_out;
	return global_out;

}

/**
 * @brief  The algorithm applied to  the color channels independently 
 *
 * @param data input array
 * @param T the threshold to compute the dark regions
 * @param a the gradient amplification factor 
 * @param alpha Power function for the gradient
 * @param nx x-dimension of the array
 * @param ny y-dimension of the array
 * @param channel the number of channels of the input image
 * @param s the percentage of pixels saturated to minimum/maximum value
 *
 * @return data_norm  simplest color balance of the input with s% of saturation
 * @return local_out the result of the local option of the algorithm
 * @return global_out the result of the global option of the algoritm
 */


double *poisson_rgb(double *data_norm, double *local_out, double *global_out, 
					double *data, double T, double a, double alpha, size_t nx, size_t ny, int channel, float s)
{
	
	int dim, n;
	
	dim=(int)nx*(int)ny;
	
	for(n=0;n<channel;n++) simplest_color_balance(data_norm+n*dim, data+n*dim, dim, s);
	
	for(n=0;n<channel;n++){
		poisson_local(local_out+n*dim, data_norm+n*dim, T, a, nx,ny);
		simplest_color_balance(local_out+n*dim, local_out+n*dim, dim, s);
	}
	
	for(n=0;n<channel;n++){
		poisson_global(global_out+n*dim, data_norm+n*dim, alpha, nx,ny);
		simplest_color_balance(global_out+n*dim, global_out+n*dim, dim, s);
	}
	
	return data_norm;
	return local_out;
	return global_out;	
   

}
