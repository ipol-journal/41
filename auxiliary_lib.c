/*
 *  auxiliary_lib.c
 *  
 *
 * Created by Catalina Sbert Juan on 13/11/12.
 *
 * Copyright 2012 IPOL Image Processing On Line http://www.ipol.im/
 *
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



/**
 *@brief Transform an unsigned char array into a double array
 *
 * @param data_in input unsigned char array
 * @return data double array
 */

double *int_double(double *data, unsigned char *data_in, size_t dim, int channel)
{
	int l,n;
	
	for(l=0;l< (int) dim; l++){
		for(n=0;n<channel;n++)
			data[l+n*(int)dim]=(double)data_in[l+n*(int)dim];
	}
	
	return data;
}


/**
 *@brief Transform a double array into an unsigned char array
 *
 * @param data_in input unsigned char array
 * @return data double array
 */

unsigned char *double_int(unsigned char *data, double *data_in, size_t dim, int channel)
{
	int l,n;
	
	for(l=0;l< (int) dim; l++){
		for(n=0;n<channel;n++)
			data[l+n*(int)dim]=(unsigned char)data_in[l+n*(int)dim];
	}
	
	return data;
}
/**
 * @brief Computes the gray intensity value of a color image
*
* @f$ gray= (R+ G+ B)/3 \f$
*
* @param data_in input color image
* @param dim size of the image
* 
* @return gray output gray image
*
*/

double *gray_intensity(double *gray, double *data_in, size_t dim)
{
    int l;
	
	for(l=0; l< (int) dim; l++) 
		gray[l]=(data_in[l]+data_in[l+dim]+data_in[l+2*dim])/3.;
	
	return gray;
	
	
}

 
 
 
 
 
/**
 *  @brief compute the R G B components of the output image from its gray level 
 *
 *   Given a color image C=(R, G, B), given its gray level
 *
 * @f$ gray= (R+ G+ B)/3 \f$
 *
 * Given a modified gray image gray1
 *
 * This function computes an output color image C1=(R1,G1,B1) which each channel is proportinal 
 * to the input channel and whose gray level is gray1, 
 *
 * @f$ R1=\frac{gray1}{gray} R    G1=\frac{gray1}{gray} G    B1= \frac{gray1}{gray} B \f$
 *
 *Note that we make a restriction and does not permit a factor  greater than 3
 *
 * @param data_out Output color image
 * @param data input color image
 * @param gray gray level of the input color image
 * @param gray1 modified gray image
 * @param dim size of the image
 *
 * @return data_out
 */

double *color(double *data_out, double *data, double *gray, double *gray1, size_t dim)
{
	
	double *ptr_red, *ptr_green, *ptr_blue;
	double *ptr_end,*ptr_gray, *ptr_gray1;
	double *ptr_in_red, *ptr_in_green, *ptr_in_blue;
	double A,B;
	
	/* sanity check*/
	if (NULL == data_out || NULL == data || NULL== gray || NULL== gray1)
	{
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
	}
	
	
	
	ptr_red=data_out;
	ptr_green=data_out+dim;
	ptr_blue=data_out+2*dim;
	ptr_gray=gray;
	ptr_gray1=gray1;
	ptr_end=ptr_gray+dim;
	ptr_in_red=data; ptr_in_green=data+dim; ptr_in_blue=data+2*dim;
	while(ptr_gray< ptr_end){
		if(*ptr_gray <= 1.) *ptr_gray=1.;
		A=*ptr_gray1/ *ptr_gray;
		if(A > 3.) A=3.;
		*ptr_red=A* (*ptr_in_red);
		*ptr_green=A*(*ptr_in_green);
		*ptr_blue=A*(*ptr_in_blue);
		if( *ptr_red > 255. || *ptr_green > 255. || *ptr_blue> 255.){
			B=*ptr_in_red;
			if(*ptr_in_green > B) B=*ptr_in_green;
			if( *ptr_in_blue > B) B=*ptr_in_blue;
			A= 255. /B;
			*ptr_red=A* (*ptr_in_red);
			*ptr_green=A*(*ptr_in_green);
			*ptr_blue=A*(*ptr_in_blue);
		}
		ptr_gray++; ptr_gray1++; ptr_red++; ptr_green++; ptr_blue++;
		ptr_in_red++;  ptr_in_green++; ptr_in_blue++;
	}
	return data_out;
}
int myComparisonFunction(const void *x, const void *y) {
	
    /* x and y are pointers to doubles.
	 Returns  -1 if x < y
	           0 if x == y
			  +1 if x > y*/
	
    double dx, dy;
	
    dx = *(double *)x;
    dy = *(double *)y;
	
    if (dx < dy) {
        return -1;
    } else if (dx > dy) {
        return +1;
    }
    return 0;
}

double *simplest_color_balance(double *data_out, double *data, size_t dim, float s)
{
	
	
	int per;
	int l;
	double min,max, scale;
	double *sortdata;
	
	sortdata= (double*) malloc(dim*sizeof(double));
	
	memcpy(sortdata, data, dim*sizeof(double));
	
	qsort(sortdata, dim, sizeof sortdata[0], &myComparisonFunction);
	
	per = (int) (s*dim/100);
	min=sortdata[per];
	max= sortdata[dim-1-per];
	
	
	
	if(max  <= min)  for(l=0;l< (int) dim; l++) data_out[l]=max;
	else{
		scale=255./(max-min);
		for(l=0;l< (int)dim; l++){
			if(data[l] < min) data_out[l]=0.;
			else if(data[l]> max) data_out[l]=255.;
			else data_out[l]=scale*(data[l]-min);
		}
	}
	
	return data_out;
	
	
	
}

