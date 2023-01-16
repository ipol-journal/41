/*
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

/**
 * @mainpage Selective Contrast Adjustment
 *
 * README.txt:
 * @verbinclude README.txt
 */

/**
 * @file poisson_lca.c
 * @brief command-line interface
* 
* Given an image, the algorithm computes the dark regions 
* of the gray intensity levels (by threshold), then it computes 
* the new gradient vector field of the gray amplifying 
* the gradient vector in the dark regions and finally
* it recovers the new image by solving the poisson 
* equation using the Fast Fourier Transform. 
* 
* For obtainig the resulting image, the solution
* is normalized to [0,255] and computed the color proportionally.
*
* The input is also normalized by the same normalization process.
*
* @author Catalina Sbert <catalina.sbert@uib.es/>
* @author Ana Belén Petro <anabelen.petro@uib.es/>
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "poisson_lca_lib.h"
#include "auxiliary_lib.h"
#include "io_png.h"

#include "parser.h"

using namespace std;


/**
 * @brief main function call
 */
int main(int argc, char **argv)
{
        double T;  /* threshold dark region */
        double a;  /* amplification factor for the gradient in dark zones*/
        double alpha; /* power factor for the gardient */ 
        size_t nx, ny, nc, dim;
        unsigned char *data_in, *output;
        double  *data, *data_norm, *out_local, *out_global;
        int channel;
        float s;
        
        
        std::vector <OptStruct *> options;
        OptStruct ot  = {"t:", 0, "50", NULL, "thresold dark regions"}; 
	    options.push_back(&ot); 
        OptStruct oA  = {"a:", 0, "3", NULL, "amplification factor for the gradient"}; 
		options.push_back(&oA);      
        OptStruct oB  = {"b:", 0, "0.8",NULL, "Power of the gradient"}; 
	    options.push_back(&oB);     
        
        
        
        std::vector<ParStruct *> pparameters;
        ParStruct pinput = {"input", NULL, "input file"}; 
	    pparameters.push_back(&pinput);
        ParStruct psim_gray = {"sim_gray", NULL, 
			"simplest color balance on the gray "}; 
	    pparameters.push_back(&psim_gray);
        ParStruct pout_local_gray={"loc_gray", NULL, 
			"result on the gray modifying gradient on dark region"}; 
	    pparameters.push_back(&pout_local_gray);
        ParStruct pout_global_gray={"glo_gray", NULL, 
			"result on the gray modifying gradient by a power function"}; 
	    pparameters.push_back(&pout_global_gray);
        ParStruct psim_rgb = {"sim_rgb", NULL,
			"simplest color balance on rgb channels"}; 
	    pparameters.push_back(&psim_rgb);
        ParStruct pout_local_rgb={"loc_rgb", NULL, 
			"result on rgb modifying gradient on dark region"};
	    pparameters.push_back(&pout_local_rgb);
        ParStruct pout_global_rgb={"glo_rgb", NULL, 
			"result on rgb modifying gradient by a power function"};
	    pparameters.push_back(&pout_global_rgb);
        
        
        if (!parsecmdline("poisson_lca",
						  "Modify the gradient and solve the Poisson equation", 
						  argc, argv, options, pparameters))
                return EXIT_FAILURE;
        
        T=atof(ot.value);
     
        
    /* threshold for the dark regions */

   /* T = atof(argv[1]);*/
    if (0. > T || 255. < T)
    {
        fprintf(stderr, "the  threshold must be in [0..255]\n");
        return EXIT_FAILURE;
    }
     
    /* amplification factor for the gradient*/
    a=atof(oA.value);
    if(0.>= a || 10. < a)
    {
        fprintf(stderr, "the  factor must be in (0..10]\n");
        return EXIT_FAILURE;
    }
    /* power of the norm of gradient*/
    alpha=atof(oB.value);
    if(0.>=alpha || 3. < alpha)
    {
        fprintf(stderr, "the  alpha must be in (0..3]\n");
        return EXIT_FAILURE;
    }
    

    /* read the PNG image into data */   
        
        if (NULL == (data_in = io_png_read_u8(pinput.value, &nx, &ny, &nc)))
    {
        fprintf(stderr, "the image could not be properly read\n");
        return EXIT_FAILURE;
    }
    dim=nx*ny;
        
        s=0.5;

    /* the image has either 1 or 3 non-alpha channels */
    if (3 <= nc)
        channel = 3;
    else
        channel = 1;
        

    /* allocate data and  data_out  */
        
    if (NULL == (out_local = (double *) malloc(nc * nx * ny * sizeof(double))) ||
        NULL == (data_norm = (double *) malloc(nc * nx * ny * sizeof(double))) ||
                NULL == (data = (double *) malloc(nc * nx * ny * sizeof(double))) ||
                NULL == (output = (unsigned char *) malloc(nc * nx * ny * sizeof(unsigned char)))||
                NULL == (out_global = (double *) malloc(nc * nx * ny * sizeof(double))))
    {
        fprintf(stderr, "allocation error\n");
        return EXIT_FAILURE;
    }

     
        
        /* transform input into a double array*/
        
        int_double(data, data_in, dim, channel);
        
        
        
        if(channel==1){ /* the input image is a gray image*/
                                
                poisson_gray(data_norm, out_local, out_global, data, T, a,
							 alpha, nx,ny, channel,s);
                
                double_int(output, data_norm, dim, channel);
                io_png_write_u8(psim_gray.value, output, nx, ny, nc);
                
                double_int(output, out_local, dim, channel);
                io_png_write_u8(pout_local_gray.value, output, nx, ny, nc);

                double_int(output, out_global, dim, channel);
                io_png_write_u8(pout_global_gray.value, output, nx, ny, nc);
                
        }
                
        else{ /* the input image is a color image*/
                
                /*------- The algorithm on the gray intensity channel -----*/
                
                
                poisson_gray(data_norm, out_local, out_global, 
							 data, T, a, alpha, nx,ny, channel, s);
                
                double_int(output, data_norm, dim, channel);
                io_png_write_u8(psim_gray.value, output, nx, ny, nc);
                
                double_int(output, out_local, dim, channel);
                io_png_write_u8(pout_local_gray.value, output, nx, ny, nc);
                
                double_int(output, out_global, dim, channel);
                io_png_write_u8(pout_global_gray.value, output, nx, ny, nc);
                
                
                
                /*---------The algorithm on RGB channel-----------*/
                        
                poisson_rgb(data_norm, out_local, out_global, data, T, a, 
							alpha, nx,ny, channel,s);
                
                double_int(output, data_norm, dim, channel);
                io_png_write_u8(psim_rgb.value, output, nx, ny, nc);
                
                double_int(output, out_local, dim, channel);
                io_png_write_u8(pout_local_rgb.value, output, nx, ny, nc);
                
                double_int(output, out_global, dim, channel);
                io_png_write_u8(pout_global_rgb.value, output, nx, ny, nc);
                
        }

    free(out_local); free(out_global); free(output); 
	free(data_norm); free(data);
        
    return EXIT_SUCCESS;
}

