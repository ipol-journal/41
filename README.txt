% Selective Contrast adjustment 

# ABOUT

* Author    : Catalina Sbert  <catalina.sbert@uib.es>
              Ana Belén Petro <anabelen.petro@uib.es>
* Copyright : (C) 2009-2012 IPOL Image Processing On Line http://www.ipol.im/
* License   : GPL v3+, see GPLv3.txt

# OVERVIEW

Given an image, the algorithm modifies the gradient vector field of the image  and finally it recovers the new image by solving the poisson equation using the Fast Fourier Transform as described in IPOL
     http://www.ipol.im/pub/algo/lmps_selective_contrast_adjustment/
The modification of the gradient is made in two ways:
1) Amplifying the gradient vector in the dark regions of the image. These dark regions are computed by threshold.
2) Applying a power function to the gradient in all the pixels of the image.

This program reads a PNG image, given T a threshold, and an amplification factor a and a power alpha, computes the new image. 
The program works on each color channel independently and on the gray intensity channel and then computes the color such that it preserves the hue.
# REQUIREMENTS

The code is written in ANSI C, and should compile on any system with
an ANSI C compiler.

The libpng header and libraries are required on the system for
compilation and execution. See http://www.libpng.org/pub/png/libpng.html

The fftw3 header and libraries are required on the system for
compilation and execution. See http://www.fftw.org/

# COMPILATION

Simply use the provided makefile, with the command `make`.

# USAGE

poisson_lca [-t t] [-a a] [-b b]  input  sim_gray  loc_gray  glo_gray  sim_rgb  loc_rgb  glo_rgb 
	-t  t	 thresold dark regions (Default: 50)
	-a  a	 amplification factor for the gradient (Default: 3)
	-b  b	 Power of the gradient (Default: 0.8)
	input	 input file
	sim_gray	 simplest color balance on the gray intensity
	loc_gray	 result on the gray modifying gradient on dark region
	glo_gray	 result on the gray modifying gradient by a power function
	sim_rgb	 simplest color balance on rgb channels
	loc_rgb	 result on rgb modifying gradient on dark region
	glo_rgb	 result on rgb modifying gradient by a power function




# ABOUT THIS FILE

Copyright 2009-2012 IPOL Image Processing On Line http://www.ipol.im/
Author: Catalina Sbert <catalina.sbert@uib.es>

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
