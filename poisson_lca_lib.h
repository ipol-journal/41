
/*
 *  poisson-lca_lib.h
 *  
 *
 *  Created by Catalina Sbert Juan on 13/11/12.
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








#ifdef __cplusplus
extern "C" {
#endif
	
double *poisson_gray(double *data_norm, double *local_out, double *global_out, double *data, double T, 
					 double a, double alpha, size_t nx, size_t ny, int channel, float s);
double *poisson_rgb(double *data_norm, double *local_out, double *global_out, double *data, double T,
					double a, double alpha, size_t nx, size_t ny, int channel, float s);

#ifdef __cplusplus
}
#endif