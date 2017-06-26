// ---------------
//   Dongzhuo Li
//   March, 2015
//
//   Modified from the fast marching eikonal solver
//   written by D.Kroon University of Twente (Oct 2010)
// ---------------

#include "math.h"
#include "common.c"

// function prototype
double CalculateDistance2D(double *T, double Fij, int *dims, int i, int j, bool usesecond, bool usecross, bool *Frozen);
double CalculateDistance3D(double *T, double Fijk, int *dims, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen);
double second_derivative(double Txm1, double Txm2, double Txp1, double Txp2);
double averageVelocityAroundSource2D(int bz, int bx, int izsrcn, int ixsrcn, int *dims, double *F);
double CalculateTConstant2D(double aveVel, int izsrcn, int ixsrcn, int zloc, int xloc);
double averageVelocityAroundSource3D(int bz, int bx, int by, int izsrcn, int ixsrcn, int iysrcn, int *dims, double *F);
double CalculateTConstant3D(double aveVel, int izsrcn, int ixsrcn, int iysrcn, int zloc, int xloc, int yloc);

void msfm2dCpp(double *T, double *F, int izsrcn, int ixsrcn, int *dims, bool usesecond, bool usecross) {

	// box around the source - Dongzhuo
	int bz = 0;
	int bx = 0;

	/* Euclidian distance image */
	double *Y; // NOT USED CURRENTLY!!!

	/* Current distance values */
	double Tt, Ty;

	/* Matrix containing the Frozen Pixels" */
	bool *Frozen;

	/* Augmented Fast Marching (For skeletonize) */
	bool Ed = 0;

	// /* Size of input image */
	// const mwSize *dims_c;
	// mwSize dims[2];

	// /* Size of  SourcePoints array */
	// const mwSize *dims_sp_c;
	// mwSize dims_sp[2];

	/* Number of pixels in image */
	int npixels;
	int dims_sp[2] = {2, 1};

	/* Neighbour list */
	int neg_free;
	int neg_pos;
	double *neg_listv;
	double *neg_listx;
	double *neg_listy;
	double *neg_listo;

	int *listprop;
	double **listval;

	/* Neighbours 4x2 */
	int ne[8] = { -1, 1, 0, 0, 0, 0, -1, 1};

	/* Loop variables  */
	int z, k, itt, q;

	/* Current location */
	int i, j, x, y;

	/* Index */
	int IJ_index, XY_index, index;

	npixels = dims[0] * dims[1];

	// compute average velocity in the small box around the soruce - Dongzhuo
	double aveVel = averageVelocityAroundSource2D(bz, bx, izsrcn, ixsrcn, dims, F);
	// cout << "aveVel = " << aveVel << endl;

	/* Pixels which are processed and have a final distance are frozen */
	Frozen = (bool *)malloc( npixels * sizeof(bool) );
	for (q = 0; q < npixels; q++) {
		Frozen[q] = 0;
		T[q] = -1;
	}
	if (Ed) {
		for (q = 0; q < npixels; q++) {
			Y[q] = -1;
		}
	}

	/*Free memory to store neighbours of the (segmented) region */
	neg_free = 100000;
	neg_pos = 0;
	neg_listx = (double *)malloc( neg_free * sizeof(double) );
	neg_listy = (double *)malloc( neg_free * sizeof(double) );
	if (Ed) {
		neg_listo = (double *)malloc( neg_free * sizeof(double) );
		for (q = 0; q < neg_free; q++) {
			neg_listo[q] = 0;
		}
	}

	/* List parameters array */
	listprop = (int *)malloc(3 * sizeof(int));
	/* Make jagged list to store a maximum of 2^64 values */
	listval = (double **)malloc( 64 * sizeof(double *) );
	/* Initialize parameter list */
	initialize_list(listval, listprop);
	neg_listv = listval[listprop[1] - 1];

	/*(There are 3 pixel classes:
	 *  - frozen (processed)
	 *  - narrow band (boundary) (in list to check for the next pixel with smallest distance)
	 *  - far (not yet used)
	 */

	/* set all starting points to distance zero and frozen  */
	/* and add all neighbours of the starting points to narrow list  */
	for (z = 0; z < dims_sp[1]; z++) {
		/*starting point  */
		x = izsrcn;
		y = ixsrcn;
		XY_index = x + y * dims[0];

		/*Set starting point to frozen and distance to zero  */
		Frozen[XY_index] = 1;
		T[XY_index] = 0;
		if (Ed) {
			Y[XY_index] = 0;
		}
	}

	for (z = 0; z < dims_sp[1]; z++) {
		/*starting point  */
		x = izsrcn;
		y = ixsrcn;
		XY_index = x + y * dims[0];


		/* Add neigbours of starting points  */
		for (k = 0; k < 4; k++) {
			/*Location of neighbour  */
			i = x + ne[k]; j = y + ne[k + 4];
			IJ_index = i + j * dims[0];

			/*Check if current neighbour is not yet frozen and inside the
			 *picture  */
			if (isntfrozen2d(i, j, dims, Frozen)) {
				Tt = (1 / (max(F[IJ_index], eps)));

				Ty = 1;
				/*Update distance in neigbour list or add to neigbour list */
				if (T[IJ_index] > 0) {
					if (neg_listv[(int)T[IJ_index]] > Tt) {
						listupdate(listval, listprop, (int)T[IJ_index], Tt);
					}
					if (Ed) {
						neg_listo[(int)T[IJ_index]] = min(neg_listo[(int)T[IJ_index]], Ty);
					}
				} else {
					/*If running out of memory at a new block  */
					if (neg_pos >= neg_free) {
						neg_free += 100000;
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					list_add(listval, listprop, Tt);
					neg_listv = listval[listprop[1] - 1];
					neg_listx[neg_pos] = i;
					neg_listy[neg_pos] = j;
					if (Ed) {
						neg_listo[neg_pos] = Ty;
					}
					T[IJ_index] = neg_pos;
					neg_pos++;
				}
			}
		}
	}
	/*Loop through all pixels of the image  */
	for (itt = 0; itt < npixels; itt++) {
		/*Get the pixel from narrow list (boundary list) with smallest
		 *distance value and set it to current pixel location  */
		index = list_minimum(listval, listprop);
		neg_listv = listval[listprop[1] - 1];

		/* Stop if pixel distance is infinite (all pixels are processed)  */
		if (IsInf(neg_listv[index])) {
			break;
		}
		x = (int)neg_listx[index]; y = (int)neg_listy[index];

		XY_index = x + y * dims[0];
		Frozen[XY_index] = 1;
		T[XY_index] = neg_listv[index];
		if (Ed) {
			Y[XY_index] = neg_listo[index];
		}


		/*Remove min value by replacing it with the last value in the array  */
		list_remove_replace(listval, listprop, index) ;
		neg_listv = listval[listprop[1] - 1];
		if (index < (neg_pos - 1)) {
			neg_listx[index] = neg_listx[neg_pos - 1];
			neg_listy[index] = neg_listy[neg_pos - 1];
			if (Ed) {
				neg_listo[index] = neg_listo[neg_pos - 1];
			}
			T[(int)(neg_listx[index] + neg_listy[index]*dims[0])] = index;
		}
		neg_pos = neg_pos - 1;


		/*Loop through all 4 neighbours of current pixel  */
		for (k = 0; k < 4; k++) {

			/*Location of neighbour  */
			i = x + ne[k]; j = y + ne[k + 4];
			IJ_index = i + j * dims[0];

			/*Check if current neighbour is not yet frozen and inside the  */
			/*picture  */
			if (isntfrozen2d(i, j, dims, Frozen)) {
				// analyticaly compute the time table near the source - Dongzhuo
				// originally
				// Tt = CalculateDistance2D(T, F[IJ_index], dims, i, j, usesecond, usecross, Frozen);
				if (abs(izsrcn - i) <= bz && abs(ixsrcn - j) <= bx) {
					Tt = CalculateTConstant2D(aveVel, izsrcn, ixsrcn, i, j);
				} else {
					// Tt = CalculateDistance2D(T, F[IJ_index], dims, i, j, usesecond, usecross, Frozen);
					// should multiply the dimension of the grid - Dongzhuo
					Tt = CalculateDistance2D(T, F[IJ_index], dims, i, j, usesecond, usecross, Frozen);

				}

				if (Ed) {
					// should multiply the dimension of the grid - Dongzhuo
					Ty = CalculateDistance2D(Y, 1, dims, i, j, usesecond, usecross, Frozen);
				}

				/*Update distance in neigbour list or add to neigbour list */
				IJ_index = i + j * dims[0];
				if ((T[IJ_index] > -1) && T[IJ_index] <= listprop[0]) {
					if (neg_listv[(int)T[IJ_index]] > Tt) {
						listupdate(listval, listprop,    (int)T[IJ_index], Tt);
					}
					if (Ed) {
						neg_listo[neg_pos] = min(neg_listo[neg_pos], Ty);
					}
				} else {
					/*If running out of memory at a new block */
					if (neg_pos >= neg_free) {
						neg_free += 100000;
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					list_add(listval, listprop, Tt);
					neg_listv = listval[listprop[1] - 1];
					neg_listx[neg_pos] = i; neg_listy[neg_pos] = j;
					if (Ed) {
						neg_listo[neg_pos] = Ty;
					}
					T[IJ_index] = neg_pos;
					neg_pos++;
				}
			}
		}

	}
	/* Free memory */
	/* Destroy parameter list */
	destroy_list(listval, listprop);
	free(neg_listx);
	free(neg_listy);
	if (Ed) {
		free(neg_listo);
	}
	free(Frozen);
}

void msfm3dCpp(double *T, double *F, int izsrcn, int ixsrcn, int iysrcn, int *dims, bool usesecond, bool usecross) {
	/* The input variables */

	// box around the source - Dongzhuo
	int bz = 0;
	int bx = 0;
	int by = 0;

	/* Euclidian distance image */
	double *Y; // NOT USED CURRENTLY!!!

	/* Current distance values */
	double Tt, Ty;

	/* Matrix containing the Frozen Pixels" */
	bool *Frozen;

	/* Augmented Fast Marching (For skeletonize) */
	bool Ed = 0;

	/* Size of input image */
	// const mwSize *dims_c;
	// mwSize dims[3];

	/* Size of  SourcePoints array */
	// const mwSize *dims_sp_c;
	// mwSize dims_sp[3];

	/* Number of pixels in image */
	int npixels;
	int dims_sp[2] = {3, 1};

	/* Neighbour list */
	int neg_free;
	int neg_pos;
	double *neg_listv;
	double *neg_listx;
	double *neg_listy;
	double *neg_listz;
	double *neg_listo;

	int *listprop;
	double **listval;

	/* Neighbours 6x3 */
	int ne[18] = { -1,  0,  0, 1, 0, 0, 0, -1,  0, 0, 1, 0, 0,  0, -1, 0, 0, 1};

	/* Loop variables */
	int s, w, itt, q;

	/* Current location */
	int i, j, k, x, y, z;

	/* Index */
	int IJK_index, XYZ_index, index;

	npixels = dims[0] * dims[1] * dims[2];

	// compute average velocity in the small box around the soruce - Dongzhuo
	double aveVel = averageVelocityAroundSource3D(bz, bx, by, izsrcn, ixsrcn, iysrcn, dims, F);
	cout << "aveVel = " << aveVel << endl;

	/* Check for proper number of input and output arguments. */
	// if (nrhs < 3) {
	//  mexErrMsgTxt("2 to 4 inputs are required.");
	// }
	// if (nlhs == 1) {
	//  Ed = 0;
	// } else if (nlhs == 2) {
	//  Ed = 1;
	// } else {
	//  mexErrMsgTxt("One or Two outputs required");
	// }

	/* Check data input types */
	// if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) {
	//  mexErrMsgTxt("Speed image must be of class double");
	// }
	// if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) {
	//  mexErrMsgTxt("SourcePoints must be of class double");
	// }

	// if ((nrhs > 2) && (mxGetClassID(prhs[2]) != mxLOGICAL_CLASS)) {
	//  mexErrMsgTxt("UseSecond must be of class boolean / logical");
	// }

	// if ((nrhs > 3) && (mxGetClassID(prhs[3]) != mxLOGICAL_CLASS)) {
	//  mexErrMsgTxt("UseCross must be of class boolean / logical");
	// }

	/* Get the sizes of the input image volume */
	// if (mxGetNumberOfDimensions(prhs[0]) == 3) {
	//  dims_c = mxGetDimensions(prhs[0]);
	//  dims[0] = dims_c[0]; dims[1] = dims_c[1]; dims[2] = dims_c[2];
	//  npixels = dims[0] * dims[1] * dims[2];
	// } else {
	//  mexErrMsgTxt("Speed image must be 3d.");
	// }

	/* Get the sizes of the  SourcePoints */
	// dims_sp_c = mxGetDimensions(prhs[1]);

	// if (dims_sp_c[0] != 3) {
	//  mexErrMsgTxt("SourcePoints must be a 3xn matrix.");
	// }
	// dims_sp[0] = dims_sp_c[0]; dims_sp[1] = dims_sp_c[1]; dims_sp[2] = dims_sp_c[2];

	/* Get pointers/data from  to each input. */
	// F = (double *)mxGetPr(prhs[0]);
	// SourcePoints = (double *)mxGetPr(prhs[1]);
	// if (nrhs > 2) {
	//  useseconda = (bool *)mxGetPr(prhs[2]);
	//  usesecond = useseconda[0];
	// }
	// if (nrhs > 3) {
	//  usecrossa = (bool *)mxGetPr(prhs[3]);
	//  usecross = usecrossa[0];
	// }

	/* Create the distance output array */
	// plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	/* Assign pointer to output. */
	/*Distance image, also used to store the index of narrowband pixels */
	/*during marching process */
	// T = mxGetPr(plhs[0]);
	// if (Ed) {
	//  plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	//  Y = mxGetPr(plhs[1]);
	// }

	/* Pixels which are processed and have a final distance are frozen */
	Frozen = (bool *)malloc( npixels * sizeof(bool) );
	for (q = 0; q < npixels; q++) {
		Frozen[q] = 0;
		T[q] = -1;
	}
	if (Ed) {
		for (q = 0; q < npixels; q++) {
			Y[q] = -1;
		}
	}



	/*Free memory to store neighbours of the (segmented) region */
	neg_free = 100000;
	neg_pos = 0;

	neg_listx = (double *)malloc( neg_free * sizeof(double) );
	neg_listy = (double *)malloc( neg_free * sizeof(double) );
	neg_listz = (double *)malloc( neg_free * sizeof(double) );
	if (Ed) {
		neg_listo = (double *)malloc( neg_free * sizeof(double) );
		for (q = 0; q < neg_free; q++) {
			neg_listo[q] = 0;
		}
	}

	/* List parameters array */
	listprop = (int *)malloc(3 * sizeof(int));
	/* Make jagged list to store a maximum of 2^64 values */
	listval = (double **)malloc( 64 * sizeof(double *) );

	/* Initialize parameter list */
	initialize_list(listval, listprop);
	neg_listv = listval[listprop[1] - 1];


	/*(There are 3 pixel classes: */
	/*  - frozen (processed) */
	/*  - narrow band (boundary) (in list to check for the next pixel with smallest distance) */
	/*  - far (not yet used) */
	/* set all starting points to distance zero and frozen */
	/* and add all neighbours of the starting points to narrow list */

	for (s = 0; s < dims_sp[1]; s++) {
		/*starting point */
		x = izsrcn;
		y = ixsrcn;
		z = iysrcn;

		XYZ_index = mindex3(x, y, z, dims[0], dims[1]);


		Frozen[XYZ_index] = 1;
		T[XYZ_index] = 0;
		if (Ed) {
			Y[XYZ_index] = 0;
		}
	}

	for (s = 0; s < dims_sp[1]; s++) {
		/*starting point */
		x = izsrcn;
		y = ixsrcn;
		z = iysrcn;

		XYZ_index = mindex3(x, y, z, dims[0], dims[1]);


		for (w = 0; w < 6; w++) {
			/*Location of neighbour */
			i = x + ne[w];
			j = y + ne[w + 6];
			k = z + ne[w + 12];

			IJK_index = mindex3(i, j, k, dims[0], dims[1]);

			/*Check if current neighbour is not yet frozen and inside the */

			/*picture */
			if (isntfrozen3d(i, j, k, dims, Frozen)) {
				Tt = (1 / (max(F[IJK_index], eps)));
				/*Update distance in neigbour list or add to neigbour list */
				if (T[IJK_index] > 0) {
					if (neg_listv[(int)T[IJK_index]] > Tt) {
						listupdate(listval, listprop, (int)T[IJK_index], Tt);
					}
				} else {
					/*If running out of memory at a new block */
					if (neg_pos >= neg_free) {
						neg_free += 100000;
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						neg_listz = (double *)realloc(neg_listz, neg_free * sizeof(double) );
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					list_add(listval, listprop, Tt);
					neg_listv = listval[listprop[1] - 1];
					neg_listx[neg_pos] = i;
					neg_listy[neg_pos] = j;
					neg_listz[neg_pos] = k;
					T[IJK_index] = neg_pos;
					neg_pos++;
				}
			}
		}
	}

	/*Loop through all pixels of the image */
	for (itt = 0; itt < (npixels); itt++) { /* */
		/*Get the pixel from narrow list (boundary list) with smallest */
		/*distance value and set it to current pixel location */
		index = list_minimum(listval, listprop);
		neg_listv = listval[listprop[1] - 1];
		/* Stop if pixel distance is infinite (all pixels are processed) */
		if (IsInf(neg_listv[index])) {
			break;
		}

		/*index=minarray(neg_listv, neg_pos); */
		x = (int)neg_listx[index]; y = (int)neg_listy[index]; z = (int)neg_listz[index];
		XYZ_index = mindex3(x, y, z, dims[0], dims[1]);
		Frozen[XYZ_index] = 1;
		T[XYZ_index] = neg_listv[index];
		if (Ed) {
			Y[XYZ_index] = neg_listo[index];
		}

		/*Remove min value by replacing it with the last value in the array */
		list_remove_replace(listval, listprop, index) ;
		neg_listv = listval[listprop[1] - 1];
		if (index < (neg_pos - 1)) {
			neg_listx[index] = neg_listx[neg_pos - 1];
			neg_listy[index] = neg_listy[neg_pos - 1];
			neg_listz[index] = neg_listz[neg_pos - 1];
			if (Ed) {
				neg_listo[index] = neg_listo[neg_pos - 1];
			}
			T[(int)mindex3((int)neg_listx[index], (int)neg_listy[index], (int)neg_listz[index], dims[0], dims[1])] = index;
		}
		neg_pos = neg_pos - 1;

		/*Loop through all 6 neighbours of current pixel */
		for (w = 0; w < 6; w++) {
			/*Location of neighbour */
			i = x + ne[w]; j = y + ne[w + 6]; k = z + ne[w + 12];
			IJK_index = mindex3(i, j, k, dims[0], dims[1]);

			/*Check if current neighbour is not yet frozen and inside the */
			/*picture */
			if (isntfrozen3d(i, j, k, dims, Frozen)) {
				// analyticaly compute the time table near the source - Dongzhuo
				// originally
				// Tt = CalculateDistance3D(T, F[IJK_index], dims, i, j, k, usesecond, usecross, Frozen);
				if (abs(izsrcn - i) <= bz && abs(ixsrcn - j) <= bx && abs(iysrcn - i) <= by) {
					Tt = CalculateTConstant3D(aveVel, izsrcn, ixsrcn, iysrcn, i, j, k);
				} else {
					// Tt = CalculateDistance3D(T, F[IJK_index], dims, i, j, k, usesecond, usecross, Frozen);
					// should multiply the dimension of the grid - Dongzhuo
					Tt = CalculateDistance3D(T, F[IJK_index], dims, i, j, k, usesecond, usecross, Frozen);
				}

				if (Ed) {
					// should multiply the dimension of the grid - Dongzhuo
					Ty = CalculateDistance3D(Y, 1, dims, i, j, k, usesecond, usecross, Frozen);
				}

				/*Update distance in neigbour list or add to neigbour list */
				IJK_index = mindex3(i, j, k, dims[0], dims[1]);
				if ((T[IJK_index] > -1) && T[IJK_index] <= listprop[0]) {
					if (neg_listv[(int)T[IJK_index]] > Tt) {
						listupdate(listval, listprop, (int)T[IJK_index], Tt);
					}
				} else {
					/*If running out of memory at a new block */
					if (neg_pos >= neg_free) {
						neg_free += 100000;
						neg_listx = (double *)realloc(neg_listx, neg_free * sizeof(double) );
						neg_listy = (double *)realloc(neg_listy, neg_free * sizeof(double) );
						neg_listz = (double *)realloc(neg_listz, neg_free * sizeof(double) );
						if (Ed) {
							neg_listo = (double *)realloc(neg_listo, neg_free * sizeof(double) );
						}
					}
					list_add(listval, listprop, Tt);
					neg_listv = listval[listprop[1] - 1];
					neg_listx[neg_pos] = i; neg_listy[neg_pos] = j; neg_listz[neg_pos] = k;
					if (Ed) {
						neg_listo[neg_pos] = Ty;
					}

					T[IJK_index] = neg_pos;
					neg_pos++;
				}
			}
		}

	}
	/* Free memory */
	/* Destroy parameter list */
	destroy_list(listval, listprop);
	free(neg_listx);
	free(neg_listy);
	free(neg_listz);
	free(Frozen);
}

double CalculateDistance2D(double *T, double Fij, int *dims, int i, int j, bool usesecond, bool usecross, bool *Frozen) {
	/* Derivatives */
	double Tm[4] = {0, 0, 0, 0};
	double Tm2[4] = {0, 0, 0, 0};
	double Coeff[3];

	/* local derivatives in distance image */
	double Tpatch_2_3, Txm2, Tpatch_4_3, Txp2;
	double Tpatch_3_2, Tym2, Tpatch_3_4, Typ2;
	double Tpatch_2_2, Tr1m2, Tpatch_4_4, Tr1p2;
	double Tpatch_2_4, Tr2m2, Tpatch_4_2, Tr2p2;

	/* Return values root of polynomial */
	double ansroot[2] = {0, 0};

	/* Loop variables  */
	int q, t;

	/* Derivative checks */
	bool ch1, ch2;

	/* Order derivatives in a certain direction */
	int Order[4] = {0, 0, 0, 0};

	/* Current location */
	int in, jn;

	/* Constant cross term */
	const double c1 = 0.5;

	double Tt, Tt2;
	/*Get First order derivatives (only use frozen pixel)  */
	in = i - 1; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
		Tpatch_2_3 = T[mindex2(in, jn, dims[0])];
	} else {
		Tpatch_2_3 = INF;
	}
	in = i + 0; jn = j - 1; if (isfrozen2d(in, jn, dims, Frozen)) {
		Tpatch_3_2 = T[mindex2(in, jn, dims[0])];
	} else {
		Tpatch_3_2 = INF;
	}
	in = i + 0; jn = j + 1; if (isfrozen2d(in, jn, dims, Frozen)) {
		Tpatch_3_4 = T[mindex2(in, jn, dims[0])];
	} else {
		Tpatch_3_4 = INF;
	}
	in = i + 1; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
		Tpatch_4_3 = T[mindex2(in, jn, dims[0])];
	} else {
		Tpatch_4_3 = INF;
	}
	if (usecross) {
		in = i - 1; jn = j - 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_2_2 = T[mindex2(in, jn, dims[0])];
		} else {
			Tpatch_2_2 = INF;
		}
		in = i - 1; jn = j + 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_2_4 = T[mindex2(in, jn, dims[0])];
		} else {
			Tpatch_2_4 = INF;
		}
		in = i + 1; jn = j - 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_4_2 = T[mindex2(in, jn, dims[0])];
		} else {
			Tpatch_4_2 = INF;
		}
		in = i + 1; jn = j + 1; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tpatch_4_4 = T[mindex2(in, jn, dims[0])];
		} else {
			Tpatch_4_4 = INF;
		}
	}
	/*The values in order is 0 if no neighbours in that direction  */
	/*1 if 1e order derivatives is used and 2 if second order  */
	/*derivatives are used  */
	Order[0] = 0; Order[1] = 0; Order[2] = 0; Order[3] = 0;
	/*Make 1e order derivatives in x and y direction  */
	Tm[0] = min( Tpatch_2_3 , Tpatch_4_3); if (IsFinite(Tm[0])) {
		Order[0] = 1;
	}
	Tm[1] = min( Tpatch_3_2 , Tpatch_3_4); if (IsFinite(Tm[1])) {
		Order[1] = 1;
	}
	/*Make 1e order derivatives in cross directions  */
	if (usecross) {
		Tm[2] = min( Tpatch_2_2 , Tpatch_4_4); if (IsFinite(Tm[2])) {
			Order[2] = 1;
		}
		Tm[3] = min( Tpatch_2_4 , Tpatch_4_2); if (IsFinite(Tm[3])) {
			Order[3] = 1;
		}
	}

	/*Make 2e order derivatives  */
	if (usesecond) {
		/*Get Second order derivatives (only use frozen pixel) */
		in = i - 2; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
			Txm2 = T[mindex2(in, jn, dims[0])];
		} else {
			Txm2 = INF;
		}
		in = i + 2; jn = j + 0; if (isfrozen2d(in, jn, dims, Frozen)) {
			Txp2 = T[mindex2(in, jn, dims[0])];
		} else {
			Txp2 = INF;
		}
		in = i + 0; jn = j - 2; if (isfrozen2d(in, jn, dims, Frozen)) {
			Tym2 = T[mindex2(in, jn, dims[0])];
		} else {
			Tym2 = INF;
		}
		in = i + 0; jn = j + 2; if (isfrozen2d(in, jn, dims, Frozen)) {
			Typ2 = T[mindex2(in, jn, dims[0])];
		} else {
			Typ2 = INF;
		}
		if (usecross) {
			in = i - 2; jn = j - 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr1m2 = T[mindex2(in, jn, dims[0])];
			} else {
				Tr1m2 = INF;
			}
			in = i - 2; jn = j + 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr2m2 = T[mindex2(in, jn, dims[0])];
			} else {
				Tr2m2 = INF;
			}
			in = i + 2; jn = j - 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr2p2 = T[mindex2(in, jn, dims[0])];
			} else {
				Tr2p2 = INF;
			}
			in = i + 2; jn = j + 2; if (isfrozen2d(in, jn, dims, Frozen)) {
				Tr1p2 = T[mindex2(in, jn, dims[0])];
			} else {
				Tr1p2 = INF;
			}
		}

		Tm2[0] = 0; Tm2[1] = 0; Tm2[2] = 0; Tm2[3] = 0;
		/*pixels with a pixeldistance 2 from the center must be */
		/*lower in value otherwise use other side or first order */
		ch1 = (Txm2 < Tpatch_2_3) && IsFinite(Tpatch_2_3); ch2 = (Txp2 < Tpatch_4_3) && IsFinite(Tpatch_4_3);
		if (ch1 && ch2) {
			Tm2[0] = min( (4.0 * Tpatch_2_3 - Txm2) / 3.0 , (4.0 * Tpatch_4_3 - Txp2) / 3.0);  Order[0] = 2;
		} else if (ch1) {
			Tm2[0] = (4.0 * Tpatch_2_3 - Txm2) / 3.0; Order[0] = 2;
		} else if (ch2) {
			Tm2[0] = (4.0 * Tpatch_4_3 - Txp2) / 3.0; Order[0] = 2;
		}

		ch1 = (Tym2 < Tpatch_3_2) && IsFinite(Tpatch_3_2); ch2 = (Typ2 < Tpatch_3_4) && IsFinite(Tpatch_3_4);

		if (ch1 && ch2) {
			Tm2[1] = min( (4.0 * Tpatch_3_2 - Tym2) / 3.0 , (4.0 * Tpatch_3_4 - Typ2) / 3.0); Order[1] = 2;
		} else if (ch1) {
			Tm2[1] = (4.0 * Tpatch_3_2 - Tym2) / 3.0; Order[1] = 2;
		} else if (ch2) {
			Tm2[1] = (4.0 * Tpatch_3_4 - Typ2) / 3.0; Order[1] = 2;
		}
		if (usecross) {
			ch1 = (Tr1m2 < Tpatch_2_2) && IsFinite(Tpatch_2_2); ch2 = (Tr1p2 < Tpatch_4_4) && IsFinite(Tpatch_4_4);
			if (ch1 && ch2) {
				Tm2[2] = min( (4.0 * Tpatch_2_2 - Tr1m2) / 3.0 , (4.0 * Tpatch_4_4 - Tr1p2) / 3.0); Order[2] = 2;
			} else if (ch1) {
				Tm2[2] = (4.0 * Tpatch_2_2 - Tr1m2) / 3.0; Order[2] = 2;
			} else if (ch2) {
				Tm2[2] = (4.0 * Tpatch_4_4 - Tr1p2) / 3.0; Order[2] = 2;
			}

			ch1 = (Tr2m2 < Tpatch_2_4) && IsFinite(Tpatch_2_4); ch2 = (Tr2p2 < Tpatch_4_2) && IsFinite(Tpatch_4_2);
			if (ch1 && ch2) {
				Tm2[3] = min( (4.0 * Tpatch_2_4 - Tr2m2) / 3.0 , (4.0 * Tpatch_4_2 - Tr2p2) / 3.0); Order[3] = 2;
			} else if (ch1) {
				Tm2[3] = (4.0 * Tpatch_2_4 - Tr2m2) / 3.0; Order[3] = 2;
			} else if (ch2) {
				Tm2[3] = (4.0 * Tpatch_4_2 - Tr2p2) / 3.0; Order[3] = 2;
			}
		}
	}
	/*Calculate the distance using x and y direction */
	Coeff[0] = 0; Coeff[1] = 0; Coeff[2] = -1 / (max(pow2(Fij), eps));

	for (t = 0; t < 2; t++) {
		switch (Order[t]) {
		case 1:
			Coeff[0] += 1; Coeff[1] += -2 * Tm[t]; Coeff[2] += pow2(Tm[t]);
			break;
		case 2:
			Coeff[0] += (2.2500); Coeff[1] += -2.0 * Tm2[t] * (2.2500); Coeff[2] += pow2(Tm2[t]) * (2.2500);
			break;
		}
	}
	roots(Coeff, ansroot);
	Tt = max(ansroot[0], ansroot[1]);
	/*Calculate the distance using the cross directions */
	if (usecross) {
		/* Original Equation */
		/*    Coeff[0]=0; Coeff[1]=0; Coeff[2]=-1/(max(pow2(Fij),eps)) */
		Coeff[0] += 0; Coeff[1] += 0; Coeff[2] += -1 / (max(pow2(Fij), eps));
		for (t = 2; t < 4; t++) {
			switch (Order[t]) {
			case 1:
				Coeff[0] += c1; Coeff[1] += -2.0 * c1 * Tm[t]; Coeff[2] += c1 * pow2(Tm[t]);
				break;
			case 2:
				Coeff[0] += c1 * 2.25; Coeff[1] += -2 * c1 * Tm2[t] * (2.25); Coeff[2] += pow2(Tm2[t]) * c1 * 2.25;
				break;
			}
		}
		if (Coeff[0] > 0) {
			roots(Coeff, ansroot);
			Tt2 = max(ansroot[0], ansroot[1]);
			/*Select minimum distance value of both stensils */
			Tt = min(Tt, Tt2);
		}
	}
	/*Upwind condition check, current distance must be larger */
	/*then direct neighbours used in solution */
	/*(Will this ever happen?) */
	if (usecross) {
		for (q = 0; q < 4; q++) {
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 4)] + (1 / (max(Fij, eps)));
			}
		}
	} else {
		for (q = 0; q < 2; q++) {
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 2)] + (1 / (max(Fij, eps)));
			}
		}
	}
	return Tt;
}

double CalculateDistance3D(double *T, double Fijk, int *dims, int i, int j, int k, bool usesecond, bool usecross, bool *Frozen) {

	/* Loop variables */
	int q, t;

	/* Current location */
	int in, jn, kn;

	/* Derivatives */
	double Tm[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double Tm2[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double Coeff[3];

	/* local derivatives in distance image */
	double Txm1, Txm2, Txp1, Txp2;
	double Tym1, Tym2, Typ1, Typ2;
	double Tzm1, Tzm2, Tzp1, Tzp2;
	/* local cross derivatives in distance image */
	double Tr2t1m1, Tr2t1m2, Tr2t1p1, Tr2t1p2;
	double Tr2t2m1, Tr2t2m2, Tr2t2p1, Tr2t2p2;
	double Tr2t3m1, Tr2t3m2, Tr2t3p1, Tr2t3p2;
	double Tr3t1m1, Tr3t1m2, Tr3t1p1, Tr3t1p2;
	double Tr3t2m1, Tr3t2m2, Tr3t2p1, Tr3t2p2;
	double Tr3t3m1, Tr3t3m2, Tr3t3p1, Tr3t3p2;
	double Tr4t1m1, Tr4t1m2, Tr4t1p1, Tr4t1p2;
	double Tr4t2m1, Tr4t2m2, Tr4t2p1, Tr4t2p2;
	double Tr4t3m1, Tr4t3m2, Tr4t3p1, Tr4t3p2;
	double Tr5t1m1, Tr5t1m2, Tr5t1p1, Tr5t1p2;
	double Tr5t2m1, Tr5t2m2, Tr5t2p1, Tr5t2p2;
	double Tr5t3m1, Tr5t3m2, Tr5t3p1, Tr5t3p2;
	double Tr6t1m1, Tr6t1m2, Tr6t1p1, Tr6t1p2;
	double Tr6t2m1, Tr6t2m2, Tr6t2p1, Tr6t2p2;
	double Tr6t3m1, Tr6t3m2, Tr6t3p1, Tr6t3p2;

	double Tt, Tt2;


	/* Return values root of polynomial */
	double ansroot[2] = {0, 0};

	/* Order derivatives in a certain direction */
	int Order[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	/* Neighbours 4x2 */
	int ne[18] = { -1,  0,  0, 1, 0, 0, 0, -1,  0, 0, 1, 0, 0,  0, -1, 0, 0, 1};

	/* Stencil constants */
	double G1[18] = {1, 1, 1, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.3333333333333, 0.3333333333333, 0.5, 0.3333333333333, 0.3333333333333};
	double G2[18] = {2.250, 2.250, 2.250, 2.250, 1.125, 1.125, 2.250, 1.125, 1.125, 2.250, 1.125, 1.125, 1.125, 0.750, 0.750, 1.125, 0.750, 0.750};


	/*Get First order derivatives (only use frozen pixel) */
	in = i - 1; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Txm1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Txm1 = INF;
	}
	in = i + 1; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Txp1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Txp1 = INF;
	}
	in = i + 0; jn = j - 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Tym1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Tym1 = INF;
	}
	in = i + 0; jn = j + 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Typ1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Typ1 = INF;
	}
	in = i + 0; jn = j + 0; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Tzm1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Tzm1 = INF;
	}
	in = i + 0; jn = j + 0; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
		Tzp1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
	} else {
		Tzp1 = INF;
	}

	if (usecross) {
		Tr2t1m1 = Txm1;
		Tr2t1p1 = Txp1;
		in = i - 0; jn = j - 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t2m1 = INF;
		}
		in = i + 0; jn = j + 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t2p1 = INF;
		}
		in = i - 0; jn = j - 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t3m1 = INF;
		}
		in = i + 0; jn = j + 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr2t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr2t3p1 = INF;
		}
		Tr3t1m1 = Tym1;
		Tr3t1p1 = Typ1;
		in = i - 1; jn = j + 0; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t2m1 = INF;
		}
		in = i + 1; jn = j + 0; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t2p1 = INF;
		}
		in = i - 1; jn = j - 0; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t3m1 = INF;
		}
		in = i + 1; jn = j + 0; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr3t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr3t3p1 = INF;
		}
		Tr4t1m1 = Tzm1;
		Tr4t1p1 = Tzp1;
		in = i - 1; jn = j - 1; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t2m1 = INF;
		}
		in = i + 1; jn = j + 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t2p1 = INF;
		}
		in = i - 1; jn = j + 1; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t3m1 = INF;
		}
		in = i + 1; jn = j - 1; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr4t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr4t3p1 = INF;
		}
		Tr5t1m1 = Tr3t3m1;
		Tr5t1p1 = Tr3t3p1;
		in = i - 1; jn = j - 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t2m1 = INF;
		}
		in = i + 1; jn = j + 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t2p1 = INF;
		}
		in = i - 1; jn = j + 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t3m1 = INF;
		}
		in = i + 1; jn = j - 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr5t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr5t3p1 = INF;
		}
		Tr6t1m1 = Tr3t2p1;
		Tr6t1p1 = Tr3t2m1;
		in = i - 1; jn = j - 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t2m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t2m1 = INF;
		}
		in = i + 1; jn = j + 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t2p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t2p1 = INF;
		}
		in = i - 1; jn = j + 1; kn = k - 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t3m1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t3m1 = INF;
		}
		in = i + 1; jn = j - 1; kn = k + 1; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tr6t3p1 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tr6t3p1 = INF;
		}
	}

	/*The values in order is 0 if no neighbours in that direction */
	/*1 if 1e order derivatives is used and 2 if second order */
	/*derivatives are used */


	/*Make 1e order derivatives in x and y direction */
	Tm[0] = min( Txm1 , Txp1); if (IsFinite(Tm[0])) {
		Order[0] = 1;
	} else {
		Order[0] = 0;
	}
	Tm[1] = min( Tym1 , Typ1); if (IsFinite(Tm[1])) {
		Order[1] = 1;
	} else {
		Order[1] = 0;
	}
	Tm[2] = min( Tzm1 , Tzp1); if (IsFinite(Tm[2])) {
		Order[2] = 1;
	} else {
		Order[2] = 0;
	}

	/*Make 1e order derivatives in cross directions */
	if (usecross) {
		Tm[3] = Tm[0]; Order[3] = Order[0];
		Tm[4] = min( Tr2t2m1 , Tr2t2p1); if (IsFinite(Tm[4])) {
			Order[4] = 1;
		} else {
			Order[4] = 0;
		}
		Tm[5] = min( Tr2t3m1 , Tr2t3p1); if (IsFinite(Tm[5])) {
			Order[5] = 1;
		} else {
			Order[5] = 0;
		}

		Tm[6] = Tm[1]; Order[6] = Order[1];
		Tm[7] = min( Tr3t2m1 , Tr3t2p1); if (IsFinite(Tm[7])) {
			Order[7] = 1;
		} else {
			Order[7] = 0;
		}
		Tm[8] = min( Tr3t3m1 , Tr3t3p1); if (IsFinite(Tm[8])) {
			Order[8] = 1;
		} else {
			Order[8] = 0;
		}

		Tm[9] = Tm[2]; Order[9] = Order[2];
		Tm[10] = min( Tr4t2m1 , Tr4t2p1); if (IsFinite(Tm[10])) {
			Order[10] = 1;
		} else {
			Order[10] = 0;
		}
		Tm[11] = min( Tr4t3m1 , Tr4t3p1); if (IsFinite(Tm[11])) {
			Order[11] = 1;
		} else {
			Order[11] = 0;
		}

		Tm[12] = Tm[8]; Order[12] = Order[8];
		Tm[13] = min( Tr5t2m1 , Tr5t2p1); if (IsFinite(Tm[13])) {
			Order[13] = 1;
		} else {
			Order[13] = 0;
		}
		Tm[14] = min( Tr5t3m1 , Tr5t3p1); if (IsFinite(Tm[14])) {
			Order[14] = 1;
		} else {
			Order[14] = 0;
		}

		Tm[15] = Tm[7]; Order[15] = Order[7];
		Tm[16] = min( Tr6t2m1 , Tr6t2p1); if (IsFinite(Tm[16])) {
			Order[16] = 1;
		} else {
			Order[16] = 0;
		}
		Tm[17] = min( Tr6t3m1 , Tr6t3p1); if (IsFinite(Tm[17])) {
			Order[17] = 1;
		} else {
			Order[17] = 0;
		}
	}

	/*Make 2e order derivatives */
	if (usesecond) {
		/*Get Second order derivatives (only use frozen pixel) */
		/*Get First order derivatives (only use frozen pixel) */
		in = i - 2; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Txm2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Txm2 = INF;
		}
		in = i + 2; jn = j + 0; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Txp2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Txp2 = INF;
		}
		in = i + 0; jn = j - 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tym2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tym2 = INF;
		}
		in = i + 0; jn = j + 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Typ2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Typ2 = INF;
		}
		in = i + 0; jn = j + 0; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tzm2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tzm2 = INF;
		}
		in = i + 0; jn = j + 0; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
			Tzp2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
		} else {
			Tzp2 = INF;
		}

		if (usecross) {
			Tr2t1m2 = Txm2;
			Tr2t1p2 = Txp2;
			in = i - 0; jn = j - 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t2m2 = INF;
			}
			in = i + 0; jn = j + 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t2p2 = INF;
			}
			in = i - 0; jn = j - 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t3m2 = INF;
			}
			in = i + 0; jn = j + 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr2t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr2t3p2 = INF;
			}
			Tr3t1m2 = Tym2;
			Tr3t1p2 = Typ2;
			in = i - 2; jn = j + 0; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t2m2 = INF;
			}
			in = i + 2; jn = j + 0; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t2p2 = INF;
			}
			in = i - 2; jn = j - 0; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t3m2 = INF;
			}
			in = i + 2; jn = j + 0; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr3t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr3t3p2 = INF;
			}
			Tr4t1m2 = Tzm2;
			Tr4t1p2 = Tzp2;
			in = i - 2; jn = j - 2; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t2m2 = INF;
			}
			in = i + 2; jn = j + 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t2p2 = INF;
			}
			in = i - 2; jn = j + 2; kn = k - 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t3m2 = INF;
			}
			in = i + 2; jn = j - 2; kn = k + 0; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr4t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr4t3p2 = INF;
			}
			Tr5t1m2 = Tr3t3m2;
			Tr5t1p2 = Tr3t3p2;
			in = i - 2; jn = j - 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t2m2 = INF;
			}
			in = i + 2; jn = j + 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t2p2 = INF;
			}
			in = i - 2; jn = j + 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t3m2 = INF;
			}
			in = i + 2; jn = j - 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr5t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr5t3p2 = INF;
			}
			Tr6t1m2 = Tr3t2p2;
			Tr6t1p2 = Tr3t2m2;
			in = i - 2; jn = j - 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t2m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr6t2m2 = INF;
			}
			in = i + 2; jn = j + 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t2p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr6t2p2 = INF;
			}
			in = i - 2; jn = j + 2; kn = k - 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t3m2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr6t3m2 = INF;
			}
			in = i + 2; jn = j - 2; kn = k + 2; if (isfrozen3d(in, jn, kn, dims, Frozen)) {
				Tr6t3p2 = T[mindex3(in, jn, kn, dims[0], dims[1])];
			} else {
				Tr6t3p2 = INF;
			}
		}


		/*pixels with a pixeldistance 2 from the center must be */
		/*lower in value otherwise use other side or first order */

		Tm2[0] = second_derivative(Txm1, Txm2, Txp1, Txp2); if (IsInf(Tm2[0])) {
			Tm2[0] = 0;
		} else {
			Order[0] = 2;
		}
		Tm2[1] = second_derivative(Tym1, Tym2, Typ1, Typ2); if (IsInf(Tm2[1])) {
			Tm2[1] = 0;
		} else {
			Order[1] = 2;
		}
		Tm2[2] = second_derivative(Tzm1, Tzm2, Tzp1, Tzp2); if (IsInf(Tm2[2])) {
			Tm2[2] = 0;
		} else {
			Order[2] = 2;
		}

		if (usecross) {
			Tm2[3] = Tm2[0]; Order[3] = Order[0];
			Tm2[4] = second_derivative(Tr2t2m1, Tr2t2m2, Tr2t2p1, Tr2t2p2); if (IsInf(Tm2[4])) {
				Tm2[4] = 0;
			} else {
				Order[4] = 2;
			}
			Tm2[5] = second_derivative(Tr2t3m1, Tr2t3m2, Tr2t3p1, Tr2t3p2); if (IsInf(Tm2[5])) {
				Tm2[5] = 0;
			} else {
				Order[5] = 2;
			}

			Tm2[6] = Tm2[1]; Order[6] = Order[1];
			Tm2[7] = second_derivative(Tr3t2m1, Tr3t2m2, Tr3t2p1, Tr3t2p2); if (IsInf(Tm2[7])) {
				Tm2[7] = 0;
			} else {
				Order[7] = 2;
			}
			Tm2[8] = second_derivative(Tr3t3m1, Tr3t3m2, Tr3t3p1, Tr3t3p2); if (IsInf(Tm2[8])) {
				Tm2[8] = 0;
			} else {
				Order[8] = 2;
			}

			Tm2[9] = Tm2[2]; Order[9] = Order[2];
			Tm2[10] = second_derivative(Tr4t2m1, Tr4t2m2, Tr4t2p1, Tr4t2p2); if (IsInf(Tm2[10])) {
				Tm2[10] = 0;
			} else {
				Order[10] = 2;
			}
			Tm2[11] = second_derivative(Tr4t3m1, Tr4t3m2, Tr4t3p1, Tr4t3p2); if (IsInf(Tm2[11])) {
				Tm2[11] = 0;
			} else {
				Order[11] = 2;
			}

			Tm2[12] = Tm2[8]; Order[12] = Order[8];
			Tm2[13] = second_derivative(Tr5t2m1, Tr5t2m2, Tr5t2p1, Tr5t2p2); if (IsInf(Tm2[13])) {
				Tm2[13] = 0;
			} else {
				Order[13] = 2;
			}
			Tm2[14] = second_derivative(Tr5t3m1, Tr5t3m2, Tr5t3p1, Tr5t3p2); if (IsInf(Tm2[14])) {
				Tm2[14] = 0;
			} else {
				Order[14] = 2;
			}

			Tm2[15] = Tm2[7]; Order[15] = Order[7];
			Tm2[16] = second_derivative(Tr6t2m1, Tr6t2m2, Tr6t2p1, Tr6t2p2); if (IsInf(Tm2[16])) {
				Tm2[16] = 0;
			} else {
				Order[16] = 2;
			}
			Tm2[17] = second_derivative(Tr6t3m1, Tr6t3m2, Tr6t3p1, Tr6t3p2); if (IsInf(Tm2[17])) {
				Tm2[17] = 0;
			} else {
				Order[17] = 2;
			}
		}

	}

	/*Calculate the distance using x and y direction */
	Coeff[0] = 0; Coeff[1] = 0; Coeff[2] = -1 / (max(pow2(Fijk), eps));

	for (t = 0; t < 3; t++) {
		switch (Order[t]) {
		case 1:
			Coeff[0] += G1[t]; Coeff[1] += -2.0 * Tm[t] * G1[t]; Coeff[2] += pow2(Tm[t]) * G1[t];
			break;
		case 2:
			Coeff[0] += G2[t]; Coeff[1] += -2.0 * Tm2[t] * G2[t]; Coeff[2] += pow2(Tm2[t]) * G2[t];
			break;
		}
	}


	roots(Coeff, ansroot);
	Tt = max(ansroot[0], ansroot[1]);

	/*Calculate the distance using the cross directions */
	if (usecross) {

		for (q = 1; q < 6; q++) {
			/* Original Equation */
			/*    Coeff[0]=0; Coeff[1]=0; Coeff[2]=-1/(max(pow2(Fijk),eps)) */
			Coeff[0] += 0; Coeff[1] += 0; Coeff[2] += -1 / (max(pow2(Fijk), eps));

			for (t = q * 3; t < ((q + 1) * 3); t++) {
				switch (Order[t]) {
				case 1:
					Coeff[0] += G1[t]; Coeff[1] += -2.0 * Tm[t] * G1[t]; Coeff[2] += pow2(Tm[t]) * G1[t];
					break;
				case 2:
					Coeff[0] += G2[t]; Coeff[1] += -2.0 * Tm2[t] * G2[t]; Coeff[2] += pow2(Tm2[t]) * G2[t];
					break;
				}
			}
			/*Select maximum root solution and minimum distance value of both stensils */
			if (Coeff[0] > 0) {
				roots(Coeff, ansroot);
				Tt2 = max(ansroot[0], ansroot[1]);
				Tt = min(Tt, Tt2);
			}
		}
	}

	/*Upwind condition check, current distance must be larger */
	/*then direct neighbours used in solution */
	/*(Will this ever happen?) */
	if (usecross) {
		for (q = 0; q < 18; q++) {
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 18)] + (1 / (max(Fijk, eps)));
			}
		}
	} else {
		for (q = 0; q < 3; q++) {
			if (IsFinite(Tm[q]) && (Tt < Tm[q])) {
				Tt = Tm[minarray(Tm, 3)] + (1 / (max(Fijk, eps)));
			}
		}
	}

	return Tt;
}

double second_derivative(double Txm1, double Txm2, double Txp1, double Txp2) {
	bool ch1, ch2;
	double Tm;
	Tm = INF;
	ch1 = (Txm2 < Txm1) && IsFinite(Txm1); ch2 = (Txp2 < Txp1) && IsFinite(Txp1);
	if (ch1 && ch2) {
		Tm = min( (4.0 * Txm1 - Txm2) / 3.0 , (4.0 * Txp1 - Txp2) / 3.0);
	} else if (ch1) {
		Tm = (4.0 * Txm1 - Txm2) / 3.0;
	} else if (ch2) {
		Tm = (4.0 * Txp1 - Txp2) / 3.0;
	}
	return Tm;
}




// calculate the average velocity near the source, note that the indexs are different from the original code. Here
// we have the z, x, y system
// Dongzhuo June 5 2015
double averageVelocityAroundSource2D(int bz, int bx, int izsrcn, int ixsrcn, int *dims, double *F) {
	double count = 0.0;
	double velSum = 0.0;
	for (int i = 0; i <= 2 * bx; i++ ) {
		for (int k = 0; k <= 2 * bz; k++) {
			if (k >= 0 && k <= dims[0] && i >= 0 && i <= dims[1]) {
				velSum = velSum + F[mindex2(izsrcn + (k - bz), ixsrcn + (i - bx), dims[0])];
				count ++;
			}
		}
	}
	double aveVel = velSum / count;
	return aveVel;
}


// calculate time table near source based upon the average veocity in the box
// Dongzhuo June 5 2015
double CalculateTConstant2D(double aveVel, int izsrcn, int ixsrcn, int zloc, int xloc) {
	return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2)) / aveVel;
}

// calculate the average velocity near the source, note that the indexs are different from the original code. Here
// we have the z, x, y system
// Dongzhuo June 5 2015
double averageVelocityAroundSource3D(int bz, int bx, int by, int izsrcn, int ixsrcn, int iysrcn, int *dims, double *F) {
	double count = 0.0;
	double velSum = 0.0;
	for (int j = 0; j <= 2 * by; j++) {
		for (int i = 0; i <= 2 * bx; i++ ) {
			for (int k = 0; k <= 2 * bz; k++) {
				if (k >= 0 && k <= dims[0] && i >= 0 && i <= dims[1] && j >= 0 && j <= dims[2]) {
					velSum = velSum + F[mindex3(izsrcn + (k - bz), ixsrcn + (i - bx), iysrcn + (j - by), dims[0], dims[1])];
					count ++;
				}
			}
		}
	}
	double aveVel = velSum / count;
	return aveVel;
}

// calculate time table near source based upon the average veocity in the box
// Dongzhuo June 5 2015
double CalculateTConstant3D(double aveVel, int izsrcn, int ixsrcn, int iysrcn, int zloc, int xloc, int yloc) {
	return sqrt(pow((zloc - izsrcn), 2) + pow((xloc - ixsrcn), 2) + pow((yloc - iysrcn), 2)) / aveVel;
}