// ---------------
//   Dongzhuo Li & Ping Tong
//   May, 2015
// ---------------


#include <cmath>

void localDimension(double xsrc, double ysrc, double zsrc, double xrec, double yrec, double zrec, double dxyz_loc, int &dim_x, int &dim_y, int &dim_z) {
	double distance_x = sqrt((xsrc - xrec) * (xsrc - xrec) + (ysrc - yrec) * (ysrc - yrec));
	double distance_z = fabs(zsrc - zrec);

	dim_x = int(fmax(distance_x / dxyz_loc + 10, 10));
	if (dim_x > 1000) dim_x = 1000;

	dim_z = int(fmax(distance_z / dxyz_loc + 10, 10));
	if (dim_z > 1000) dim_z = 1000;

	dim_y = 50;

}