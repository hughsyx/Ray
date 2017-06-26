// ---------------
//   Dongzhuo Li & Ping Tong
//   May, 2015
// ---------------
#include <iostream>
#include <cmath>

double Vglb2loc(double *velModel, int nx, int ny, int nz, double xmin, double ymin, double zmin, \
                double xmax, double ymax, double zmax, double dxyz_glb, double x_glb, double y_glb, double z_glb);

using namespace std;

void coordinateTransform(int nx, int ny, int nz, double xmin, double ymin, double zmin, \
                         double dxyz_glb, double xsrc, double ysrc, double zsrc, double xrec, double yrec, \
                         double zrec, int dim_x, int dim_y, int dim_z, double dxyz_loc, double *velModel, \
                         double &loc_xsrc, double &loc_ysrc, double &loc_zsrc, double &loc_xrec, double &loc_yrec, \
                         double &loc_zrec, int &ixsrcn, int &iysrcn, int &izsrcn, \
                         double  &xmin_loc, double &ymin_loc, double &zmin_loc, \
                         double &theta, double *locVelModel) {

#define LOCVEL(ix,iy,iz) locVelModel[iy*dim_z*dim_x+ix*dim_z+iz]
	// double  PI = acos(-1.0);

	loc_xsrc = 0.0;
	loc_ysrc = 0.0;
	loc_zsrc = 0.0;
	loc_xrec = sqrt((xsrc - xrec) * (xsrc - xrec) + (ysrc - yrec) * (ysrc - yrec));
	loc_yrec = 0.0;
	loc_zrec = zrec - zsrc;

	// izsrcn = 4;
	ixsrcn = 4; // in local coordinate, the source is at the 5th grid point in x direction
	iysrcn = 25; // similar

	xmin_loc = -4.0 * dxyz_loc;
	ymin_loc = -25.0 * dxyz_loc;

	double x_local, y_local, z_local, x_glb, y_glb, z_glb;
	double xmax = xmin + nx * dxyz_glb;
	double ymax = ymin + ny * dxyz_glb;
	double zmax = zmin + nz * dxyz_glb;

	theta = atan2(yrec - ysrc, xrec - xsrc); // be careful about  last term ==0

	for (int i = 0; i < dim_x; i++) {
		x_local = (i - 4) * dxyz_loc; // the 5th grid is the source
		for (int j = 0; j < dim_y; j++) {
			y_local = (j - 25) * dxyz_loc;

			for (int k = 0; k < dim_z; k++) {
				if (loc_zrec >= 0) {
					z_local = (k - 4) * dxyz_loc;
					izsrcn = 4;
					zmin_loc = -4.0 * dxyz_loc;
				} else {
					z_local = (k - dim_z + 4) * dxyz_loc;
					izsrcn = dim_z - 4;
					zmin_loc = (4 - dim_z) * dxyz_loc;
				}

				x_glb = cos(theta) * x_local - sin(theta) * y_local + xsrc;
				y_glb = sin(theta) * x_local + cos(theta) * y_local + ysrc;
				z_glb = z_local + zsrc;

				// cout << "x_glb = " << x_glb << "y_glb = " << y_glb << "z_glb = " << endl;
				LOCVEL(i, j, k) = Vglb2loc(velModel, nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax, dxyz_glb, x_glb, y_glb, z_glb);
			}
		}
	}
}


double Vglb2loc(double *velModel, int nx, int ny, int nz, double xmin, double ymin, double zmin, \
                double xmax, double ymax, double zmax, double dxyz_glb, double x_glb, double y_glb, double z_glb) {
#define VEL(ix,iy,iz) velModel[iy*nz*nx+ix*nz+iz]

	int i1, j1, k1, i2, j2, k2;
	double wx, wy, wz;

	if (x_glb < xmin || y_glb < ymin || z_glb < zmin || x_glb > xmax || y_glb > ymax || z_glb > zmax) {
		return 0.1;
	}

	i1 = floor((x_glb - xmin) / dxyz_glb);
	i2 = i1 + 1;

	j1 = floor((y_glb - ymin) / dxyz_glb);
	j2 = j1 + 1;

	k1 = floor((z_glb - zmin) / dxyz_glb);
	k2 = k1 + 1;

	wx = 1.0 - (x_glb - i1 * dxyz_glb - xmin) / dxyz_glb;
	wy = 1.0 - (y_glb - i1 * dxyz_glb - ymin) / dxyz_glb;
	wz = 1.0 - (z_glb - i1 * dxyz_glb - zmin) / dxyz_glb;

	return VEL(i1, j1, k1) * wx * wy * wz \
	       + VEL(i2, j1, k1) * (1 - wx) * wy * wz \
	       + VEL(i2, j2, k1) * (1 - wx) * (1 - wy) * wz \
	       + VEL(i1, j2, k1) * wx * (1 - wy) * wz \
	       + VEL(i1, j1, k2) * wx * wy * (1 - wz) \
	       + VEL(i2, j1, k2) * (1 - wx) * wy * (1 - wz) \
	       + VEL(i2, j2, k2) * (1 - wx) * (1 - wy) * (1 - wz) \
	       + VEL(i1, j2, k2) * wx * (1 - wy) * (1 - wz);
}

