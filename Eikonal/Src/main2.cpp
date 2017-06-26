// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include "par_model_source_receiver.h"
#include "eikonal.h"
#include "raytracing.h"

using namespace std;

int main() {
	int dims[3] = {0, 0, 0};
	int dim_x, dim_y, dim_z;
	// parameters in Par_file
	int dimension, nx, ny, nz, invnx, invny, invnz, lengthG, ixsrcn, iysrcn, izsrcn;
	double dx, dy, dz, xmin, ymin, zmin, invdx, invdy, invdz, invxmin, invymin, invzmin;
	double dxyz_loc, loc_xsrc, loc_ysrc, loc_zsrc, loc_xrec, loc_yrec, loc_zrec, theta, xmin_loc, ymin_loc, zmin_loc;
	string velModel_file, source_file, receiver_file;
	double *velModel, *T;
	double *locVelModel;
	bool usesecond = true;
	bool usecross = true;
	ofstream timeMap, rayPointsOutPut;
	vector <double *> SourcePoints;
	map < int, vector< vector <double> > > ReceiverPoints;
	vector< vector<double> > rayPoints;

	// ------------------------ Initializing ---------------------------------------
	dim_x = dim_y = dim_z = 0;
	dimension = nx = ny = nz = invnx = invny = invnz = ixsrcn = iysrcn = izsrcn = 0;
	dx = dy = dz = xmin = ymin = zmin = invdx = invdy = invdz = invxmin = invymin = invzmin = 0.0;
	loc_xsrc = loc_ysrc = loc_zsrc = loc_xrec = loc_yrec = loc_zrec = theta = xmin_loc = ymin_loc = zmin_loc = 0.0;


	// ------------------------ Read Parameters & Load Model ----------------------
	read_parafile(dimension, nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, invnx, invny, invnz, \
	              invdx, invdy, invdz, invxmin, invymin, invzmin, velModel_file, source_file, receiver_file);

	dxyz_loc = dx / 10; // Change this later!!!

	if (dimension == 2) {
		dims[0] = nz; dims[1] = nx;
		lengthG = (invnx - 1) * (invnz - 1);
		velModel = (double *)malloc(dims[0] * dims[1] * sizeof(double));
	}
	if (dimension == 3) {
		dims[0] = nz; dims[1] = nx; dims[2] = ny;
		lengthG = (invnx - 1) * (invny - 1) * (invnz - 1);
		velModel = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double));
	}

	velModelLoad(dimension, velModel_file, velModel);

	sourceLoad(dimension, source_file, SourcePoints);

	receiverLoad(dimension, receiver_file, ReceiverPoints);

	rayPointsOutPut.open("./Out/raypoints.txt");

	// ----------------------- Time Map Calculations & Ray Tracing--------------------
	if (dimension == 2) {
		T = (double *)malloc(dims[0] * dims[1] * sizeof(double));
		for (int i = 0; i < SourcePoints.size(); i++) {
			izsrcn = round((SourcePoints.at(i)[1] - zmin) / dz);
			ixsrcn = round((SourcePoints.at(i)[2] - xmin) / dx);
			cout << "Time Map Computing: Source " << i + 1 << endl;
			msfm2dCpp(T, velModel, izsrcn, ixsrcn, dims, usesecond, usecross);

			for (int j = 0; j < ReceiverPoints[int(SourcePoints.at(i)[0])].size(); j++) {
				raytracing(dimension, nx, 1, nz, dx, dy, dz, xmin, 0, zmin, \
				           invnx, 1, invnz, invdx, invdy, invdz, invxmin, invymin, invzmin, \
				           SourcePoints.at(i)[2], 0, SourcePoints.at(i)[1], \
				           ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2), \
				           0, \
				           ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1), \
				           T, rayPoints);
				cout << "Raytracing: " << "Source " << i + 1 << ", " << "Receiver " << j + 1 << endl;

				rayPointsOutPut << int(SourcePoints.at(i)[0]) << " " << ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0) \
				                << " " << rayPoints.size() << endl;
				for (int k = 0; k < rayPoints.size(); k++) {
					rayPointsOutPut << rayPoints.at(k).at(0) << "  " << rayPoints.at(k).at(2) << endl;

					// rayPointsOutPut << int(SourcePoints.at(i)[0])  << " "  \
					//                 << ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0) << " " \
					//                 << rayPoints.at(k).at(0) << " " << rayPoints.at(k).at(2) << endl;
				}
				rayPoints.clear();
			}
		}
	}

	if (dimension == 3) {

		for (int i = 0; i < SourcePoints.size(); i++) {
			for (int j = 0; j < ReceiverPoints[int(SourcePoints.at(i)[0])].size(); j++) {
				localDimension(SourcePoints.at(i)[2], SourcePoints.at(i)[3], SourcePoints.at(i)[1], \
				               ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2), \
				               ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(3), \
				               ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1), \
				               dxyz_loc, dim_x, dim_y, dim_z);
				T = (double *)malloc(dim_x * dim_y * dim_z * sizeof(double) );
				locVelModel = (double *)malloc(dim_x * dim_y * dim_z * sizeof(double) );
				// coordinate transformation
				coordinateTransform(nx, ny, nz, xmin, ymin, zmin, dx, SourcePoints.at(i)[2], SourcePoints.at(i)[3], SourcePoints.at(i)[1], \
				                    ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2), \
				                    ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(3), \
				                    ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1), \
				                    dim_x, dim_y, dim_z, dxyz_loc, velModel, loc_xsrc, loc_ysrc, loc_zsrc, loc_xrec, loc_yrec, loc_zrec, \
				                    ixsrcn, iysrcn, izsrcn, xmin_loc, ymin_loc, zmin_loc, \
				                    theta, locVelModel);
				dims[0] = dim_z;
				dims[1] = dim_x;
				dims[2] = dim_y;

				cout << dim_x << dim_y << dim_z << endl;

				cout << "Time Map Computing: Source " << i + 1 << endl;

				msfm3dCpp(T, locVelModel, izsrcn, ixsrcn, iysrcn, dims, usesecond, usecross);

				raytracing(dimension, dim_x, dim_y, dim_z, dxyz_loc, dxyz_loc, dxyz_loc, xmin_loc, ymin_loc, zmin_loc, \
				           invnx, invny, invnz, invdx, invdy, invdz, invxmin, invymin, invzmin, \
				           loc_xsrc, loc_ysrc, loc_zsrc, loc_xrec, loc_yrec, loc_zrec, T, rayPoints);

				coordinateBack(rayPoints, theta, SourcePoints.at(i)[2], SourcePoints.at(i)[3], SourcePoints.at(i)[1]);

				cout << "Raytracing: " << "Source " << i + 1 << ", " << "Receiver " << j + 1 << endl;

				// rayPointsOutPut << int(SourcePoints.at(i)[0]) << " " << ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0) \
				//                 << " " << rayPoints.size() << endl;
				for (int k = 0; k < rayPoints.size(); k++) {
					// rayPointsOutPut << rayPoints.at(k).at(0) << " " << rayPoints.at(k).at(1) << " " << rayPoints.at(k).at(2) << endl;

					rayPointsOutPut << int(SourcePoints.at(i)[0])  << " "  \
					                << ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0) << " " \
					                << rayPoints.at(k).at(0) << " " << rayPoints.at(k).at(1) << " " << rayPoints.at(k).at(2) << endl;
				}
				rayPoints.clear();
				free(T);
				free(locVelModel);
			}
		}
	}

	// -------------------------------- Clean Up ------------------------------------
	rayPointsOutPut.close();

	return 0;
}
