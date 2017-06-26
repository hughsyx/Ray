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
	// parameters in Par_file
	int dimension, nx, ny, nz, invnx, invny, invnz, lengthG, \
	ixsrcn, iysrcn, izsrcn, ixrecn, iyrecn, izrecn, iglob, i1, j1, k1;
	double dx, dy, dz, xmin, ymin, zmin, invdx, invdy, invdz, invxmin, invymin, invzmin;
	string velModel_file, source_file, receiver_file;
	double *velModel, *T;
	bool usesecond = true;
	bool usecross = false;
	ofstream timeMap, rayPointsOutPut;
	vector <double *> SourcePoints;
	vector <double *> AllTravelTimeMaps;
	map < int, vector< vector <double> > > ReceiverPoints;
	vector< vector<double> > rayPoints;
	vector<double> vtemp;
	map<int, map < int, vector< vector <double> > > > AllOutPuts;


	// ------------------------ Read Parameters & Load Model ----------------------
	read_parafile(dimension, nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, invnx, invny, invnz, \
	              invdx, invdy, invdz, invxmin, invymin, invzmin, velModel_file, source_file, receiver_file);

	if (dimension == 2) {
		dims[0] = nz; dims[1] = nx;
		lengthG = (invnx - 1) * (invnz - 1);
		velModel = (double *)malloc(dims[0] * dims[1] * sizeof(double));
		T        = (double *)malloc(dims[0] * dims[1] * sizeof(double));
	}
	if (dimension == 3) {
		dims[0] = nz; dims[1] = nx; dims[2] = ny;
		lengthG = (invnx - 1) * (invny - 1) * (invnz - 1);
		velModel = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double));
		T        = (double *)malloc(dims[0] * dims[1] * dims[2] * sizeof(double) );
	}

	velModelLoad(dimension, velModel_file, velModel);

	sourceLoad(dimension, source_file, SourcePoints);

	receiverLoad(dimension, receiver_file, ReceiverPoints);

	// --------------------------- Time Map Calculations ----------------------------
	// output time map
	//timeMap.open("./Out/timeMap.txt");
	rayPointsOutPut.open("./Out/raypoints.txt");
	for (int i = 0; i < SourcePoints.size(); i++) {
		if (dimension == 2) {
			izsrcn = round((SourcePoints.at(i)[1] - zmin) / dz);
			ixsrcn = round((SourcePoints.at(i)[2] - xmin) / dx);
			cout << "Time Map Computing: Source " << i + 1 << endl;
			msfm2dCpp(T, velModel, izsrcn, ixsrcn, dims, usesecond, usecross);
			// ------------------------------- Ray Tracing -------------------------------
			for (int j = 0; j < ReceiverPoints[int(SourcePoints.at(i)[0])].size(); j++) {
				izrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1) - zmin) / dz);
				ixrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2) - xmin) / dx);
				//iglob  = dims[0]*ixsrcn+izsrcn;
				//cout << "iglob = " << iglob << " T " << dx*T[iglob] << endl;
				iglob  = dims[0] * ixrecn + izrecn;
				//cout << "iglob = " << iglob << " T " << dx*T[iglob] << endl;
				//cout << "izrecn = " << izrecn << "ixrecn = " << ixrecn << endl;
				raytracing(dimension, nx, 1, nz, dx, dy, dz, xmin, 0, zmin, \
				           invnx, 1, invnz, invdx, invdy, invdz, invxmin, invymin, invzmin, \
				           SourcePoints.at(i)[2], 0, SourcePoints.at(i)[1], \
				           ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2), \
				           0, \
				           ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1), \
				           T, rayPoints);
				cout << "Raytracing: " << "Source " << i + 1 << ", " << "Receiver " << j + 1 << endl;

				rayPointsOutPut << int(SourcePoints.at(i)[0]) << " " << ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0) \
				                << " " << rayPoints.size() << " " << dx *T[iglob] << endl;
				for (int k = 0; k < rayPoints.size(); k++) {
					rayPointsOutPut << rayPoints.at(k).at(0) << "  " << rayPoints.at(k).at(2) << endl;

					// rayPointsOutPut << int(SourcePoints.at(i)[0])  << " "  \
					// << ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0) << " " \
					// << rayPoints.at(k).at(0) << " " << rayPoints.at(k).at(2) << endl;
				}
				rayPoints.clear();
			}
		}
		if (dimension == 3) {
			izsrcn = round((SourcePoints.at(i)[1] - zmin) / dz);
			ixsrcn = round((SourcePoints.at(i)[2] - xmin) / dx);
			iysrcn = round((SourcePoints.at(i)[3] - ymin) / dy);
			cout << "izsrcn = " << izsrcn << "ixsrcn = " << ixsrcn << "iysrcn" << iysrcn << endl;

			cout << "Time Map Computing: Source " << i + 1 << endl;
			msfm3dCpp(T, velModel, izsrcn, ixsrcn, iysrcn, dims, usesecond, usecross);
			// output time map
			//for (int i = 0; i < dims[0]*dims[1]*dims[2]; i++) {
			//j1 = int(i/(dims[0]*dims[1]));
			//i1 = int((i-j1*dims[0]*dims[1])/dims[0]);
			//k1 = i-j1*dims[0]*dims[1]-i1*dims[0];
			//timeMap << xmin+i1*dx <<" "<<  ymin+j1*dy <<" "<<  zmin+k1*dz <<" "<< T[i] << endl;
			//}
			// ------------------------------- Ray Tracing ----------------------------------
			for (int j = 0; j < ReceiverPoints[int(SourcePoints.at(i)[0])].size(); j++) {
				izrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1) - zmin) / dz);
				ixrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2) - xmin) / dx);
				iyrecn = round((ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(3) - ymin) / dy);
				//iglob  = dims[0]*dims[1]*iysrcn+dims[0]*ixsrcn+izsrcn;
				//cout << "iglob = " << iglob << " T " << dx*T[iglob] << endl;
				iglob  = dims[0] * dims[1] * iyrecn + dims[0] * ixrecn + izrecn;
				//cout << "iglob = " << iglob << " T " << dx*T[iglob] << endl;
				//cout << "izrecn = " << izrecn << "ixrecn = " << ixrecn << "iysrcn = " << iyrecn << endl;
				raytracing(dimension, nx, ny, nz, dx, dy, dz, xmin, ymin, zmin, \
				           invnx, invny, invnz, invdx, invdy, invdz, invxmin, invymin, invzmin, \
				           SourcePoints.at(i)[2], SourcePoints.at(i)[3], SourcePoints.at(i)[1], \
				           ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(2), \
				           ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(3), \
				           ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(1), \
				           T, rayPoints);
				cout << "Raytracing: " << "Source " << i + 1 << ", " << "Receiver " << j + 1 << endl;
				vtemp.push_back(int(SourcePoints.at(i)[0]));
				vtemp.push_back(ReceiverPoints[int(SourcePoints.at(i)[0])].at(j).at(0));
				vtemp.push_back(rayPoints.size());
				vtemp.push_back(dx * T[iglob]);
				AllOutPuts[i][j].push_back(vtemp);
				vtemp.clear();
				for (int k = 0; k < rayPoints.size(); k++) {
					vtemp.push_back(rayPoints.at(k).at(0));
					vtemp.push_back(rayPoints.at(k).at(1));
					vtemp.push_back(rayPoints.at(k).at(2));
					AllOutPuts[i][j].push_back(vtemp);
					vtemp.clear();
				}
				rayPoints.clear();
			}
		}
	}
	// output to ascii file
	for (int i = 0; i < SourcePoints.size(); i++) {
		for (int j = 0; j < ReceiverPoints[int(SourcePoints.at(i)[0])].size(); j++) {
			rayPointsOutPut << int(AllOutPuts[i][j].at(0).at(0)) << " " << AllOutPuts[i][j].at(0).at(1) \
			                << " " << AllOutPuts[i][j].at(0).at(2) << " " << AllOutPuts[i][j].at(0).at(3) << endl;
			for (int k = 1; k < AllOutPuts[i][j].size(); k++ ) {
				rayPointsOutPut << AllOutPuts[i][j].at(k).at(0) << " " << AllOutPuts[i][j].at(k).at(1) \
				                << " " << AllOutPuts[i][j].at(k).at(2) << endl;
			}
		}
	}


	// output time map
	// timeMap.close();
	// -------------------------------- Clean Up ------------------------------------
	rayPointsOutPut.close();
	free(velModel);
	for (int i = 0; i < SourcePoints.size(); i++) {
		delete SourcePoints.at(i);
	}

	return 0;
}
