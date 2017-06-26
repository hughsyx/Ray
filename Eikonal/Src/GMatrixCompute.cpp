// ---------------
//   Dongzhuo Li
//   April, 2015
// ---------------
// 2D temporarily

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "armadillo"

using namespace std;
using namespace arma;

void GMatrixCompute(int dimension, double invnx, double invnz, double invdx, \
                    double invdz, double invxmin, double invzmin, double dx, mat &Gmatrix) {

	string sourceIndex, receiverIndex, nPoints_string, x , y , z;
	int nPoints;
	double dist = dx / 5.0; // the constant distance between two adjacent ray points
	ifstream rayfile;
	rayfile.open("./Out/raypoints.txt");
	string line;
	stringstream ss;
	mat rayPoints, rayPointsUpLeft, rayPointsUpRight, rayPointsDownLeft, rayPointsDownRight;

	if (dimension == 2) {
		// for every ray
		while (getline(rayfile, line)) {
			ss.clear();
			ss << line;
			ss >> sourceIndex >> receiverIndex >> nPoints_string;
			nPoints = stoi(nPoints_string);
			rayPoints.set_size(nPoints, 2);
			rayPointsUpLeft.set_size(nPoints, 2);
			rayPointsUpRight.set_size(nPoints, 2);
			rayPointsDownLeft.set_size(nPoints, 2);
			rayPointsDownRight.set_size(nPoints, 2);
			// read in every point on this ray
			cout << nPoints << endl;
			for (int i = 0; i < nPoints; i++) {
				getline(rayfile, line);
				ss.clear();
				ss << line;
				ss >> z >> x;
				rayPoints(i, 0) = stod(z);
				rayPoints(i, 1) = stod(x);
			}
			// transform the ray points coordinates from real to grid
			rayPoints.col(0) = (rayPoints.col(0) - invzmin) / invdz;
			rayPoints.col(1) = (rayPoints.col(1) - invxmin) / invdx;
			rayPointsUpLeft.col(0) = floor(rayPoints.col(0));
			rayPointsUpLeft.col(1) = floor(rayPoints.col(1));
			rayPointsUpRight.col(0) = floor(rayPoints.col(0));
			rayPointsUpRight.col(1) = ceil(rayPoints.col(1));
			rayPointsDownLeft.col(0) = ceil(rayPoints.col(0));
			rayPointsDownLeft.col(1) = floor(rayPoints.col(1));
			rayPointsDownRight.col(0) = ceil(rayPoints.col(0));
			rayPointsDownRight.col(1) = ceil(rayPoints.col(1));
		}
		rayPoints.print("rayPoints");
		rayPointsUpLeft.print("rayPointsUpLeft");
	}
}