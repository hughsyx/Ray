// ---------------
//   Dongzhuo Li & Ping Tong
//   May, 2015
// ---------------

#include <vector>
#include <cmath>

using namespace std;

void coordinateBack(vector< vector<double> > &rayPoints, double theta, double xsrc, double ysrc, double zsrc) {
	double x_local, y_local, z_local, x_glb, y_glb, z_glb;
	for (int i = 0; i < rayPoints.size(); i++) {
		x_local = rayPoints.at(i).at(0);
		y_local = rayPoints.at(i).at(1);
		z_local = rayPoints.at(i).at(2);
		x_glb = cos(theta) * x_local - sin(theta) * y_local + xsrc;
		y_glb = sin(theta) * x_local + cos(theta) * y_local + ysrc;
		z_glb = z_local + zsrc;
		rayPoints.at(i).at(0) = x_glb;
		rayPoints.at(i).at(1) = y_glb;
		rayPoints.at(i).at(2) = z_glb;
	}
}