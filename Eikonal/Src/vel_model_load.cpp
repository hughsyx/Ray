// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

void velModelLoad(int dimension, string &velModel_file, double *velModel) {
	int i = 0;
	string a, b , c , d;

	ifstream vfile;
	vfile.open(velModel_file);
	string line;
	stringstream ss;

	if (dimension == 2) {
		while (getline(vfile, line)) {
			ss.clear();
			ss << line;
			ss >> a >> b >> c;
			velModel[i] = stod(c);
			i++;
		}
	}

	if (dimension == 3) {
		while (getline(vfile, line)) {
			ss.clear();
			ss << line;
			ss >> a >> b >> c >> d;
			velModel[i] = stod(d);
			i++;
		}
	}

	vfile.close();
}
