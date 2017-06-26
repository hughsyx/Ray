// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

void sourceLoad(int dimension, string &source_file, vector<double * > &SourcePoints) {
	int i = 0;
	string index, x , y , z;

	ifstream sfile;
	sfile.open(source_file);
	string line;
	stringstream ss;

	if (dimension == 2) {
		while (getline(sfile, line)) {
			ss.clear();
			ss << line;
			ss >> index >> z >> x;
			SourcePoints.push_back(new double[3]);
			SourcePoints.at(i)[0] = stod(index);
			SourcePoints.at(i)[1] = stod(z);
			SourcePoints.at(i)[2] = stod(x);
			i++;
		}
	}

	if (dimension == 3) {
		while (getline(sfile, line)) {
			ss.clear();
			ss << line;
			ss >> index >> z >> x >> y;
			SourcePoints.push_back(new double[4]);
			SourcePoints.at(i)[0] = stod(index);
			SourcePoints.at(i)[1] = stod(z);
			SourcePoints.at(i)[2] = stod(x);
			SourcePoints.at(i)[3] = stod(y);
			i++;
		}
	}

	sfile.close();
}