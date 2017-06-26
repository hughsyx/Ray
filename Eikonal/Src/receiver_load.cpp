// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

using namespace std;

void receiverLoad(int dimension, string &receiver_file, map <int, vector<vector<double> > > &ReceiverPoints) {
	string sourceIndex, receiverIndex, x , y , z;
	vector <double> newReceiverPoint;

	ifstream rfile;
	rfile.open(receiver_file);
	string line;
	stringstream ss;

	if (dimension == 2) {
		while (getline(rfile, line)) {
			ss.clear();
			ss << line;
			ss >> sourceIndex >> receiverIndex >> z >> x;
			newReceiverPoint.push_back(stod(receiverIndex));
			newReceiverPoint.push_back(stod(z));
			newReceiverPoint.push_back(stod(x));
			ReceiverPoints[stoi(sourceIndex)].push_back(newReceiverPoint);
			newReceiverPoint.clear();
		}
	}

	if (dimension == 3) {
		while (getline(rfile, line)) {
			ss.clear();
			ss << line;
			ss >> sourceIndex >> receiverIndex >> z >> x >> y;
			newReceiverPoint.push_back(stod(receiverIndex));
			newReceiverPoint.push_back(stod(z));
			newReceiverPoint.push_back(stod(x));
			newReceiverPoint.push_back(stod(y));
			ReceiverPoints[stoi(sourceIndex)].push_back(newReceiverPoint);
			newReceiverPoint.clear();
		}
	}

	rfile.close();
}