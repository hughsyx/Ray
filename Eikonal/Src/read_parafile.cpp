// ---------------
//   Dongzhuo Li
//   March, 2015
// ---------------

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void read_parafile(int &dimension, int &nx, int &ny, int &nz, double &dx, double &dy, double &dz, double &xmin, \
                   double &ymin, double &zmin, int &invnx, int &invny, int &invnz, double &invdx, double &invdy, double &invdz, \
                   double &invxmin, double &invymin, double &invzmin, string &velModel_file, string &source_file, string &receiver_file) {

	string line;
	ifstream parafile;

	parafile.open("./Par/Par_file.txt");

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			dimension = stoi(line);
			// cout << "dimension = " << dimension << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			nx = stoi(line);
			// cout << "nx = " << nx << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			ny = stoi(line);
			// cout << "ny = " << ny << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			nz = stoi(line);
			// cout << "nz = " << nz << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			dx = stod(line);
			// cout << "dx = " << dx << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			dy = stod(line);
			// cout << "dy = " << dy << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			dz = stod(line);
			// cout << "dz = " << dz << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			xmin = stod(line);
			// cout << "xmin = " << xmin << endl;
			break;
		}
	}


	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			ymin = stod(line);
			// cout << "ymin = " << ymin << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			zmin = stod(line);
			// cout << "zmin = " << zmin << endl;
			break;
		}
	}


	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invnx = stoi(line);
			// cout << "invnx = " << invnx << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invny = stoi(line);
			// cout << "invny = " << invny << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invnz = stoi(line);
			// cout << "invnz = " << invnz << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invdx = stod(line);
			// cout << "invdx = " << invdx << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invdy = stod(line);
			// cout << "invdy = " << invdy << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invdz = stod(line);
			// cout << "invdz = " << invdz << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invxmin = stod(line);
			// cout << "invxmin = " << invxmin << endl;
			break;
		}
	}


	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invymin = stod(line);
			// cout << "invymin = " << invymin << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			invzmin = stod(line);
			// cout << "invzmin = " << invzmin << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			velModel_file = line;
			// cout << "velocity model =" << velModel_file << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			source_file = line;
			// cout << "source file = " << source_file << endl;
			break;
		}
	}

	while (!parafile.eof()) {
		getline(parafile, line);
		if (line != "" && line.at(0) != '#' && !isspace(line.at(0))) {
			receiver_file = line;
			// cout << "receiver file = " << receiver_file << endl;
			break;
		}
	}


	parafile.close();
}