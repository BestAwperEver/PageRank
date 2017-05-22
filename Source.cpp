#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <regex>
using namespace std;

typedef vector<vector<int>> links;

std::vector<std::string> split(const string& input, const string& regex = " ") {
	std::regex re(regex);
	std::sregex_token_iterator
		first{ input.begin(), input.end(), re, -1 },
		last;
	return { first, last };
}

inline double dot_product(vector<double>& X, const vector<int>& IncomingLinks, const vector<int>& LinksCount) {
	double res = 0;
	for each (int j in IncomingLinks) {
		res += X[j] / LinksCount[j];
	}
	return res;
}

vector<double> power_iteration(const links& IncomingLinks, const vector<int>& LinksCount, double d, double eps) {
	vector<double> X(LinksCount.size());
	bool coveraged = false;
	vector<bool> cov_vector(X.size(), false);
	int iteration = 0;
	while (!coveraged) {
		coveraged = true;
		cout << "Iteration " << iteration << ':' << endl;
		for (int i = 0; i < X.size(); ++i) {
			if (cov_vector[i] == false) {
				double new_Xi = 1 - d + d * dot_product(X, IncomingLinks[i], LinksCount);
				if (abs(new_Xi - X[i]) > eps) coveraged = false;
				else cov_vector[i] = true;
				X[i] = new_Xi;
			}
			if (i % 100000 == 0) cout << i << endl;
		}
		cout << endl;
		++iteration;
	}
	return X;
}

bool cmd_processing(int argc, char *argv[], bool& test, string& filename, double& epsilon, double& dump_factor) {
	if (argc < 3 || argc > 5) return false;

	string arg;

	for (int i = 1; i < argc; ++i) {
		arg = argv[i];
		int ind = arg.find('=');
		if (ind == string::npos) return false;
		if (arg.substr(0, ind) == "filename") {
			filename = arg.substr(ind + 1);
		} else if (arg.substr(0, ind) == "precision") {
			epsilon = atof(arg.substr(ind + 1).c_str());
		} else if (arg.substr(0, ind) == "test") {
			if (strcmp(argv[i] + 5, "false") == 0) {
				test = false;
			} else if (strcmp(argv[i] + 5, "true") == 0) {
				test = true;
			} else return false;
		} else if (arg.substr(0, ind) == "dumping_factor") {
			dump_factor = atof(arg.substr(ind + 1).c_str());
		}
	}

	return true;
}

int main(int argc, char *argv[]) {

	bool test = false;
	string filename;
	double epsilon = 1, dump_factor = .85;

	if (cmd_processing(argc, argv, test, filename, epsilon, dump_factor) == false) {
		cout << "Command line parameters:\n"
			"filename=X precision=X [test={true|false}] [dumping_factor=X]\n"
			"Check readme.txt for further information" << endl;
		system("pause");
		return 0;
	}

	ifstream f(filename);

	if (!f.is_open()) {
		cout << "Can't open file " << filename << endl;
		system("pause");
		return 0;
	}

	int N;
	int min, max;

	if (test == false) {

		f >> min >> max;

		if (min > max) swap(min, max);

		while (f.peek() + 1) {
			int i;
			f >> i;
			if (i < min) min = i;
			else if (i > max) max = i;
		}

		f.clear();
		f.seekg(0, ios::beg);

		N = max - min + 1;
	} else {
		f >> N;
	}

	links IncomingLinks(N);
	vector<int> LinksCount(N, 0);
	vector<string> DictIntString(N);
	map<string, int> DictStringInt;
	int k = 0;

	if (test == false) {
		while (f.peek() + 1) {
			int donor_index, recepient_index;
			f >> donor_index >> recepient_index;
			IncomingLinks[recepient_index - min].push_back(donor_index - min);
			++LinksCount[donor_index - min];
		}
	} else {
		for (string line; getline(f, line); ) {
			int donor_index, recepient_index;
			auto sites = split(line);
			if (sites.size() < 2) continue;

			auto it = DictStringInt.find(sites[0]);
			if (it == DictStringInt.end()) {
				DictStringInt[sites[0]] = k;
				DictIntString[k] = sites[0];
				donor_index = k++;
			} else {
				donor_index = it->second;
			}

			for (int i = 1; i < sites.size(); ++i) {
				auto it = DictStringInt.find(sites[i]);
				if (it == DictStringInt.end()) {
					DictStringInt[sites[i]] = k;
					DictIntString[k] = sites[i];
					recepient_index = k++;
				} else {
					recepient_index = it->second;
				}
				IncomingLinks[recepient_index].push_back(donor_index);
				++LinksCount[donor_index];
			}
		}
	}

	f.close();

	auto X = power_iteration(IncomingLinks, LinksCount, dump_factor, epsilon);

	if (test == false) {
		ofstream o("out.txt");
		for (int i = 0; i < X.size(); ++i) {
			o << i << '\t' << X[i] << endl;
		}
		o.close();
	} else {
		for (int i = 0; i < X.size(); ++i) {
			cout << DictIntString[i] << ": " << X[i] << endl;
		}
	}

	system("pause");

	return 0;
}