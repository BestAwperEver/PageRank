#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <regex>
using namespace std;

typedef map<pair<int, int>, double> edges;

typedef vector<list<int>> links;

std::vector<std::string> split(const string& input, const string& regex = " ") {
	std::regex re(regex);
	std::sregex_token_iterator
		first{ input.begin(), input.end(), re, -1 },
		last;
	return { first, last };
}

double dot_product(vector<double> X, list<int> IncomingLinks, vector<int> LinksCount) {
	double res = 0;
	for each (int j in IncomingLinks) {
		res += X[j] / LinksCount[j];
	}
	return res;
}

vector<double> power_iteration(links IncomingLinks, vector<int> LinksCount, double d) {
	vector<double> X(LinksCount.size());
	for (int i = 0; i < 20; ++i) {
		cout << "Iteration " << i << ": " << endl;
		for (int i = 0; i < X.size(); ++i) {
			X[i] = 1 - d + d * dot_product(X, IncomingLinks[i], LinksCount);
			cout << X[i] << endl;
		}
		cout << endl;
	}
	return X;
}

int main() {
	ifstream f("in.txt");

	int N;
	f >> N;

	links IncomingLinks(N);
	vector<int> LinksCount(N, 0);
	vector<string> DictIntString(N);
	map<string, int> DictStringInt;
	int k = 0;

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

	f.close();

	auto X = power_iteration(IncomingLinks, LinksCount, .85);

	for (int i = 0; i < X.size(); ++i) {
		cout << DictIntString[i] << ": " << X[i] << endl;
	}

	system("pause");

	return 0;
}