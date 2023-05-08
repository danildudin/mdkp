#include <iostream>
#include <vector>
#include <algorithm>
#include <glpk.h>
#include <map>
#include <set>
#include <cmath>
using namespace std;

void init(vector<vector<int>> &a, vector<int> &b, vector<int> &c) {
	int n, m;

	cin >> n >> m;
	a.resize(m, vector<int>(n));
	b.resize(m);
	c.resize(n);

	for (int j = 0; j < n; j++) {
		cin >> c[j];
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cin >> a[i][j];
		}
		cin >> b[i];
	}
}

bool is_feasible(const vector<vector<int>> &a, const vector<int> &b, const vector<int> &res) {
	for (int i = 0; i < b.size(); i++) {
		int sum = 0;
		for (auto id : res) {
			sum += a[i][id];
		}

		if (sum > b[i]) return false;
	}

	return true;
}

int get_z(const vector<int> &c, vector<int> &res) {
	int sum = 0;
	for (auto id : res) {
		sum += c[id];
	}

	return sum;
}

void print(const vector<int> &c, vector<int> &res) {
	cout << get_z(c, res) << endl;
	cout << res.size() << endl;
	for (auto id : res) {
		cout << id << " ";
	}
	cout << endl;
}

void solve(const vector<vector<int>> &a, const vector<int> &b, const vector<int> &c, vector<int> &cur, vector<int> &res, int id) {
	if (id == c.size()) {
		if (is_feasible(a, b, cur) && get_z(c, cur) > get_z(c, res)) {
			res = cur;
		}
		return;
	}

	solve(a, b, c, cur, res, id + 1);
	cur.push_back(id);
	solve(a, b, c, cur, res, id + 1);
	cur.pop_back();
}

int main() {
	vector<vector<int>> a;
	vector<int> b, c;
	init(a, b, c);

	vector<int> res, cur;
	solve(a, b, c, cur, res, 0);

	// cout << "feasible:\t" << is_feasible(a, b, res) << endl;
	print(c, res);
}
