#include <iostream>
#include <vector>
#include <algorithm>
#include <glpk.h>

using namespace std;

// bool is_feasible(const vector<pair<int, int>> &p, const vector<std::vector<int>> &w, const vector<int> &c, const vector<int> &res) {
// 	for (int i = 0; i < m; i++) {
// 		int cur = 0;
// 		for (int j = 0; j < res.size(); j++) {
// 			cur += w[res[j].second][i];
// 		}

// 		if (cur > c[i]) {
// 			return false;
// 		}
// 	}

// 	return true;
// }

// void search_tree(const vector<pair<int, int>> &p, const vector<std::vector<int>> &w, const vector<int> &c, pair<int, vector<int>> &cur_res, res, int k, int pos) {
// 	if (pos == p.size()) {
// 		return;
// 	}
// 	if (cur_res.second.size() == k) {
// 		if (cur_res.first > res.first) {
// 			res = cur_res;
// 		}
// 		return;
// 	}

// 	cur_res.first += p[pos].first;
// 	cur_res.second.push_back(pos);
// 	if (is_feasible(p, w, c, cur_res)) {
// 		search_tree(p, w, c, cur_res, res, k, pos + 1);
// 	}
// 	cur_res.first -= p[pos].first;
// 	cur_res.second.pop_back();

// 	search_tree(p, w, c, cur_res, res, k, pos + 1);
// }

// int find_k_max(const vector<pair<int, int>> &p, const vector<std::vector<int>> &w, const vector<int> &c) {
// 	vector<int> sol;
// 	for (int j = p.size() - 1; j >= 0; j--) {
// 		sol.push_back(j);
// 		if (!is_feasible(p, w, c, sol)) {
// 			return p.size() - j - 1;
// 		}
// 	}

// 	return p.size();
// }

// int find_k_min(const vector<pair<int, int>> &p, const vector<std::vector<int>> &w, const vector<int> &c, int lb) {
// 	vector<int> sol;
// 	int z = 0;
// 	for (int j = 0; j < p.size(); j++) {
// 		sol.push_back(j);
// 		z += p[j];
// 		if (!is_feasible(p, w, c, sol)) {
// 			return p.size() + 1;
// 		}
// 		if (z > lb) {
// 			return j + 1;
// 		}
// 	}

// 	return p.size() + 1;
// }

// pair<int, vector<int>> solve_restricted_core(const vector<pair<int, int>> p, const vector<std::vector<int>> &w, const vector<int> &c, int lb) {
// 	sort(p.begin(), p.end(), [&](const pair<int, int> &a, &b) { a.first < b.first });

// 	int k_min = find_k_min(p, w, c, lb);
// 	int k_max = find_k_max(p, w, c, lb);
// 	pair<int, vector<int>> res;
// 	for (int k = k_min; k <= k_max; k++) {
// 		pair<int, vector<int>> cur;
// 		search_tree(p, w, c, cur, res, k, 0);
// 	}
// }

// int main() {
// 	int n, m;
// 	vector<pair<int, int>> p;
// 	vector<std::vector<int>> w;
// 	vector<int> c;

// 	cin >> n >> m;
// 	p.resize(n);
// 	w.resize(m, std::vector<int>(n));
// 	c.resize(m);

// 	for (int i = 0; i < n; i++) {
// 		cin >> p[i].first;
// 		p[i].second = i;
// 	}
// 	for (int i = 0; i < m; i++) {
// 		for (int j = 0; j < n; j++) {
// 			cin >> w[i][j];
// 		}
// 	}
// 	for (int i = 0; i < m; i++) {
// 		cin >> c[i];
// 	}

// 	auto res = solve_restricted_core(p, w, c);

// }

/*
10 2
84 34 31 14 67 65 86 98 50 7
20 12 7 75 93 21 75 67 34 28 190
41 51 24 40 84 70 34 41 49 27 250
*/

bool is_feasible(const vector<vector<int>> &a, const vector<int> &b, const vector<double> &x) {
	for (int j = 0; j < x.size(); j++) {
		if (x[j] < 0 || x[j] > 1) {
			return false;
		}
	}

	for (int i = 0; i < a.size(); i++) {
		double w = 0;
		for (int j = 0; j < a[i].size(); j++) {
			w += x[j] * double(a[i][j]);
		}

		if (w > b[i]) {
			return false;
		}
	}

	return true;
}

void solve_lp_relaxation(const vector<vector<int>> &a, const vector<int> &b, const vector<int> &c) {
	int m = a.size();
	int n = a[0].size();
	
	int* ia = new int[1 + n * m];
	int* ja = new int[1 + n * m];
	double* ar = new double[1 + n * m];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			ia[i * n + j + 1] = i + 1;
			ja[i * n + j + 1] = j + 1;
			ar[i * n + j + 1] = a[i][j];
		}
	}

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, m);
	for (int i = 0; i < m; i++) {
		glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, b[i]);
	}

	glp_add_cols(lp, n);
	for (int j = 0; j < n; j++) {
		glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j + 1, c[j]);
	}

	glp_load_matrix(lp, n * m, ia, ja, ar);
	glp_simplex(lp, NULL);

	vector<double> res(n);

	cout << "solution:" << endl;
	cout << "z: " << glp_get_obj_val(lp) << endl;
	cout << "x, reduced cost" << endl;
	for (int j = 0; j < n; j++) {
		res[j] = glp_get_col_prim(lp, j + 1);
		cout << res[j] << "\t" << glp_get_col_dual(lp, j + 1) << endl;
	}
	cout << "is_feasible: " << is_feasible(a, b, res) << endl;
}

int main() {
	int n, m;
	vector<int> c;
	vector<std::vector<int>> a;
	vector<int> b;

	cin >> n >> m;
	c.resize(n);
	a.resize(m, std::vector<int>(n));
	b.resize(m);

	for (int i = 0; i < c.size(); ++i) {
		cin >> c[i];
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cin >> a[i][j];
		}
		cin >> b[i];
	}

	solve_lp_relaxation(a, b, c);

	return 0;
}
