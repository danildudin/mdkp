
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
	Problem p;
	cin >> p.n >> p.m;
	p.resize();


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

#include <iostream>
#include <vector>
#include <algorithm>
#include <glpk.h>
#include <map>

using namespace std;

/*
10 2
84 34 31 14 67 65 86 98 50 7
20 12 7 75 93 21 75 67 34 28 190
41 51 24 40 84 70 34 41 49 27 250
*/

class Item {
public:
	Item(id, int c, vector<int> w): id(id), c(c), w(std::move(w)) {}

public:
	int id;
	int c;
	vector<int> w;
};

class Solution {
public:
	Solution(int m): z(0) {
		w.resize(m);
	};

	bool is_feasible(const vector<int>& b) const {
		for (int i = 0; i < w.size(); i++) {
			if (w[i] > b[i]) {
				return false;
			}
		}
		return true;
	}

	void add(const Item& item) {
		if (x.find(item.id) != x.end()) return;

		x.insert(item.id);
		for (int i = 0; i < w.size(); i++) {
			w[i] += item.w[i];
		}
	}

	void remove(const Item& item) {
		if (x.find(item.id) == x.end()) return;

		x.erase(item.id);
		for (int i = 0; i < w.size(); i++) {
			w[i] -= item.w[i];
		}
	}

public:
	int z;
	set<int> x;
	vector<int> w;
}

class LPSolution {
public:
	double z;
	vector<double> x;
	vector<double> rc;
};

class Problem {
public:
	init() {
		cin >> n >> m;

		items.resize(n);
		for (int i = 0; i < n; i++) {
			cin >> items[i].c;
			items[i].id = i;
			items[i].w.resize(m);
		}

		b.resize(m);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cin >> items[j][i];
			}
			cin >> b[i];
		}
	}

	LPSolution solve_lp_relaxation() const {
		int* ia = new int[1 + n * m];
		int* ja = new int[1 + n * m];
		double* ar = new double[1 + n * m];
	
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				ia[i * n + j + 1] = i + 1;
				ja[i * n + j + 1] = j + 1;
				ar[i * n + j + 1] = item[j].w[i];
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

		LPSolution res;
		res.z = glp_get_obj_val(lp);
		res.x.resize(n);
		res.rc.resize(n);
		for (int j = 0; j < n; j++) {
			res.x[j] = glp_get_col_prim(lp, j + 1);
			res.rc[j] = glp_get_col_dual(lp, j + 1);
		}
		glp_delete_prob(lp);

		return res;
	}

public:
	int n, m;
	vector<Item> items;
	vector<int> b;
};

Solution solve(const Problem &p, int lb) {
	p_new, fixed = variable_fixing(p, lb);

	

	return res;
}

Solution coral(const Problem &p) {
	lp_res = p.solve_lp_relaxation();
	int ub = trunk(lp_res.z);

	vector<int> sorted(p.n);
	for (int i = 0; i < sorted.size(); i++) {
		sorted[i] = i;
	}
	sort(sorted.begin(), sorted.end(), [&](int a, b) { return abs(lp_res[a].rc) < abs(lp_res[b].rc) });


	set<int> n0, n1;
	Problem core;
	core.n = 0;
	core.m = p.m;
	core.b = p.b;

	for (int j = 0; j < sorted.size(); j++) {
		if (j < m) {
			core.n++;
			core.items.push_back(p.items[sorted[j]]);
		} else if (lp_res[sorted[j]] == 1.0) {
			n1.insert(sorted[j]);
			for (int i = 0; i < core.m; i++) {
				core.b[i] -= p.items[sorted[j]].w[i];
			}
		} else {
			n0.insert(sorted[j]);
		}
	}

	Solution cur_res;
	for (auto id : n1) {
		cur_res.add(p.items[id]);
	}
	auto core_res = solve(core, 0);

	Solution res = unite(cur_res, core_res);
	for (int j = m; j < sorted.size(); j++) {
		if (abs(lp_res[sorted[j]]) >= ub - res.z) {
			return res;
		}

		if (n1.find(sorted[j]) != n1.end()) {
			cur_res.remove(p.items[sorted[j]]);
			n1.erase(sorted[j]);
			n0.insert(sorted[j]);
			for (int i = 0; i < p.m; i++) {
				core.b[i] += p.items[sorted[j]].w[i];
			}
		} else {
			cur_res.add(p.items[sorted[j]]);
			n1.insert(sorted[j]);
			n0.erase(sorted[j]);
			for (int i = 0; i < p.m; i++) {
				core.b[i] -= p.items[sorted[j]].w[i];
			}
		}

		core_res = solve(core, res.z - cur_res.z);
		auto cur = unite(cur_res, core_res);
		if (res.z < cur.z) {
			res = cur;
		}

		if (n1.find(sorted[j]) != n1.end()) {
			cur_res.remove(p.items[sorted[j]]);
			for (int i = 0; i < p.m; i++) {
				core.b[i] += p.items[sorted[j]].w[i];
			}
		}
		n0.erase(sorted[j]);
		n1.erase(sorted[j]);
		cur_res.remove(p.items[sorted[j]]);
	}

	return res;
}

int main() {
	Problem p;
	p.init();


}