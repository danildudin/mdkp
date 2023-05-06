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

enum XType {IN, OUT, CORE};

class LPSolution {
public:
	struct Item {
		int id;
		double x, rc;
	};

	double cost;
	vector<Item> items;
};

class Mkp {
public:
	Mkp(const Problem& p): problem(p), cost(0), b(problem.b) {
		weights.resize(problem.m);
		for (int i = 0; i < problem.n; i++) {
			core.insert(i);
		}
	}

	void set(int id, XType type) {
		if (type == IN) {
			if (n1.find(id) != n1.end()) return;

			cost += problem.c[id];
			for (int i = 0; i < problem.m; i++) {
				weights[i] += problem.a[i][id];
			}

			n1.insert(id);
			n0.erase(id);
			core.erase(id);
		} else if (type == OUT) {
			if (n0.find(id) != n0.end()) return;

			if (n1.find(id) != n1.end()) {
				cost -= problem.c[id];
				for (int i = 0; i < problem.m; i++) {
					weights[i] -= problem.a[i][id];
				}
			}

			n0.insert(id);
			n1.erase(id);
			core.erase(id);
		} else {
			if (core.find(id) != core.end()) return;

			if (n1.find(id) != n1.end()) {
				cost -= problem.c[id];
				for (int i = 0; i < problem.m; i++) {
					weights[i] -= problem.a[i][id];
				}
			}

			n0.erase(id);
			n1.erase(id);
			core.insert(id);
		}
	}

	bool is_n0(int id) {
		return n0.find(id) != n0.end();
	}

	bool is_n1(int id) {
		return n1.find(id) != n1.end();
	}

	bool is_core(int id) {
		return core.find(id) != core.end();
	}

	bool is_feasible() const {
		for (int i = 0; i < problem.m; i++) {
			if (weights[i] > b[i]) {
				return false;
			}
		}

		return true;
	}

	Mkp subproblem() const {
		Mkp res(problem);
		res.core = core;
		for (int i = 0; i < b.size(); i++) {
			res.b[i] = b[i] - weights[i];
		}

		return res;
	}

	LPSolution lp_relaxation() const {
		int m = problem.m;
		int n = core.size();
		vector<int> ids(core.begin(), core.end());

		int* ia = new int[1 + n * m];
		int* ja = new int[1 + n * m];
		double* ar = new double[1 + n * m];
	
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				ia[i * n + j + 1] = i + 1;
				ja[i * n + j + 1] = j + 1;
				ar[i * n + j + 1] = problem.a[i][ids[j]];
			}
		}

		glp_prob *lp;
		lp = glp_create_prob();
		glp_set_obj_dir(lp, GLP_MAX);

		glp_add_rows(lp, m);
		for (int i = 0; i < m; i++) {
			glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, b[i] - weights[i]);
		}

		glp_add_cols(lp, n);
		for (int j = 0; j < n; j++) {
			glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
			glp_set_obj_coef(lp, j + 1, problem.c[ids[j]]);
		}

		glp_load_matrix(lp, n * m, ia, ja, ar);
		glp_simplex(lp, NULL);

		LPSolution res;
		res.cost = glp_get_obj_val(lp);
		for (int j = 0; j < n; j++) {
			res.items.emplace_back(ids[j], glp_get_col_prim(lp, j + 1), glp_get_col_dual(lp, j + 1));
		}
		sort(res.items.begin(), res.items.end(), [&](const LPSolution::Item &a, &b) { return abs(a.rc) < abs(b.rc) });

		glp_delete_prob(lp);
		return res;
	}

public:
	const Problem& problem;
	int cost;
	set<int> n0, n1, core;
	vector<int> weights;
	vector<int> b;
};

class Problem {
public:
	void init() {
		cin >> n >> m;
		
		a.resize(m, vector<int>(n));
		b.resize(m);
		c.resize(n);
		for (int j = 0; j < n; j++) {
			cin >> c[i];
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				cin >> a[i][j];
			}
			cin >> b[i];
		}
	}

public:
	int n, m;
	vector<vector<int>> a;
	vector<int> c, b;
}

int find_k(const Mkp &mkp, int lb, bool need_min) {
	int m = mkp.problem.m;
	int n = mkp.core.size();
	vector<int> ids(mkp.core.begin(), mkp.core.end());

	int* ia = new int[1 + n * (m + 1)];
	int* ja = new int[1 + n * (m + 1)];
	double* ar = new double[1 + n * (m + 1)];

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			ia[i * n + j + 1] = i + 1;
			ja[i * n + j + 1] = j + 1;
			ar[i * n + j + 1] = mkp.problem.a[i][ids[j]];
		}
	}
	for (int j = 0; j < n; j++) {
		ia[m * n + j + 1] = m + 1;
		ja[m * n + j + 1] = j + 1;
		ar[m * n + j + 1] = problem.c[ids[j]];
	}

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, need_min ? GLP_MIN : GLP_MAX);

	glp_add_rows(lp, m + 1);
	for (int i = 0; i < m; i++) {
		glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, mkp.b[i]);
	}
	glp_set_row_bnds(lp, m, GLP_LO, lb + 1, 0.0);

	glp_add_cols(lp, n);
	for (int j = 0; j < n; j++) {
		glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j + 1, 1);
	}

	glp_load_matrix(lp, n * (m + 1), ia, ja, ar);
	glp_simplex(lp, NULL);

	int res = need_min ? ceil(glp_get_obj_val(lp)) : trunk(glp_get_obj_val(lp));
	glp_delete_prob(lp);

	return res;
}

void search_tree(Mkp &cur, int k, int &lb, Mkp &res) {
	if (cur.n1.size() == k) {
		if (cur.cost > res.cost) {
			res = cur;
		}
		return;
	}
	if (cur.core.empty()) {
		return;
	}

	int id = sorted[j];
	cur.set(id, IN);
	search_tree(...)
	cur.set(id, OUT);
	search_tree(...);
	cur.set(id, CORE);
}

Mkp solve_restricted_core_problem(Mkp mkp, int lb) {
	vector<int> sorted(mkp.core.begin(), mkp.core.end());
	sort(sorted.begin(), sorted.end(), [&](int a, b) { return mkp.problem.c[a] > mkp.problem.c[b] });

	Mkp res;
	int k_min = find_k(mkp, lb, true);
	int k_max = find_k(mkp, lb, false);
	for (int k = k_min; k <= k_max; k++) {
		search_tree(mkp, lb, sorted, res);
	}

	return res;
}

pair<Mkp, int> variable_fixing(Mkp cur, int lb) {
	auto lp_res = cur.lp_relaxation();
	int ub = trunk(lp_res.cost);
	int cur_cost = 0;

	while (!lp_res.empty()) {
		const auto& item = lp_res.items.back();
		if (abs(item.rc) < ub - lb) break;

		if (item.x == 1.0) {
			cur.set(item.id, IN);
			cur_cost += cur.problem.c[item.id];
		} else {
			cur.set(item.id, OUT);
		}
		lp_res.pop_back(items);
	}

	while (cur.core.size() > 30) {
		const auto& item = lp_res.items.back();

		auto tmp = cur;
		int tmp_cost = cur_cost;

		if (item.rc > 0) {
			tmp.set(item.id, OUT);
		} else {
			tmp.set(item.id, IN);
			tmp_cost += tmp.problem.c[item.id];
		}
		for (const auto &t : lp_res.items) {
			if (abs(t.rc) > ub - lb - abs(item.rc)) {
				if (t.x == 1.0) {
					tmp.set(t, IN);
					tmp_cost += tmp.problem.c[t.id];
				} else {
					tmp.set(t, OUT);
				}	
			}
		}

		int z;
		if (tmp.core.size() > 30) {
			tie(tmp, z) = variable_fixing(tmp, lb - tmp_cost);
		} else {
			tmp = solve_restricted_core_problem(tmp.subproblem(), lb - tmp_cost);
			z = tmp.cost;
		}
		if (z + tmp_cost > lb) {
			lb = z + tmp_cost;
		} else {
			if (item.rc > 0) {
				cur.set(item.id, IN);
				cur_cost += cur.problem.c[item.id];
			} else {
				cur.set(item.id, OUT);
			}
			lp_res.pop_back();
		}
	}

	return cur, lb;
}

Mkp solve_optimal(Mkp cur, int lb) {
	int cur_cost = cur.cost;
	int z;
	tie(cur, z) = variable_fixing(cur, lb - cur.cost);
	lb = max(lb, z + cur_cost);

	auto core = solve_restricted_core_problem(cur.subproblem(), lb - cur.cost);
	cur.merge(core);

	cur = solve_restricted_core_problem(cur, lb);
}

Mkp coral(Problem p) {
	Mkp cur(p), res(p);
	lp_res = cur.lp_relaxation();

	int ub = trunk(lp_res.cost);
	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (item.x == 1.0) {
			cur.set(item.id, IN);
		} else {
			cur.set(item.id, OUT);
		}
	}

	auto cur_res = solve_optimal(cur, res.cost);
	if (cur_res.cost > res.cost) {
		res = cur_res;
	}

	for (int j = m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (abs(item.rc) >= ub - res.cost) break;

		if (cur.is_n0(item.id)) {
			cur.set(item.id, OUT);
		} else {
			cur.set(item.id, IN);
		}

		cur_res = solve_optimal(cur, res.cost);
		if (cur_res.cost > res.cost) {
			res = cur_res;
		}

		cur.set(item.id, CORE);
	}

	return res;
}

int main() {

}