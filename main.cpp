#include <iostream>
#include <vector>
#include <algorithm>
#include <glpk.h>
#include <map>
#include <set>
#include <cmath>
using namespace std;

/*
10 2
84 34 31 14 67 65 86 98 50 7
20 12 7 75 93 21 75 67 34 28 190
41 51 24 40 84 70 34 41 49 27 250
*/

const int CORE_SIZE = 4;
const double EPS = 1e-7;

enum XType {IN, OUT, CORE};

class LPSolution {
public:
	struct Item {
		int id;
		double x, rc;
	};

	void print() {

	}

public:
	double cost;
	vector<Item> items;
	vector<Item*> items_map;
};

class Problem {
public:
	void init() {
		cin >> n >> m;
		
		a.resize(m, vector<int>(n));
		a_sum.resize(n);
		b.resize(m);
		c.resize(n);
		for (int j = 0; j < n; j++) {
			cin >> c[j];
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cin >> a[i][j];
				a_sum[j] += a[i][j];
			}
			cin >> b[i];
		}
	}

	void print() const {
		cout << "n:\t" << n << endl;
		cout << "m:\t" << m << endl;
		
		cout << "c:\t";
		for (auto val : c) {
			cout << val << "\t";
		}
		cout << endl;

		cout << "b:\t";
		for (auto val : b) {
			cout << val << "\t";
		}
		cout << endl;

		cout << "a_sum:\t";
		for (auto val : a_sum) {
			cout << val << "\t";
		}
		cout << endl;

		cout << "a:" << endl;
		for (const auto &row : a) {
			for (auto val : row) {
				cout << val << "\t";
			}
			cout << endl;
		}
		cout << endl;
	}

public:
	int n, m;
	vector<vector<int>> a;
	vector<int> c, b, a_sum;
};

class Mkp {
public:
	Mkp(Problem& p): problem(p), cost(0), b(problem.b) {
		weights.resize(problem.m);
		for (int i = 0; i < problem.n; i++) {
			core.insert(i);
		}
	}

	Mkp& operator=(const Mkp& other) {
		problem = other.problem;
		cost = other.cost;
		w_sum = other.w_sum;
		b_sum = other.b_sum;
		n0 = other.n0;
		n1 = other.n1;
		core = other.core;
		weights = other.weights;
		b = other.b;

		return *this;
	}

	void set_x(int id, XType type) {
		if (type == IN) {
			if (n1.find(id) != n1.end()) return;

			cost += problem.c[id];
			w_sum = 0;
			for (int i = 0; i < problem.m; i++) {
				weights[i] += problem.a[i][id];
				w_sum += weights[i];
			}

			n1.insert(id);
			n0.erase(id);
			core.erase(id);
		} else if (type == OUT) {
			if (n0.find(id) != n0.end()) return;

			if (n1.find(id) != n1.end()) {
				cost -= problem.c[id];
				w_sum = 0;
				for (int i = 0; i < problem.m; i++) {
					weights[i] -= problem.a[i][id];
					w_sum += weights[i];
				}
			}

			n0.insert(id);
			n1.erase(id);
			core.erase(id);
		} else {
			if (core.find(id) != core.end()) return;

			if (n1.find(id) != n1.end()) {
				cost -= problem.c[id];
				w_sum = 0;
				for (int i = 0; i < problem.m; i++) {
					weights[i] -= problem.a[i][id];
					w_sum += weights[i];
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

	void merge(const Mkp &other) {
		for (auto id : other.n1) {
			set_x(id, IN);
		}
		for (auto id : other.n0) {
			set_x(id, OUT);
		}
		for (auto id : other.core) {
			set_x(id, CORE);
		}
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
		res.items_map.resize(problem.n);
		res.cost = glp_get_obj_val(lp);
		for (int j = 0; j < n; j++) {
			res.items.emplace_back(ids[j], glp_get_col_prim(lp, j + 1), glp_get_col_dual(lp, j + 1));
			res.items_map[ids[j]] = &res.items.back();
		}
		sort(res.items.begin(), res.items.end(), [&](const LPSolution::Item &a, const LPSolution::Item &b) { return abs(a.rc) < abs(b.rc); });

		glp_delete_prob(lp);
		return res;
	}

	void print() const {
		cout << "z:\t" << cost << endl;
		cout << "feasible:\t" << is_feasible() << endl;
		cout << "wid:\t";
		for (int i = 0; i < weights.size(); i++) {
			cout << i << "\t";
		}
		cout << endl;

		cout << "w:\t";
		for (auto w : weights) {
			cout << w << "\t";
		}
		cout << endl;

		cout << "id:\t";
		for (int j = 0; j < problem.n; j++) {
			cout << j << "\t";
		}
		cout << endl;

		cout << "x:\t";
		for (int j = 0; j < problem.n; j++) {
			if (n1.find(j) != n1.end()) {
				cout << "1\t"; 
			} else if (n0.find(j) != n0.end()) {
				cout << "0\t"; 
			} else if (core.find(j) != core.end()) {
				cout << "N\t";
			} else {
				cout << "E";
			}
		}
		cout << endl;
	}

public:
	Problem& problem;
	int cost, w_sum, b_sum;
	set<int> n0, n1, core;
	vector<int> weights;
	vector<int> b;
};

struct CoreData {
	vector<int> sorted, w_sorted;
	LPSolution lp_res;
	int k;
};

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
		ar[m * n + j + 1] = mkp.problem.c[ids[j]];
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

	int res = need_min ? ceil(glp_get_obj_val(lp)) : trunc(glp_get_obj_val(lp));
	glp_delete_prob(lp);

	return res;
}

bool conditon1(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	int sum = 0;
	for (int j = pos + 1; j < min(data.sorted.size(), pos + 1 + data.k - mkp.n1.size()); j++) {
		sum += mkp.problem.c[data.sorted[j]];
	}
	return !(sum < lb);
}

bool conditon2(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	int wsq = 0;
	int wsq_cnt = 0;
	for (auto id : data.w_sorted) {
		if (mkp.core.find(id) != mkp.core.end()) continue;

		wsq += mkp.problem.a_sum[id];
		wsq_cnt++;
		if (wsq_cnt == data.k - mkp.n1.size()) break;
	}
	return !(mkp.w_sum + mkp.problem.a_sum[data.sorted[pos]] + wsq > mkp.b_sum);
}

bool conditon3(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	double threshold = data.lp_res.cost - lb;
	for (auto id : mkp.n1) {
		threshold -= abs(data.lp_res.items_map[id]->rc);
	}
	return !(abs(data.lp_res.items_map[data.sorted[pos]]->rc) > threshold);
}

bool conditon4(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	for (int i = 0; i < mkp.problem.m; i++) {
		int wsq = 0;
		int wsq_cnt = 0;
		for (auto id : data.w_sorted) {
			if (mkp.core.find(id) != mkp.core.end()) continue;

			wsq += mkp.problem.a[i][id];
			wsq_cnt++;
			if (wsq_cnt == data.k - mkp.n1.size()) break;
		}

		if (mkp.weights[i] + mkp.problem.a[i][data.sorted[pos]] + wsq > mkp.b[i]) return false;
	}

	return true;
}

bool check_conditons(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	return conditon1(data, mkp, res, lb, pos)
		&& conditon2(data, mkp, res, lb, pos)
		&& conditon3(data, mkp, res, lb, pos)
		&& conditon4(data, mkp, res, lb, pos);
}

void search_tree(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	if (!mkp.is_feasible()) return;

	int cur_id = data.sorted[pos];
	if (mkp.n1.size() == data.k) {
		if (mkp.cost > res.cost) {
			res = mkp;
			lb = max(lb, res.cost);
		}
		return;
	}
	if (mkp.core.empty()) return;
	if (mkp.n1.size() + data.sorted.size() - pos < data.k) return;

	if (check_conditons(data, mkp, res, lb, pos)) {
		mkp.set_x(cur_id, IN);
		search_tree(data, mkp, res, lb, pos + 1);
	}
	mkp.set_x(cur_id, OUT);
	search_tree(data, mkp, res, lb, pos + 1);
	mkp.set_x(cur_id, CORE);
}

Mkp solve_restricted_core_problem(Mkp mkp, int lb) {
	if (mkp.core.empty()) return mkp;

	CoreData data{
		vector<int>(mkp.core.begin(), mkp.core.end()),
		vector<int>(mkp.core.begin(), mkp.core.end()),
		mkp.lp_relaxation(),
	};
	sort(data.sorted.begin(), data.sorted.end(), [&](int a, int b) { return mkp.problem.c[a] > mkp.problem.c[b]; });
	sort(data.w_sorted.begin(), data.w_sorted.end(), [&](int a, int b) { return mkp.problem.a_sum[a] < mkp.problem.a_sum[b]; });

	Mkp res = mkp;
	int k_min = find_k(mkp, lb, true);
	int k_max = find_k(mkp, lb, false);
	for (data.k = k_min; data.k <= k_max; data.k++) {
		search_tree(data, mkp, res, lb, 0);
	}

	return res;
}

bool double_equal(double a, double b) {
	return abs(a - b) < EPS;
}

pair<Mkp, int> variable_fixing(Mkp mkp, int lb) {
	if (!mkp.is_feasible()) return make_pair<Mkp, int>(Mkp(mkp.problem), 0);

	auto lp_res = mkp.lp_relaxation();
	int ub = trunc(lp_res.cost);
	int fixed_cost = 0;

	while (!lp_res.items.empty()) {
		const auto& item = lp_res.items.back();
		if (abs(item.rc) < ub - lb) break;

		if (double_equal(item.x, 1.0)) {
			mkp.set_x(item.id, IN);
			fixed_cost += mkp.problem.c[item.id];
		} else {
			mkp.set_x(item.id, OUT);
		}
		lp_res.items.pop_back();
	}

	while (mkp.core.size() > CORE_SIZE) {
		const auto& item = lp_res.items.back();

		auto tmp = mkp;
		int tmp_cost = fixed_cost;
		if (item.rc > 0) {
			tmp.set_x(item.id, OUT);
		} else {
			tmp.set_x(item.id, IN);
			tmp_cost += tmp.problem.c[item.id];
		}
		for (const auto &t : lp_res.items) {
			if (abs(t.rc) > ub - lb - abs(item.rc)) {
				if (double_equal(t.x, 1.0)) {
					tmp.set_x(t.id, IN);
					tmp_cost += tmp.problem.c[t.id];
				} else {
					tmp.set_x(t.id, OUT);
				}	
			}
		}

		int z;
		if (tmp.core.size() > CORE_SIZE) {
			tie(tmp, z) = variable_fixing(tmp, lb - tmp_cost);
		} else {
			tmp = solve_restricted_core_problem(tmp.subproblem(), lb - tmp_cost);
			z = tmp.cost;
		}
		if (tmp.is_feasible() && z + tmp_cost > lb) {
			lb = z + tmp_cost;
		} else {
			if (item.rc > 0) {
				mkp.set_x(item.id, IN);
				fixed_cost += mkp.problem.c[item.id];
			} else {
				mkp.set_x(item.id, OUT);
			}
			lp_res.items.pop_back();
		}
	}

	return make_pair<Mkp, int>(std::move(mkp), std::move(lb));
}

Mkp solve_optimal(Mkp mkp, int lb) {
	if (!mkp.is_feasible()) return Mkp(mkp.problem);

	int z = 0;
	int cur_cost = mkp.cost;
	tie(mkp, z) = variable_fixing(mkp, lb - mkp.cost);
	lb = max(lb, z + cur_cost);

	auto core = solve_restricted_core_problem(mkp.subproblem(), lb - mkp.cost);
	if (core.is_feasible()) {
		mkp.merge(core);
	}

	return mkp;
}

Mkp coral(Problem problem) {
	Mkp mkp(problem), res(problem);
	const auto lp_res = mkp.lp_relaxation();

	int ub = trunc(lp_res.cost);
	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (item.x == 1.0) {
			mkp.set_x(item.id, IN);
		} else {
			mkp.set_x(item.id, OUT);
		}
	}

	auto cur_res = solve_optimal(mkp, res.cost);
	if (cur_res.cost > res.cost) {
		res = cur_res;
	}

	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (abs(item.rc) >= ub - res.cost) break;

		if (mkp.is_n0(item.id)) {
			mkp.set_x(item.id, OUT);
		} else {
			mkp.set_x(item.id, IN);
		}

		cur_res = solve_optimal(mkp, res.cost);
		if (cur_res.cost > res.cost) {
			res = cur_res;
		}

		mkp.set_x(item.id, CORE);
	}

	return res;
}

int main() {
	Problem problem;
	problem.init();

	Mkp res = coral(problem);

	problem.print();
	res.print();
}
