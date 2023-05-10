#include <iostream>
#include <vector>
#include <algorithm>
#include <glpk.h>
#include <map>
#include <set>
#include <cmath>
#include <sstream>
using namespace std;

/*
10 2
84 34 31 14 67 65 86 98 50 7
20 12 7 75 93 21 75 67 34 28 190
41 51 24 40 84 70 34 41 49 27 250
*/

const int CORE_SIZE = 5;
const double EPS = 1e-7;

template<class T>
void debug_cout(string field_name, const T &value) {
	return;
	cout << "\"+++" << field_name << "\": " << value << "," << endl;
}

enum XType {IN, OUT, CORE};

std::ostream& operator<<(std::ostream& os, const XType xtype) {
	switch (xtype) {
	case IN: os << "\"IN\"";
	case OUT: os << "\"OUT\"";
	case CORE: os << "\"CORE\"";
	default: os << "\"ERR\"";
	}

	return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const vector<T>& arr) {
    os << "[";
	for (int i = 0; i < arr.size(); i++) {
		os << arr[i];
		if (i + 1 < arr.size()) {
			os << ", ";
		}
	}
	os << "]";

    return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const set<T>& s) {
	os << "[";
	for (auto it = s.begin(); it != s.end();) {
		os << *it;
		if (++it != s.end()) {
			os << ", ";
		}
	}
	os << "]";

	return os;
}

template<class K, class V>
std::ostream& operator<<(std::ostream& os, const map<K, V>& m) {
	os << "{";
	for (auto it = m.begin(); it != m.end();) {
		os << "\"" << it.first << "\": " << it.second;
		if (++it != m.end()) {
			os << ", ";
		}
	}
	os << "}";

	return os;
}

class LPSolution {
public:
	struct Item {
		int id;
		double x, rc;
	};

public:
	double cost;
	vector<Item> items;
	vector<Item*> items_map;
};

std::ostream& operator<<(std::ostream& os, const LPSolution::Item* item) {
	if (!item) {
		os << "\"NULL\"";
	} else {
    	os << "{\"ptr\": 1, \"id\": " << item->id << ", \"x\": " << item->x << ", \"rc\": " << item->rc << "}";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const LPSolution::Item& item) {
    os << "{\"id\": " << item.id << ", \"x\": " << item.x << ", \"rc\": " << item.rc << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const LPSolution& res) {
    os << "{\"cost\": " << res.cost << ", \"items\": " << res.items;
    // os  << ", \"items_map\": " << res.items_map;
    os << "}";
    return os;
}

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

public:
	int n, m;
	vector<vector<int>> a;
	vector<int> c, b, a_sum;
};

std::ostream& operator<<(std::ostream& os, const Problem& p) {
    os << "{";
    os << "\"n\": " << p.n << ", \"m\": " << p.m;
    os << ", \"c\": " << p.c;
    os << ", \"b\": " << p.b;
    os << ", \"a_sum\": " << p.a_sum;
    os << ", \"a\": " << p.a;
    os << "}";
    return os;
}

class Mkp {
public:
	Mkp(Problem& p): problem(p), cost(0), w_sum(0), b_sum(0), b(problem.b) {
		for (auto val : b) {
			b_sum += val;
		}

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
		res.cost = 0;
		res.w_sum = 0;
		res.b_sum = 0;
		for (int i = 0; i < b.size(); i++) {
			res.b[i] = b[i] - weights[i];
		}
		for (auto val : res.b) {
			res.b_sum += val;
		}
		res.core = core;

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

		glp_smcp params;
		glp_init_smcp(&params);
		params.msg_lev = GLP_MSG_ERR;

		glp_simplex(lp, &params);

		LPSolution res;
		res.items_map.resize(problem.n);
		res.cost = glp_get_obj_val(lp);
		for (int j = 0; j < n; j++) {
			res.items.emplace_back(ids[j], glp_get_col_prim(lp, j + 1), glp_get_col_dual(lp, j + 1));
		}
		sort(res.items.begin(), res.items.end(), [&](const LPSolution::Item &a, const LPSolution::Item &b) { return abs(a.rc) < abs(b.rc); });
		
		for (int j = 0; j < n; j++) {
			res.items_map[res.items[j].id] = &res.items[j];
		}

		glp_delete_prob(lp);
		return res;
	}

public:
	Problem& problem;
	int cost, w_sum, b_sum;
	set<int> n0, n1, core;
	vector<int> weights;
	vector<int> b;
};

std::ostream& operator<<(std::ostream& os, const Mkp& mkp) {
    os << "{";
    os << "\"is_feasible()\": " << mkp.is_feasible();
    os << ", \"cost\": " << mkp.cost;
    os << ", \"w_sum\": " << mkp.w_sum;
    os << ", \"b_sum\": " << mkp.b_sum;
    os << ", \"n0\": " << mkp.n0;
    os << ", \"n1\": " << mkp.n1;
    os << ", \"core\": " << mkp.core;
    os << ", \"weights\": " << mkp.weights;
    os << ", \"b\": " << mkp.b;
    // os << ", \"problem\": " << mkp.problem;
    os << "}";

    return os;
}

struct CoreData {
	vector<int> sorted, w_sum_sorted;
	vector<vector<int>> w_sorted;
	LPSolution lp_res;
	int k, start_cost;
};

std::ostream& operator<<(std::ostream& os, const CoreData& data) {
    os << "{";
    os << "\"k\": " << data.k;
    os << ", \"start_cost\": " << data.start_cost;
    os << ", \"sorted\": " << data.sorted;
    os << ", \"w_sum_sorted\": " << data.w_sum_sorted;
    os << ", \"lp_res\": " << data.lp_res;
    os << "}";

    return os;
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
		ar[m * n + j + 1] = mkp.problem.c[ids[j]];
	}

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, need_min ? GLP_MIN : GLP_MAX);

	glp_add_rows(lp, m + 1);
	for (int i = 0; i < m; i++) {
		glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, mkp.b[i] - mkp.weights[i]);
	}
	glp_set_row_bnds(lp, m + 1, GLP_LO, lb - mkp.cost + 1, 0.0);

	glp_add_cols(lp, n);
	for (int j = 0; j < n; j++) {
		glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j + 1, 1);
	}

	glp_load_matrix(lp, n * (m + 1), ia, ja, ar);

	glp_smcp params;
	glp_init_smcp(&params);
	params.msg_lev = GLP_MSG_ERR;

	glp_simplex(lp, &params);

	int k = need_min ? ceil(glp_get_obj_val(lp)) : trunc(glp_get_obj_val(lp));
	glp_delete_prob(lp);

	return mkp.n1.size() + k;
}

bool condition1(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	int sum = mkp.cost + mkp.problem.c[data.sorted[pos]];
	for (int j = pos + 1; j < min(data.sorted.size(), pos + 1 + data.k - mkp.n1.size()); j++) {
		sum += mkp.problem.c[data.sorted[j]];
	}
	return !(sum < lb);
}

bool condition2(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	int sum = mkp.w_sum + mkp.problem.a_sum[data.sorted[pos]];
	int cnt = 0;

	for (auto id : data.w_sum_sorted) {
		if (cnt >= data.k - mkp.n1.size() - 1) break;
		if (id == data.sorted[pos] || mkp.core.find(id) == mkp.core.end()) continue;

		sum += mkp.problem.a_sum[id];
		cnt++;
	}

	return !(sum > mkp.b_sum);
}

bool condition3(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	double threshold = data.lp_res.cost - mkp.cost + data.start_cost;
	for (auto id : mkp.n1) {
		if (data.lp_res.items_map[id] && data.lp_res.items_map[id]->rc < 0) {
			threshold -= abs(data.lp_res.items_map[id]->rc);
		}
	}
	for (auto id : mkp.n0) {
		if (data.lp_res.items_map[id] && data.lp_res.items_map[id]->rc > 0) {
			threshold -= abs(data.lp_res.items_map[id]->rc);
		}
	}
	return !(abs(data.lp_res.items_map[data.sorted[pos]]->rc) > threshold);
}

bool condition4(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	for (int i = 0; i < mkp.problem.m; i++) {
		int sum = mkp.weights[i] + mkp.problem.a[i][data.sorted[pos]];
		int cnt = 0;
		for (auto id : data.w_sorted[i]) {
			if (cnt >= data.k - mkp.n1.size() - 1) break;
			if (id == data.sorted[pos] || mkp.core.find(id) == mkp.core.end()) continue;

			sum += mkp.problem.a[i][id];
			cnt++;
		}

		if (sum > mkp.b[i]) return false;
	}

	return true;
}

bool check_conditions(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	return condition1(data, mkp, res, lb, pos)
		&& condition2(data, mkp, res, lb, pos)
		&& condition3(data, mkp, res, lb, pos)
		&& condition4(data, mkp, res, lb, pos);
}

void search_tree(const CoreData &data, Mkp &mkp, Mkp &res, int &lb, int pos) {
	if (!mkp.is_feasible()) return;

	int cur_id = data.sorted[pos];
	if (mkp.n1.size() >= data.k) {
		if (mkp.cost > res.cost) {
			res = mkp;
			lb = max(lb, res.cost);
		}
		return;
	}
	if (mkp.core.empty()) return;
	if (mkp.n1.size() + data.sorted.size() - pos < data.k) return;

	if (check_conditions(data, mkp, res, lb, pos)) {
		mkp.set_x(cur_id, IN);
		search_tree(data, mkp, res, lb, pos + 1);
	}
	mkp.set_x(cur_id, OUT);
	search_tree(data, mkp, res, lb, pos + 1);
	mkp.set_x(cur_id, CORE);
}

Mkp solve_restricted_core_problem(Mkp mkp, int lb) {
	if (!mkp.is_feasible() || mkp.core.empty()) return std::move(mkp);

	CoreData data;
	data.start_cost = mkp.cost;
	data.lp_res = mkp.lp_relaxation();
	data.sorted = data.w_sum_sorted = vector<int>(mkp.core.begin(), mkp.core.end());
	sort(data.sorted.begin(), data.sorted.end(), [&](int a, int b) { return mkp.problem.c[a] > mkp.problem.c[b]; });
	sort(data.w_sum_sorted.begin(), data.w_sum_sorted.end(), [&](int a, int b) { return mkp.problem.a_sum[a] < mkp.problem.a_sum[b]; });
	for (int i = 0; i < mkp.problem.m; i++) {
		auto cur = data.sorted;
		sort(cur.begin(), cur.end(), [&](int a, int b) { return mkp.problem.a[i][a] < mkp.problem.a[i][b]; });
		data.w_sorted.push_back(std::move(cur));
	}

	Mkp res = mkp;
	int k_min = find_k(mkp, lb , true);
	int k_max = find_k(mkp, lb, false);
	for (data.k = k_min; data.k <= k_max; data.k++) {
		search_tree(data, mkp, res, lb, 0);
	}

	return res;
}

bool double_equal(double a, double b) {
	return abs(a - b) < EPS;
}

Mkp variable_fixing(Mkp mkp, Mkp& res, string depth) {
	if (!mkp.is_feasible()) return mkp;

	auto lp_res = mkp.lp_relaxation();
	int ub = trunc(lp_res.cost);
	int start_cost = mkp.cost;

	while (!lp_res.items.empty()) {
		const auto& item = lp_res.items.back();
		if (!(abs(item.rc) > ub - res.cost + start_cost)) break;

		if (double_equal(item.x, 1.0)) {
			mkp.set_x(item.id, IN);
		} else {
			mkp.set_x(item.id, OUT);
		}
		lp_res.items.pop_back();
	}

	while (mkp.core.size() > CORE_SIZE) {
		const auto& item = lp_res.items.back();

		auto tmp = mkp;
		if (item.rc > 0){
			tmp.set_x(item.id, OUT);
		} else {
			tmp.set_x(item.id, IN);
		}
		for (int j = lp_res.items.size() - 2; j >= 0; j--) {
			const auto &t = lp_res.items[j];
			if (abs(t.rc) > ub - res.cost + start_cost - abs(item.rc)) {
				if (double_equal(t.x, 1.0)) {
					tmp.set_x(t.id, IN);
				} else {
					tmp.set_x(t.id, OUT);
				}	
			} else {
				break;
			}
		}

		if (tmp.core.size() > CORE_SIZE) {
			tmp = variable_fixing(std::move(tmp), res, depth + "\t");	
		}
		tmp = solve_restricted_core_problem(tmp, res.cost);

		if (tmp.is_feasible() && tmp.cost > res.cost) {
			res = tmp;
		} else {
			if (tmp.is_n1(item.id)) {
				mkp.set_x(item.id, OUT);
			} else {
				mkp.set_x(item.id, IN);
			}
			lp_res.items.pop_back();
		}
	}

	return mkp;
}

void solve_optimal(Mkp mkp, Mkp& res) {
	mkp = variable_fixing(mkp, res, "");
	mkp = solve_restricted_core_problem(mkp, res.cost);

	if (mkp.is_feasible() && mkp.cost > res.cost) {
		res = mkp;
	}
}

Mkp coral(Problem problem) {
	Mkp mkp(problem), res(problem);
	const auto lp_res = mkp.lp_relaxation();

	int ub = trunc(lp_res.cost);
	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (double_equal(item.x, 1.0)) {
			mkp.set_x(item.id, IN);
		} else {
			mkp.set_x(item.id, OUT);
		}
	}

	solve_optimal(mkp, res);

	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (abs(item.rc) >= ub - res.cost) break;

		if (mkp.is_n0(item.id)) {
			mkp.set_x(item.id, IN);
		} else {
			mkp.set_x(item.id, OUT);
		}

		solve_optimal(mkp, res);
		mkp.set_x(item.id, CORE);
	}

	return res;
}

int main() {
	Problem problem;
	problem.init();
	// cout << "{" << endl;
	debug_cout("problem", problem);

	Mkp res = coral(problem);
	// Mkp mkp(problem), res(problem);
	// debug_cout("res", res);
	// debug_cout("res.subproblem()", res.subproblem());

	// solve_optimal(mkp, res);
	// auto res = solve_restricted_core_problem(Mkp(problem), 0);

	debug_cout("res", res);
	// cout << "}" << endl;


	cout << res.cost << endl;
	cout << res.n1.size() << endl;
	for (auto id : res.n1) {
		cout << id << " ";
	}
	cout << endl;
}
