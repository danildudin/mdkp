#include <iostream>
#include <vector>
#include <algorithm>
#include <glpk.h>
#include <map>
#include <set>
#include <cmath>
#include <sstream>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
using namespace std;

const int CORE_SIZE = 10;
const double EPS = 1e-9;

struct Metadata {
	long long duration_time = 0;
	long long k_diff_sum = 0;
	long long k_cnt = 0;
	long long solve_optimal = 0;
	long long solve_restricted_core_problem = 0;
	long long variable_fixing = 0;
	long long lp_relaxation = 0;
} metadata;

std::ostream& operator<<(std::ostream& os, const Metadata& md) {
	os << "{";
	os << "\"duration_time\": " << md.duration_time;
	os << ", \"solve_optimal\": " << md.solve_optimal;
	os << ", \"solve_restricted_core_problem\": " << md.solve_restricted_core_problem;
	os << ", \"variable_fixing\": " << md.variable_fixing;
	os << ", \"lp_relaxation\": " << md.lp_relaxation;
	os << ", \"k_diff_sum\": " << md.k_diff_sum;
	os << ", \"k_cnt\": " << md.k_cnt;
	os << ", \"k_diff_sum / k_cnt\": " << (md.k_cnt == 0 ? 0 : md.k_diff_sum / md.k_cnt);
	os << "}";

	return os;
}

template<class T>
void debug_cout(string field_name, const T &value) {
	return;
	cout << "\"+++" << field_name << "\": " << value << "," << endl;
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

struct LPSolution {
	struct Item {
		int id;
		double x, rc;
	};

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
    os  << ", \"items_map\": " << res.items_map;
    os << "}";
    return os;
}

struct Problem {
	void init() {
		int real_ans;
		cin >> n >> m >> real_ans;

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
		}
		for (int i = 0; i < m; i++) {
			cin >> b[i];
		}
	}

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

enum XType {N1, N0, CORE, LAST = CORE};

std::ostream& operator<<(std::ostream& os, const XType xtype) {
	switch (xtype) {
	case N1: os << "\"N1\"";
	case N0: os << "\"N0\"";
	case CORE: os << "\"CORE\"";
	default: os << "\"ERR\"";
	}

	return os;
}

class Mdkp {
public:
	Mdkp(Problem& p): problem(p), cost(0), w_sum(0), x(problem.n, CORE), n_size(LAST + 1), b_sum(0), b(problem.b), weights(problem.m) {
		n_size[CORE] = problem.n;
		for (auto val : b) {
			b_sum += val;
		}
	}

	Mdkp(const Mdkp& other):
		problem(other.problem),
		w_sum(other.w_sum),
		b_sum(other.b_sum),
		n_size(other.n_size),
		x(other.x),
		b(other.b),
		weights(other.weights)
	{
		cost.store(other.cost);
	}

	const Mdkp& operator=(const Mdkp& other) {
		problem = other.problem;
		cost.store(other.cost);
		w_sum = other.w_sum;
		b_sum = other.b_sum;
		n_size = other.n_size;
		x = other.x;
		b = other.b;
		weights = other.weights;

		return *this;
	}

	vector<int> get_list(XType type) const {
		vector<int> res;
		for (int j = 0; j < x.size(); j++) {
			if (x[j] == type) res.push_back(j);
		}
		return res;
	}

	void set_x(int id, XType type) {
		if (x[id] == type) return;

		if (type == N1) {
			cost += problem.c[id];
			w_sum = 0;
			for (int i = 0; i < problem.m; i++) {
				weights[i] += problem.a[i][id];
				w_sum += weights[i];
			}
		} else if (x[id] == N1) {
			cost -= problem.c[id];
			w_sum = 0;
			for (int i = 0; i < problem.m; i++) {
				weights[i] -= problem.a[i][id];
				w_sum += weights[i];
			}
		}

		n_size[x[id]]--;
		n_size[type]++;
		x[id] = type;
	}

	bool is_n0(int id) {
		return x[id] == N0;
	}

	bool is_n1(int id) {
		return x[id] == N1;
	}

	bool is_core(int id) {
		return x[id] == CORE;
	}

	bool is_feasible() const {
		for (int i = 0; i < problem.m; i++) {
			if (weights[i] > b[i]) {
				return false;
			}
		}

		return true;
	}

	LPSolution lp_relaxation() const {
		metadata.lp_relaxation++;

		int m = problem.m;
		int n = n_size[CORE];
		vector<int> ids = get_list(CORE);

		int ia[1 + n * m];
		int ja[1 + n * m];
		double ar[1 + n * m];

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
	atomic_int cost;
	int w_sum, b_sum;
	vector<XType> x;
	vector<int> b, weights, n_size;
};

mutex mu;
bool compare_and_set(Mdkp& res, const Mdkp& other) {
	scoped_lock lock(mu);

	if (res.cost >= other.cost) return false;
	res = other;
	return true;
}

std::ostream& operator<<(std::ostream& os, const Mdkp& mdkp) {
    os << "{";
    os << "\"is_feasible()\": " << mdkp.is_feasible();
    os << ", \"cost\": " << mdkp.cost;
    os << ", \"n_size\": " << mdkp.n_size;
    os << ", \"x\": " << mdkp.x;
    os << ", \"w_sum\": " << mdkp.w_sum;
    os << ", \"b_sum\": " << mdkp.b_sum;
    os << ", \"weights\": " << mdkp.weights;
    os << ", \"b\": " << mdkp.b;
    os << ", \"problem\": " << mdkp.problem;
    os << "}";

    return os;
}

struct CoreData {
	CoreData(const Mdkp &mdkp) {
		k = 0;
		start_cost = mdkp.cost;
		lp_res = mdkp.lp_relaxation();

		sorted = w_sum_sorted = mdkp.get_list(CORE);
		sort(sorted.begin(), sorted.end(), [&](int a, int b) { return mdkp.problem.c[a] > mdkp.problem.c[b]; });
		sort(w_sum_sorted.begin(), w_sum_sorted.end(), [&](int a, int b) { return mdkp.problem.a_sum[a] < mdkp.problem.a_sum[b]; });
		for (int i = 0; i < mdkp.problem.m; i++) {
			auto cur = sorted;
			sort(cur.begin(), cur.end(), [&](int a, int b) { return mdkp.problem.a[i][a] < mdkp.problem.a[i][b]; });
			w_sorted.push_back(std::move(cur));
		}
	}

	int k, start_cost;
	vector<int> sorted, w_sum_sorted;
	vector<vector<int>> w_sorted;
	LPSolution lp_res;
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

int find_k(const Mdkp &mdkp, int lb, bool need_min) {
	int m = mdkp.problem.m;
	int n = mdkp.n_size[CORE];
	vector<int> ids = mdkp.get_list(CORE);;

	int ia[1 + n * (m + 1)];
	int ja[1 + n * (m + 1)];
	double ar[1 + n * (m + 1)];

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			ia[i * n + j + 1] = i + 1;
			ja[i * n + j + 1] = j + 1;
			ar[i * n + j + 1] = mdkp.problem.a[i][ids[j]];
		}
	}
	for (int j = 0; j < n; j++) {
		ia[m * n + j + 1] = m + 1;
		ja[m * n + j + 1] = j + 1;
		ar[m * n + j + 1] = mdkp.problem.c[ids[j]];
	}

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, need_min ? GLP_MIN : GLP_MAX);

	glp_add_rows(lp, m + 1);
	for (int i = 0; i < m; i++) {
		glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, mdkp.b[i] - mdkp.weights[i]);
	}
	glp_set_row_bnds(lp, m + 1, GLP_LO, lb - mdkp.cost + 1, 0.0);

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

	return mdkp.n_size[N1] + k;
}

bool condition1(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
	int sum = mdkp.cost + mdkp.problem.c[data.sorted[pos]];
	for (int j = pos + 1; j < min(int(data.sorted.size()), pos + 1 + data.k - mdkp.n_size[N1]); j++) {
		sum += mdkp.problem.c[data.sorted[j]];
	}
	return !(sum < lb);
}

bool condition2(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
	int sum = mdkp.w_sum + mdkp.problem.a_sum[data.sorted[pos]];
	int cnt = 0;

	for (auto id : data.w_sum_sorted) {
		if (cnt >= data.k - mdkp.n_size[N1] - 1) break;
		if (id == data.sorted[pos] || !mdkp.is_core(id)) continue;

		sum += mdkp.problem.a_sum[id];
		cnt++;
	}

	return !(sum > mdkp.b_sum);
}

bool condition3(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
	double threshold = data.lp_res.cost - mdkp.cost + data.start_cost;
	for (auto id : mdkp.get_list(N1)) {
		if (data.lp_res.items_map[id] && data.lp_res.items_map[id]->rc < 0) {
			threshold -= abs(data.lp_res.items_map[id]->rc);
		}
	}
	for (auto id : mdkp.get_list(N0)) {
		if (data.lp_res.items_map[id] && data.lp_res.items_map[id]->rc > 0) {
			threshold -= abs(data.lp_res.items_map[id]->rc);
		}
	}
	return !(abs(data.lp_res.items_map[data.sorted[pos]]->rc) > threshold);
}

bool condition4(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
	for (int i = 0; i < mdkp.problem.m; i++) {
		int sum = mdkp.weights[i] + mdkp.problem.a[i][data.sorted[pos]];
		int cnt = 0;
		for (auto id : data.w_sorted[i]) {
			if (cnt >= data.k - mdkp.n_size[N1] - 1) break;
			if (id == data.sorted[pos] || !mdkp.is_core(id)) continue;

			sum += mdkp.problem.a[i][id];
			cnt++;
		}

		if (sum > mdkp.b[i]) return false;
	}

	return true;
}

bool check_conditions(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
	return condition1(data, mdkp, res, lb, pos)
		&& condition2(data, mdkp, res, lb, pos)
		&& condition3(data, mdkp, res, lb, pos)
		&& condition4(data, mdkp, res, lb, pos);
}

void search_tree(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
	if (!mdkp.is_feasible()) return;

	int cur_id = data.sorted[pos];
	if (mdkp.n_size[N1] >= data.k) {
		if (mdkp.cost > res.cost) {
			res = mdkp;
			lb = max(lb, res.cost.load());
		}
		return;
	}
	if (mdkp.n_size[CORE] == 0) return;
	if (mdkp.n_size[N1] + data.sorted.size() - pos < data.k) return;

	if (check_conditions(data, mdkp, res, lb, pos)) {
		mdkp.set_x(cur_id, N1);
		search_tree(data, mdkp, res, lb, pos + 1);
	}
	mdkp.set_x(cur_id, N0);
	search_tree(data, mdkp, res, lb, pos + 1);
	mdkp.set_x(cur_id, CORE);
}

Mdkp solve_restricted_core_problem(Mdkp mdkp, int lb) {
	metadata.solve_restricted_core_problem++;

	if (!mdkp.is_feasible() || mdkp.n_size[CORE] == 0) return mdkp;

	CoreData data(mdkp);

	Mdkp res = mdkp;
	int k_min = find_k(mdkp, lb, true);
	int k_max = find_k(mdkp, lb, false);
	metadata.k_diff_sum += k_max - k_min > 0 ? k_max - k_min : 0;
	metadata.k_cnt++;

	for (data.k = k_min; data.k <= k_max; data.k++) {
		search_tree(data, mdkp, res, lb, 0);
	}

	return res;
}

bool double_equal(double a, double b) {
	return abs(a - b) < EPS;
}

Mdkp variable_fixing(Mdkp mdkp, Mdkp& res, string depth) {
	metadata.variable_fixing++;

	if (!mdkp.is_feasible()) return mdkp;

	auto lp_res = mdkp.lp_relaxation();
	int ub = trunc(lp_res.cost);
	int start_cost = mdkp.cost;

	while (!lp_res.items.empty()) {
		const auto& item = lp_res.items.back();
		if (!(abs(item.rc) > ub - res.cost + start_cost)) break;

		if (double_equal(item.x, 1.0)) {
			mdkp.set_x(item.id, N1);
		} else {
			mdkp.set_x(item.id, N0);
		}
		lp_res.items.pop_back();
	}

	while (mdkp.n_size[CORE] > CORE_SIZE) {
		const auto& item = lp_res.items.back();

		auto tmp = mdkp;
		if (item.rc > 0){
			tmp.set_x(item.id, N0);
		} else {
			tmp.set_x(item.id, N1);
		}
		for (int j = lp_res.items.size() - 2; j >= 0; j--) {
			const auto &t = lp_res.items[j];
			if (!(abs(t.rc) > ub - res.cost + start_cost - abs(item.rc))) break;

			if (double_equal(t.x, 1.0)) {
				tmp.set_x(t.id, N1);
			} else {
				tmp.set_x(t.id, N0);
			}	
		}

		if (tmp.n_size[CORE] > CORE_SIZE) {
			tmp = variable_fixing(std::move(tmp), res, depth + "\t");	
		}
		tmp = solve_restricted_core_problem(std::move(tmp), res.cost);

		if (!tmp.is_feasible() || !compare_and_set(res, tmp)) {
			if (tmp.is_n1(item.id)) {
				mdkp.set_x(item.id, N0);
			} else {
				mdkp.set_x(item.id, N1);
			}
			lp_res.items.pop_back();
		}
	}

	return mdkp;
}

void solve_optimal(Mdkp mdkp, Mdkp* res) {
	metadata.solve_optimal++;

	mdkp = variable_fixing(std::move(mdkp), *res, "");
	mdkp = solve_restricted_core_problem(std::move(mdkp), res->cost);
	if (mdkp.is_feasible()) {
		compare_and_set(*res, mdkp);
	}
}

Mdkp coral(Problem problem) {
	Mdkp mdkp(problem), res(problem);
	const auto lp_res = mdkp.lp_relaxation();

	int ub = trunc(lp_res.cost);
	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (double_equal(item.x, 1.0)) {
			mdkp.set_x(item.id, N1);
		} else {
			mdkp.set_x(item.id, N0);
		}
	}

	vector<thread> threads;
	threads.emplace_back(solve_optimal, mdkp, &res);
	for (int j = problem.m; j < problem.n; j++) {
		if (threads.size() >= 31) {
			for (auto &thread : threads) {
				thread.join();
			}
			threads.clear();
		}

		const auto& item = lp_res.items[j];
		if (abs(item.rc) >= ub - res.cost) break;

		if (mdkp.is_n0(item.id)) {
			mdkp.set_x(item.id, N1);
		} else {
			mdkp.set_x(item.id, N0);
		}
		if (mdkp.is_feasible()) {
			threads.emplace_back(solve_optimal, mdkp, &res);
		}
		mdkp.set_x(item.id, CORE);
	}
	for (auto &thread : threads) {
		thread.join();
	}

	return res;
}

int main() {
	Problem problem;
	problem.init();

	metadata.duration_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	Mdkp res = coral(problem);
	metadata.duration_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - metadata.duration_time;

	cout << res.cost << endl;
	cout << res.n_size[N1] << endl;
	for (auto id : res.get_list(N1)) {
		cout << id << " ";
	}
	cout << endl;
	cout << "metadata: " << metadata << endl;
}
