#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <future>
#include <thread>
#include <vector>

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <glpk.h>

#include "CMDArgs.h"
#include "CORALSolver.h"
#include "CoreData.h"
#include "utils.h"

namespace {
	const int CORE_SIZE = 10;
	const double EPS = 1e-9;

	int find_k(const Mdkp &mdkp, int lb, bool need_min) {
		int m = mdkp.problem.m;
		int n = mdkp.n_size[Mdkp::CORE];
		std::vector<int> ids = mdkp.get_list(Mdkp::CORE);

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

		return mdkp.n_size[Mdkp::N1] + k;
	}

	bool condition1(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
		int sum = mdkp.cost + mdkp.problem.c[data.sorted[pos]];
		for (int j = pos + 1; j < std::min(int(data.sorted.size()), pos + 1 + data.k - mdkp.n_size[Mdkp::N1]); j++) {
			sum += mdkp.problem.c[data.sorted[j]];
		}
		return !(sum < lb);
	}

	bool condition2(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
		int sum = mdkp.w_sum + mdkp.problem.a_sum[data.sorted[pos]];
		int cnt = 0;

		for (auto id : data.w_sum_sorted) {
			if (cnt >= data.k - mdkp.n_size[Mdkp::N1] - 1) break;
			if (id == data.sorted[pos] || !mdkp.is_core(id)) continue;

			sum += mdkp.problem.a_sum[id];
			cnt++;
		}

		return !(sum > mdkp.b_sum);
	}

	bool condition3(const CoreData &data, Mdkp &mdkp, Mdkp &res, int &lb, int pos) {
		double threshold = data.lp_res.cost - mdkp.cost + data.start_cost;
		for (auto id : mdkp.get_list(Mdkp::N1)) {
			if (data.lp_res.items_map[id] && data.lp_res.items_map[id]->rc < 0) {
				threshold -= abs(data.lp_res.items_map[id]->rc);
			}
		}
		for (auto id : mdkp.get_list(Mdkp::N0)) {
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
				if (cnt >= data.k - mdkp.n_size[Mdkp::N1] - 1) break;
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
		if (mdkp.n_size[Mdkp::N1] >= data.k) {
			if (mdkp.cost > res.cost) {
				res = mdkp;
				lb = std::max(lb, res.cost.load());
			}
			return;
		}
		if (mdkp.n_size[Mdkp::CORE] == 0) return;
		if (mdkp.n_size[Mdkp::N1] + data.sorted.size() - pos < data.k) return;

		if (check_conditions(data, mdkp, res, lb, pos)) {
			mdkp.set_x(cur_id, Mdkp::N1);
			search_tree(data, mdkp, res, lb, pos + 1);
		}
		mdkp.set_x(cur_id, Mdkp::N0);
		search_tree(data, mdkp, res, lb, pos + 1);
		mdkp.set_x(cur_id, Mdkp::CORE);
	}

	long long get_timestamp() {
		return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	}

	Mdkp solve_restricted_core_problem(Mdkp mdkp, int lb, CORALSolver::Metadata &metadata) {
		metadata.solve_restricted_core_problem_cnt++;
		auto start_time = get_timestamp();

		if (!mdkp.is_feasible() || mdkp.n_size[Mdkp::CORE] == 0) return mdkp;

		CoreData data;
		int k_min = 0, k_max = 0;

		auto data_ts = get_timestamp();
		data = CoreData(mdkp);
		metadata.solve_restricted_core_problem_coredata_ts += get_timestamp() - data_ts;

		auto find_k_ts = get_timestamp();
		k_min = find_k(mdkp, lb, true);
		k_max = find_k(mdkp, lb, false);
		metadata.solve_restricted_core_problem_find_k_ts += get_timestamp() - find_k_ts;
		metadata.for_k_distribution[std::max(0, k_max - k_min + 1)]++;

		Mdkp res = mdkp;
		for (data.k = k_min; data.k <= k_max; data.k++) {
			metadata.search_tree_cnt++;
			auto search_tree_ts = get_timestamp();
			search_tree(data, mdkp, res, lb, 0);
			metadata.search_tree_ts += get_timestamp() - search_tree_ts;
		}

		metadata.solve_restricted_core_problem_ts += get_timestamp() - start_time;
		return res;
	}

	bool double_equal(double a, double b) {
		return abs(a - b) < EPS;
	}

	Mdkp variable_fixing(Mdkp mdkp, Mdkp& res, CORALSolver::Metadata &metadata) {
		if (!mdkp.is_feasible()) return mdkp;

		auto lp_res = mdkp.lp_relaxation();
		int ub = trunc(lp_res.cost);
		int start_cost = mdkp.cost;

		while (!lp_res.items.empty()) {
			const auto& item = lp_res.items.back();
			if (!(abs(item.rc) > ub - res.cost + start_cost)) break;

			if (double_equal(item.x, 1.0)) {
				mdkp.set_x(item.id, Mdkp::N1);
			} else {
				mdkp.set_x(item.id, Mdkp::N0);
			}
			lp_res.items.pop_back();
		}

		while (mdkp.n_size[Mdkp::CORE] > cmdargs.core_size) {
			const auto& item = lp_res.items.back();

			auto tmp = mdkp;
			if (item.rc > 0){
				tmp.set_x(item.id, Mdkp::N0);
			} else {
				tmp.set_x(item.id, Mdkp::N1);
			}
			for (int j = lp_res.items.size() - 2; j >= 0; j--) {
				const auto &t = lp_res.items[j];
				if (!(abs(t.rc) > ub - res.cost + start_cost - abs(item.rc))) break;

				if (double_equal(t.x, 1.0)) {
					tmp.set_x(t.id, Mdkp::N1);
				} else {
					tmp.set_x(t.id, Mdkp::N0);
				}	
			}

			if (tmp.n_size[Mdkp::CORE] > cmdargs.core_size) {
				tmp = variable_fixing(std::move(tmp), res, metadata);	
			}
			tmp = solve_restricted_core_problem(std::move(tmp), res.cost, metadata);

			if (!tmp.is_feasible() || !compare_and_set(res, tmp)) {
				if (tmp.is_n1(item.id)) {
					mdkp.set_x(item.id, Mdkp::N0);
				} else {
					mdkp.set_x(item.id, Mdkp::N1);
				}
				lp_res.items.pop_back();
			}
		}

		

		return mdkp;
	}

	void solve_optimal(Mdkp mdkp, Mdkp* res, CORALSolver::Metadata* metadata) {
		metadata->solve_optimal_cnt++;
		auto start_time = get_timestamp();

		metadata->variable_fixing_cnt++;
		mdkp = variable_fixing(std::move(mdkp), *res, *metadata);
		metadata->variable_fixing_ts += get_timestamp() - start_time;

		mdkp = solve_restricted_core_problem(std::move(mdkp), res->cost, *metadata);
		if (mdkp.is_feasible()) {
			compare_and_set(*res, mdkp);
		}

		metadata->solve_optimal_ts += get_timestamp() - start_time;
	}
};

void CORALSolver::solve() {
	Mdkp mdkp(problem);
	const auto lp_res = mdkp.lp_relaxation();

	int ub = trunc(lp_res.cost);
	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (double_equal(item.x, 1.0)) {
			mdkp.set_x(item.id, Mdkp::N1);
		} else {
			mdkp.set_x(item.id, Mdkp::N0);
		}
	}

	boost::asio::thread_pool pool(cmdargs.thread_count);
	metadata.reserve(problem.n);
	metadata.emplace_back();
	if (cmdargs.thread_count > 1) {
		boost::asio::post(pool, std::bind(solve_optimal, mdkp, &result, &metadata.back()));
	} else {
		solve_optimal(mdkp, &result, &metadata.back());
	}
	for (int j = problem.m; j < problem.n; j++) {
		const auto& item = lp_res.items[j];
		if (abs(item.rc) >= ub - result.cost) break;

		if (mdkp.is_n0(item.id)) {
			mdkp.set_x(item.id, Mdkp::N1);
		} else {
			mdkp.set_x(item.id, Mdkp::N0);
		}
		if (mdkp.is_feasible()) {
			if (cmdargs.thread_count > 1 && j + 1 < problem.n) {
				metadata.emplace_back();
				boost::asio::post(pool, std::bind(solve_optimal, mdkp, &result, &metadata.back()));
			} else {
				solve_optimal(mdkp, &result, &metadata.back());
			}
		}
		mdkp.set_x(item.id, Mdkp::CORE);
	}
	pool.join();
}

void CORALSolver::print_solution() {
	std::cout << result.cost << std::endl;
	std::cout << result.n_size[Mdkp::N1] << std::endl;
	for (auto id : result.get_list(Mdkp::N1)) {
		std::cout << id << " ";
	}
	std::cout << std::endl;

	std::cout << metadata << std::endl;
}

std::ostream& operator<<(std::ostream& os, const CORALSolver::Metadata& md) {
	os << "{";
	os << "\"solve_optimal_cnt\": " << md.solve_optimal_cnt;
	os << ", \"solve_optimal_ts\": " << md.solve_optimal_ts;
	os << ", \"variable_fixing_cnt\": " << md.variable_fixing_cnt;
	os << ", \"variable_fixing_ts\": " << md.variable_fixing_ts;
	os << ", \"variable_fixing_ratio\": " << double(md.variable_fixing_ts) / md.solve_optimal_ts;

	os << ", \"solve_restricted_core_problem_cnt\": " << md.solve_restricted_core_problem_cnt;
	os << ", \"solve_restricted_core_problem_ts\": " << md.solve_restricted_core_problem_ts;
	os << ", \"solve_restricted_core_problem_ratio\": " << double(md.solve_restricted_core_problem_ts) / md.solve_optimal_ts;

	os << ", \"search_tree_cnt\": " << md.search_tree_cnt;
	os << ", \"search_tree_ts\": " << md.search_tree_ts;
	os << ", \"search_tree_ratio\": " << double(md.search_tree_ts) / md.solve_optimal_ts;

	int stop = md.for_k_distribution.size() - 1;
	while (stop >= 5 && md.for_k_distribution[stop] == 0) stop--;
	os << ", \"for_k_distribution\": " << std::vector<int>(md.for_k_distribution.begin(), md.for_k_distribution.begin() + stop + 1);

	os << "}";

	return os;
}
