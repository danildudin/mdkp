#include <cmath>

#include <boost/asio.hpp>

#include "BBSolver.h"
#include "CMDArgs.h"

namespace {
	const int INF = 1000000000;

	int get_upper_bound(const Mdkp& mdkp) {
		double ub = INF;
		auto ids = mdkp.get_list(Mdkp::CORE);
		for (int i = 0; i < mdkp.problem.m; i++) {
			std::sort(ids.begin(), ids.end(), [&](int a, int b) {
				return static_cast<double>(mdkp.problem.c[a]) / mdkp.problem.a[i][a] > static_cast<double>(mdkp.problem.c[b]) / mdkp.problem.a[i][b];
			});

			double cur_ub = mdkp.cost;
			int weight = mdkp.weights[i];
			for (auto id : ids) {
				if (weight + mdkp.problem.a[i][id] <= mdkp.b[i]) {
					weight += mdkp.problem.a[i][id];
					cur_ub += mdkp.problem.c[id];
				} else {
					cur_ub += mdkp.problem.c[id] * static_cast<double>(mdkp.b[i] - weight) / mdkp.problem.a[i][id];
					break;
				}
			}
			ub = std::min(ub, std::trunc(cur_ub));
		}

		return ub;
	}

	void search_tree(Mdkp &cur, Mdkp *res, const std::vector<int> &order, int pos) {
		if (!cur.is_feasible() || get_upper_bound(cur) <= res->cost) return;
		if (pos == order.size()) {
			compare_and_set(*res, cur);
			return;
		}

		cur.set_x(order[pos], Mdkp::N1);
		search_tree(cur, res, order, pos + 1);
		cur.set_x(order[pos], Mdkp::N0);
		search_tree(cur, res, order, pos + 1);
		cur.set_x(order[pos], Mdkp::CORE);
	}

	void search_tree_parallell(Mdkp &cur, Mdkp &res, const std::vector<int> &order, int pos, boost::asio::thread_pool &pool, int stop) {
		if (!cur.is_feasible() || get_upper_bound(cur) <= res.cost) return;
		if (pos == order.size()) {
			compare_and_set(res, cur);
			return;
		}

		if (pos + 1 < stop) {
			cur.set_x(order[pos], Mdkp::N1);
			search_tree_parallell(cur, res, order, pos + 1, pool, stop);
			cur.set_x(order[pos], Mdkp::N0);
			search_tree_parallell(cur, res, order, pos + 1, pool, stop);
			cur.set_x(order[pos], Mdkp::CORE);
			return;
		}

		cur.set_x(order[pos], Mdkp::N1);
		boost::asio::post(pool, std::bind(search_tree, cur, &res, order, pos + 1));
		cur.set_x(order[pos], Mdkp::N0);
		boost::asio::post(pool, std::bind(search_tree, cur, &res, order, pos + 1));
		cur.set_x(order[pos], Mdkp::CORE);
	}
}

void BBSolver::solve() {
	Mdkp cur(problem);
	auto order = cur.get_list(Mdkp::CORE);
	std::sort(order.begin(), order.end(), [&](int a, int b) {
		return static_cast<double>(cur.problem.c[a]) / cur.problem.a_sum[a] > static_cast<double>(cur.problem.c[a]) / cur.problem.a_sum[a];
	});

	if (cmdargs.thread_count <= 1) {
		search_tree(cur, &res, order, 0);
		return;
	}

	int stop = std::ceil(std::log2(cmdargs.thread_count)) + 5;
	boost::asio::thread_pool pool(cmdargs.thread_count);
	search_tree_parallell(cur, res, order, 0, pool, stop);
	pool.join();
}
