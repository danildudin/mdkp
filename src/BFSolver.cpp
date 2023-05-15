#include <cmath>
#include <mutex>

#include <boost/asio.hpp>

#include "BFSolver.h"
#include "CMDArgs.h"

namespace {
	std::mutex mu;
	void compare_and_set(Mdkp& res, const Mdkp& other) {
		if (res.cost >= other.cost) return;

		std::scoped_lock lock(mu);
		if (res.cost < other.cost) res = other;
	}

	void search_tree(Mdkp &cur, Mdkp *res, int id) {
		if (id == cur.problem.c.size()) {
			if (cur.is_feasible()) compare_and_set(*res, cur);
			return;
		}

		search_tree(cur, res, id + 1);
		cur.set_x(id, Mdkp::N1);
		search_tree(cur, res, id + 1);
		cur.set_x(id, Mdkp::N0);
	}

	void search_tree_parallell(Mdkp &cur, Mdkp &res, int id, boost::asio::thread_pool &pool, int stop) {
		if (id == cur.problem.c.size()) {
			if (cur.is_feasible()) compare_and_set(res, cur);
			return;
		}

		if (id + 1 < stop) {
			search_tree_parallell(cur, res, id + 1, pool, stop);
			cur.set_x(id, Mdkp::N1);
			search_tree_parallell(cur, res, id + 1, pool, stop);
			cur.set_x(id, Mdkp::N0);
			return;
		}

		boost::asio::post(pool, std::bind(search_tree, cur, &res, id + 1));
		cur.set_x(id, Mdkp::N1);
		boost::asio::post(pool, std::bind(search_tree, cur, &res, id + 1));
		cur.set_x(id, Mdkp::N0);
	}
}

void BFSolver::solve() {
	Mdkp cur(problem);
	if (cmdargs.thread_count <= 1) {
		search_tree(cur, &res, 0);
		return;
	}

	int stop = std::ceil(std::log2(cmdargs.thread_count));
	boost::asio::thread_pool pool(cmdargs.thread_count);
	search_tree_parallell(cur, res, 0, pool, stop);
	pool.join();
}
