#include <iostream>

#include "DPSolver.h"

namespace {
	std::vector<int> operator-(std::vector<int> ar) {
		for (int i = 0; i < ar.size(); i++) {
			ar[i] = -ar[i];
		}

		return ar;
	}

	std::vector<int> operator+(std::vector<int> l, const std::vector<int>& r) {
		for (int i = 0; i < l.size(); i++) {
			l[i] += r[i];
		}

		return l;
	}
}

bool DPSolver::is_feasible(const std::vector<int> &weights) const {
	for (int i = 0; i < m; i++) {
		if (weights[i] > b[i]) return false;
	}
	return true;
}

DPSolver::DPItem DPSolver::set_previous(DPSolver::DPItem cur, const DPSolver::DPItem &prev) const {
	cur.cost = c[cur.cur_id] + prev.cost;
	cur.previous = &prev;
	return cur;
}

void DPSolver::init() {
	int etalon;
	std::cin >> n >> m >> etalon;

	a.resize(n, std::vector<int>(m));
	b.resize(m);
	c.resize(n);

	for (int i = 0; i < n; i++) {
		std::cin >> c[i];
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cin >> a[j][i];
		}
	}
	for (int i = 0; i < m; i++) {
		std::cin >> b[i];
	}
}

void DPSolver::solve() {
	dp.resize(n);

	if (is_feasible(a[0])) dp[0][a[0]] = DPSolver::DPItem(c[0], 0, NULL);
	for (int j = 1; j < n; j++) {
		if (is_feasible(a[j])) dp[j][a[j]] = DPSolver::DPItem(c[j], j, NULL);

		for (int k = 0; k < j; k++) {
			for (auto const& [key, val] : dp[k]) {
				auto new_key = key + a[j];
				if (is_feasible(new_key)) {
					auto cur = set_previous(DPSolver::DPItem(c[j], j, NULL), val);
					dp[j][new_key] = std::max(dp[j][new_key], cur);
				}
			}
		}
	}
}

void DPSolver::print_solution() {
	const DPSolver::DPItem* cur_val = NULL;
	for (int j = 0; j < n; j++) {
		for (auto const& [key, val] : dp[j]) {
			if (cur_val == NULL || *cur_val < val) {
				cur_val = &val;
			}
		}
	}

	int cost = cur_val->cost;
	std::vector<int> x;
	while (cur_val != NULL) {
		x.push_back(cur_val->cur_id);
		cur_val = cur_val->previous;
	}
	std::reverse(x.begin(), x.end());

	std::cout << cost << std::endl;
	std::cout << x.size() << std::endl;
	for (auto id : x) {
		std::cout << id << " ";
	}
	std::cout << std::endl;
}
