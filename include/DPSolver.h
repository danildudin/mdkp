#pragma once
#include <map>
#include <vector>

#include "ISolver.h"

class DPSolver final : public ISolver {
public:
	class DPItem {
	public:
		friend bool operator<(const DPItem& l, const DPItem& r) {
			return l.cost < r.cost || (l.cost == r.cost && l.cur_id < r.cur_id);
		}

	public:
		int cost = 0;
		int cur_id = -1;
		const DPItem* previous = NULL;
	};

	void init() override;
	void solve() override;
	void print_solution() override;

private:
	bool is_feasible(const std::vector<int> &weights) const;
	DPItem set_previous(DPItem cur, const DPItem &prev) const;

private:
	int n, m;
	std::vector<int> c, b;
	std::vector<std::vector<int>> a;
	std::vector<std::map<std::vector<int>, DPItem>> dp;
};
