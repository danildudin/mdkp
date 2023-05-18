#pragma once
#include <iostream>

#include "ISolver.h"
#include "Mdkp.h"
#include "Problem.h"

class CORALSolver final : public ISolver {
public:
	struct Metadata {
		long long solve_optimal_cnt = 0;
		long long solve_optimal_ts = 0;
		long long variable_fixing_cnt = 0;
		long long variable_fixing_ts = 0;
		long long solve_restricted_core_problem_cnt = 0;
		long long solve_restricted_core_problem_ts = 0;
		long long solve_restricted_core_problem_coredata_ts = 0;
		long long solve_restricted_core_problem_find_k_ts = 0;
		long long search_tree_ts = 0;
		long long search_tree_cnt = 0;
		std::vector<int> for_k_distribution = std::vector<int>(101);
	};

	CORALSolver(): result(problem) {}

	void init() override {
		problem.init();
	}
	void solve() override;
	void print_solution() override;

private:
	Problem problem;
	Mdkp result;
	std::vector<Metadata> metadata;
};

std::ostream& operator<<(std::ostream& os, const CORALSolver::Metadata& md);
