#pragma once
#include <iostream>

#include "ISolver.h"
#include "Mdkp.h"
#include "Problem.h"

class CORALSolver final : public ISolver {
public:
	struct Metadata {
		long long duration_time = 0;
		long long k_diff_sum = 0;
		long long solve_optimal = 0;
		long long solve_restricted_core_problem = 0;
		long long solve_restricted_core_problem_time_sum = 0;
		long long solve_restricted_core_problem_lp_relaxation_time_sum = 0;
		long long variable_fixing = 0;
		long long lp_relaxation = 0;
		// std::vector<int> k_cnt = std::vector<int>(101);
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
