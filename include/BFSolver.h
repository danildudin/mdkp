#pragma once
#include <iostream>

#include "ISolver.h"
#include "Mdkp.h"
#include "Problem.h"

class BFSolver final : public ISolver {
public:
	BFSolver(): res(problem) {}

	void init() override {
		problem.init();
	}

	void solve() override;

	void print_solution() override {
		std::cout << res.cost << std::endl;
		std::cout << res.n_size[Mdkp::N1] << std::endl;
		for (auto id : res.get_list(Mdkp::N1)) {
			std::cout << id << " ";
		}
		std::cout << std::endl;
	}

private:
	void search_tree(Mdkp& cur, int id);

private:
	Problem problem;
	Mdkp res;
};
