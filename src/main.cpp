#include <iostream>
#include <memory>

#include "CMDArgs.h"
#include "CORALSolver.h"
#include "ISolver.h"
#include "utils.h"

std::shared_ptr<ISolver> get_solver(const CMDArgs& args) {
	switch (args.solver_type) {
	case (CMDArgs::CORAL): return std::make_shared<CORALSolver>();
	default: return std::make_shared<CORALSolver>();
	}
}

int main(int argc, char **argv) {
	CMDArgs args;
	args.parse_args(argc, argv);

	auto solver = get_solver(args);
	solver->init();
	solver->solve();
	solver->print_solution();
}
