#pragma once
#include <iostream>

class CMDArgs {
public:
	enum SolverType {CORAL, LAST = CORAL};

	void parse_args(int argc, char **argv);

private:
	void set_solver_type();

public:
	SolverType solver_type = CORAL;
	int thread_count = 1;
	int core_size = 10;
};

// std::ostream& operator<<(std::ostream& os, const CMDArgs::SolverType solver_type);
// std::ostream& operator<<(std::ostream& os, const CMDArgs& args);
