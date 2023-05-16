#include <getopt.h>

#include "CMDArgs.h"
#include "utils.h"

namespace {
    void set_int(int &var) {
    	auto value = atoi(optarg);
    	if (value) var = value;
    }
}

CMDArgs cmdargs;

void CMDArgs::set_solver_type() {
	std::string name = optarg;
	if (name == "coral") {
		solver_type = CORAL;
	} else if (name == "dp") {
		solver_type = DP;
	} else if (name == "bf") {
		solver_type = BF;
	} else if (name == "bb") {
		solver_type = BB;
	}
}

void CMDArgs::parse_args(int argc, char **argv) {
	struct option long_options[] = {
        {"solver_type", required_argument, 0, 0},
        {"thread_count", required_argument, 0, 0},
        {"core_size", required_argument, 0, 0},
        {0, 0, 0, 0}
    };

	int c = 0, option_index = 0;
    while (true) {
        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;
        if (c != 0) continue;

        std::string name = long_options[option_index].name;
        if (name == "solver_type") {
            set_solver_type();
        } else if (name == "thread_count") {
            set_int(thread_count);
        } else if (name == "core_size") {
            set_int(core_size);
        }
    }
}

std::ostream& operator<<(std::ostream& os, const CMDArgs::SolverType solver_type) {
	switch (solver_type) {
	case CMDArgs::CORAL: os << "\"CORAL\"";
	default: os << "\"ERR\"";
	}

	return os;
}

std::ostream& operator<<(std::ostream& os, const CMDArgs& args) {
	os << "{";
	os << "\"solver_type\": " << args.solver_type;
	os << ", \"thread_count\": " << args.thread_count;
	os << ", \"core_size\": " << args.core_size;
	os << "}";

	return os;
}
