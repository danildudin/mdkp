#include "BFSolver.h"

void BFSolver::search_tree(Mdkp &cur, int id) {
	if (id == problem.c.size()) {
		if (cur.is_feasible() && cur.cost > res.cost) {
			res = cur;
		}
		return;
	}

	search_tree(cur, id + 1);
	cur.set_x(id, Mdkp::N1);
	search_tree(cur, id + 1);
	cur.set_x(id, Mdkp::N0);
}

void BFSolver::solve() {
	Mdkp cur(problem);
	search_tree(cur, 0);
}
