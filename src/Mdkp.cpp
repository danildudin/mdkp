#include <glpk.h>

#include "Mdkp.h"
#include "utils.h"

Mdkp::Mdkp(Problem& p):
	problem(p),
	cost(0),
	w_sum(0),
	b_sum(0),
	x(problem.n, Mdkp::CORE),
	n_size(LAST + 1),
	b(problem.b),
	weights(problem.m)
{
	n_size[Mdkp::CORE] = problem.n;
	for (auto val : b) {
		b_sum += val;
	}
}

Mdkp::Mdkp(const Mdkp& other):
	problem(other.problem),
	w_sum(other.w_sum),
	b_sum(other.b_sum),
	x(other.x),
	n_size(other.n_size),
	b(other.b),
	weights(other.weights)
{
	cost.store(other.cost);
}

const Mdkp& Mdkp::operator=(const Mdkp& other) {
	problem = other.problem;
	cost.store(other.cost);
	w_sum = other.w_sum;
	b_sum = other.b_sum;
	n_size = other.n_size;
	x = other.x;
	b = other.b;
	weights = other.weights;

	return *this;
}

std::vector<int> Mdkp::get_list(Mdkp::XType type) const {
	std::vector<int> res;
	for (int j = 0; j < x.size(); j++) {
		if (x[j] == type) res.push_back(j);
	}
	return res;
}

void Mdkp::set_x(int id, Mdkp::XType type) {
	if (x[id] == type) return;

	if (type == Mdkp::N1) {
		cost += problem.c[id];
		w_sum = 0;
		for (int i = 0; i < problem.m; i++) {
			weights[i] += problem.a[i][id];
			w_sum += weights[i];
		}
	} else if (x[id] == Mdkp::N1) {
		cost -= problem.c[id];
		w_sum = 0;
		for (int i = 0; i < problem.m; i++) {
			weights[i] -= problem.a[i][id];
			w_sum += weights[i];
		}
	}

	n_size[x[id]]--;
	n_size[type]++;
	x[id] = type;
}

bool Mdkp::is_n0(int id) {
	return x[id] == Mdkp::N0;
}

bool Mdkp::is_n1(int id) {
	return x[id] == Mdkp::N1;
}

bool Mdkp::is_core(int id) {
	return x[id] == Mdkp::CORE;
}

bool Mdkp::is_feasible() const {
	for (int i = 0; i < problem.m; i++) {
		if (weights[i] > b[i]) {
			return false;
		}
	}

	return true;
}

LPSolution Mdkp::lp_relaxation() const {
	int m = problem.m;
	int n = n_size[Mdkp::CORE];
	std::vector<int> ids = get_list(Mdkp::CORE);

	int ia[1 + n * m];
	int ja[1 + n * m];
	double ar[1 + n * m];

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			ia[i * n + j + 1] = i + 1;
			ja[i * n + j + 1] = j + 1;
			ar[i * n + j + 1] = problem.a[i][ids[j]];
		}
	}

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, m);
	for (int i = 0; i < m; i++) {
		glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, b[i] - weights[i]);
	}

	glp_add_cols(lp, n);
	for (int j = 0; j < n; j++) {
		glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, j + 1, problem.c[ids[j]]);
	}

	glp_load_matrix(lp, n * m, ia, ja, ar);

	glp_smcp params;
	glp_init_smcp(&params);
	params.msg_lev = GLP_MSG_ERR;

	glp_simplex(lp, &params);

	LPSolution res;
	res.items_map.resize(problem.n);
	res.cost = glp_get_obj_val(lp);
	for (int j = 0; j < n; j++) {
		res.items.emplace_back(ids[j], glp_get_col_prim(lp, j + 1), glp_get_col_dual(lp, j + 1));
	}
	std::sort(res.items.begin(), res.items.end(), [&](const LPSolution::Item &a, const LPSolution::Item &b) { return abs(a.rc) < abs(b.rc); });
	
	for (int j = 0; j < n; j++) {
		res.items_map[res.items[j].id] = &res.items[j];
	}

	glp_delete_prob(lp);
	return res;
}

std::ostream& operator<<(std::ostream& os, const Mdkp::XType xtype) {
	switch (xtype) {
	case Mdkp::N1: os << "\"Mdkp::N1\"";
	case Mdkp::N0: os << "\"Mdkp::N0\"";
	case Mdkp::CORE: os << "\"Mdkp::CORE\"";
	default: os << "\"ERR\"";
	}

	return os;
}

std::ostream& operator<<(std::ostream& os, const Mdkp& mdkp) {
    os << "{";
    os << "\"is_feasible()\": " << mdkp.is_feasible();
    os << ", \"cost\": " << mdkp.cost;
    os << ", \"n_size\": " << mdkp.n_size;
    os << ", \"x\": " << mdkp.x;
    os << ", \"w_sum\": " << mdkp.w_sum;
    os << ", \"b_sum\": " << mdkp.b_sum;
    os << ", \"weights\": " << mdkp.weights;
    os << ", \"b\": " << mdkp.b;
    os << ", \"problem\": " << mdkp.problem;
    os << "}";

    return os;
}