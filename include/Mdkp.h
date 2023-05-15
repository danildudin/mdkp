#pragma once
#include <iostream>
#include <atomic>

#include "Problem.h"
#include "LPSolution.h"

class Mdkp {
public:
	enum XType {N1, N0, CORE, LAST = CORE};

	Mdkp(Problem& p);
	Mdkp(const Mdkp& other);

	const Mdkp& operator=(const Mdkp& other);

	bool is_feasible() const;
	bool is_n0(int id);
	bool is_n1(int id);
	bool is_core(int id);
	void set_x(int id, XType type);
	std::vector<int> get_list(XType type) const;
	LPSolution lp_relaxation() const;

public:
	Problem& problem;
	std::atomic_int cost;
	int w_sum, b_sum;
	std::vector<XType> x;
	std::vector<int> n_size, b, weights;
};

// std::ostream& operator<<(std::ostream& os, const Mdkp::XType xtype);
// std::ostream& operator<<(std::ostream& os, const Mdkp& mdkp);
