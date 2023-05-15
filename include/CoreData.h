#pragma once
#include <iostream>
#include <vector>

#include "LPSolution.h"
#include "Mdkp.h"

class CoreData {
public:
	CoreData(const Mdkp &mdkp);

public:
	int k, start_cost;
	std::vector<int> sorted, w_sum_sorted;
	std::vector<std::vector<int>> w_sorted;
	LPSolution lp_res;
};

std::ostream& operator<<(std::ostream& os, const CoreData& data);
