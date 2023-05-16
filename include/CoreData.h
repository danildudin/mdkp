#pragma once
#include <iostream>
#include <vector>

#include "LPSolution.h"
#include "Mdkp.h"

class CoreData {
public:
	CoreData() {};
	CoreData(const Mdkp &mdkp);

public:
	int k = 0;
	int start_cost = 0;
	std::vector<int> sorted, w_sum_sorted;
	std::vector<std::vector<int>> w_sorted;
	LPSolution lp_res;
};

std::ostream& operator<<(std::ostream& os, const CoreData& data);
