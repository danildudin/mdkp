#pragma once
#include <iostream>
#include <vector>

struct LPSolution {
	struct Item {
		int id;
		double x, rc;
	};

	double cost;
	std::vector<Item> items;
	std::vector<Item*> items_map;
};

std::ostream& operator<<(std::ostream& os, const LPSolution::Item* item);
std::ostream& operator<<(std::ostream& os, const LPSolution::Item& item);
std::ostream& operator<<(std::ostream& os, const LPSolution& res);