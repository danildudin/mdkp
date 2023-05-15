#pragma once
#include <iostream>
#include <vector>

class Problem {
public:
	void init();

public:
	int n, m;
	std::vector<std::vector<int>> a;
	std::vector<int> c, b, a_sum;
};

// std::ostream& operator<<(std::ostream& os, const Problem& p);
