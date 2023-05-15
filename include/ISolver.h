#pragma once

class ISolver {
public:
	virtual void init() = 0;
	virtual void solve() = 0;
	virtual void print_solution() = 0;
};
