#include "Problem.h"
#include "utils.h"

void Problem::init() {
	int real_ans;
	std::cin >> n >> m >> real_ans;

	a.resize(m, std::vector<int>(n));
	a_sum.resize(n);
	b.resize(m);
	c.resize(n);
	for (int j = 0; j < n; j++) {
		std::cin >> c[j];
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cin >> a[i][j];
			a_sum[j] += a[i][j];
		}
	}
	for (int i = 0; i < m; i++) {
		std::cin >> b[i];
	}
}

// std::ostream& operator<<(std::ostream& os, const Problem& p) {
//     os << "{";
//     os << "\"n\": " << p.n << ", \"m\": " << p.m;
//     os << ", \"c\": " << p.c;
//     os << ", \"b\": " << p.b;
//     os << ", \"a_sum\": " << p.a_sum;
//     os << ", \"a\": " << p.a;
//     os << "}";
//     return os;
// }
