#include "CoreData.h"
#include "utils.h"

CoreData::CoreData(const Mdkp &mdkp) {
	k = 0;
	start_cost = mdkp.cost;
	lp_res = mdkp.lp_relaxation();

	sorted = w_sum_sorted = mdkp.get_list(Mdkp::CORE);
	std::sort(sorted.begin(), sorted.end(), [&](int a, int b) { return mdkp.problem.c[a] > mdkp.problem.c[b]; });
	std::sort(w_sum_sorted.begin(), w_sum_sorted.end(), [&](int a, int b) { return mdkp.problem.a_sum[a] < mdkp.problem.a_sum[b]; });
	for (int i = 0; i < mdkp.problem.m; i++) {
		auto cur = sorted;
		std::sort(cur.begin(), cur.end(), [&](int a, int b) { return mdkp.problem.a[i][a] < mdkp.problem.a[i][b]; });
		w_sorted.push_back(std::move(cur));
	}
}

// std::ostream& operator<<(std::ostream& os, const CoreData& data) {
//     os << "{";
//     os << "\"k\": " << data.k;
//     os << ", \"start_cost\": " << data.start_cost;
//     os << ", \"sorted\": " << data.sorted;
//     os << ", \"w_sum_sorted\": " << data.w_sum_sorted;
//     os << ", \"lp_res\": " << data.lp_res;
//     os << "}";

//     return os;
// }
