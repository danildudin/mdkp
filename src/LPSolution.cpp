#include "LPSolution.h"
#include "utils.h"

std::ostream& operator<<(std::ostream& os, const LPSolution::Item* item) {
	if (!item) {
		os << "\"NULL\"";
	} else {
    	os << "{\"ptr\": 1, \"id\": " << item->id << ", \"x\": " << item->x << ", \"rc\": " << item->rc << "}";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const LPSolution::Item& item) {
    os << "{\"id\": " << item.id << ", \"x\": " << item.x << ", \"rc\": " << item.rc << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const LPSolution& res) {
    os << "{\"cost\": " << res.cost << ", \"items\": " << res.items;
    os  << ", \"items_map\": " << res.items_map;
    os << "}";
    return os;
}
