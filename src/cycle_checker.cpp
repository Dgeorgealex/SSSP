#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "defs.h"
#include "graph.h"

namespace {

std::vector<NodeID> readCycle(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		ERROR("Could not open cycle file: " << filename);
	}

	std::size_t cycle_size = 0;
	if (!(file >> cycle_size)) {
		ERROR("Could not read cycle size from file: " << filename);
	}
	if (cycle_size == 0) {
		ERROR("Cycle size must be greater than 0.");
	}

	std::vector<NodeID> cycle(cycle_size);
	for (std::size_t i = 0; i < cycle_size; ++i) {
		if (!(file >> cycle[i])) {
			ERROR("Could not read cycle node " << i << " from file: " << filename);
		}
	}

	return cycle;
}

std::optional<Distance> minimumEdgeWeight(const Graph& graph, NodeID source, NodeID target) {
	std::optional<Distance> best_weight;
	for (const auto& edge : graph.getEdgesOf(source)) {
		if (edge.target != target) {
			continue;
		}
		if (!best_weight.has_value() || edge.weight < best_weight.value()) {
			best_weight = edge.weight;
		}
	}
	return best_weight;
}

void validateNoDuplicateNodes(const std::vector<NodeID>& cycle) {
	std::set<NodeID> seen;
	for (const auto& node : cycle) {
		if (seen.count(node)) {
			ERROR("Node " << node << " appears multiple times in the cycle.");
		}
		seen.insert(node);
	}
}

struct CycleStats {
	std::size_t negative_edges = 0;
	std::size_t positive_edges = 0;
	std::optional<Distance> max_positive_edge;
	Distance total_weight = 0;
};

std::optional<CycleStats> analyzeCycle(const Graph& graph, const std::vector<NodeID>& cycle) {
	CycleStats stats;
	for (std::size_t i = 0; i < cycle.size(); ++i) {
		const NodeID source = cycle[i];
		const NodeID target = cycle[(i + 1) % cycle.size()];
		auto edge_weight = minimumEdgeWeight(graph, source, target);
		if (!edge_weight.has_value()) {
			return {};
		}

		stats.total_weight += edge_weight.value();
		if (edge_weight.value() < 0) {
			++stats.negative_edges;
		} else if (edge_weight.value() > 0) {
			++stats.positive_edges;
			if (!stats.max_positive_edge.has_value() || edge_weight.value() > stats.max_positive_edge.value()) {
				stats.max_positive_edge = edge_weight.value();
			}
		}
	}
	return stats;
}

}  // namespace

int main(int argc, char** argv) {
	const std::string cycle_filename = "cycle.txt";
	std::string graph_filename = "../data/graphs/simple1.txt";
	if (argc > 1) {
		graph_filename = argv[1];
	}

	auto cycle = readCycle(cycle_filename);
	validateNoDuplicateNodes(cycle);
	auto graph = readGraph(graph_filename);
	auto stats = analyzeCycle(graph, cycle);
	if (!stats.has_value()) {
		ERROR("The cycle contains an edge that does not exist in the graph.");
	}
	auto cycle_stats = stats.value();

	if (cycle_stats.total_weight < 0) {
		std::cout << "negative cycle" << std::endl;
		std::cout << "cycle weight: " << cycle_stats.total_weight << std::endl;
		std::cout << "negative edges: " << cycle_stats.negative_edges << std::endl;
		std::cout << "positive edges: " << cycle_stats.positive_edges << std::endl;
		if (cycle_stats.max_positive_edge.has_value()) {
			std::cout << "max positive edge: " << cycle_stats.max_positive_edge.value() << std::endl;
		} else {
			std::cout << "max positive edge: none" << std::endl;
		}
		return 0;
	}

	std::cout << "not a negative cycle" << std::endl;
	std::cout << "cycle weight: " << cycle_stats.total_weight << std::endl;
	std::cout << "negative edges: " << cycle_stats.negative_edges << std::endl;
	std::cout << "positive edges: " << cycle_stats.positive_edges << std::endl;
	if (cycle_stats.max_positive_edge.has_value()) {
		std::cout << "max positive edge: " << cycle_stats.max_positive_edge.value() << std::endl;
	} else {
		std::cout << "max positive edge: none" << std::endl;
	}
	return 1;
}