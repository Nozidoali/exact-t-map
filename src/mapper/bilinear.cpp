#include "mapper/bilinear.hpp"

namespace exact_t {

BilinearBuild build_bilinear(CutMatch const& match,
                              std::vector<xag_node> const& input_nodes,
                              std::set<xag_node> const& unique_outputs) {
    BilinearBuild result{
        BilinearFunction(static_cast<uint32_t>(input_nodes.size()),
                         static_cast<uint32_t>(unique_outputs.size())),
        {}, {}};

    for (size_t i = 0; i < input_nodes.size(); ++i)
        result.node_to_input[input_nodes[i]] = static_cast<uint32_t>(i);

    uint32_t out_idx = 0;
    for (xag_node out : unique_outputs)
        result.node_to_output[out] = out_idx++;

    for (size_t i = 0; i < match.anfs.size() && i < match.outputs.size(); ++i) {
        auto it_out = result.node_to_output.find(match.outputs[i]);
        if (it_out == result.node_to_output.end()) continue;
        for (auto const& p : match.anfs[i].degree2) {
            auto it_i = result.node_to_input.find(p.first);
            auto it_j = result.node_to_input.find(p.second);
            if (it_i != result.node_to_input.end() &&
                it_j != result.node_to_input.end())
                result.bf.add_triplet(it_i->second, it_j->second, it_out->second);
        }
    }
    return result;
}

}  // namespace exact_t
