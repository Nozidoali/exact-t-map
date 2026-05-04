#include "network/reader.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace exact_t {

namespace {

std::string get_extension(std::string const& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return "";
    return path.substr(dot);
}

std::string read_file(std::string const& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open file: " + path);
    return {std::istreambuf_iterator<char>(ifs),
            std::istreambuf_iterator<char>()};
}

std::string trim(std::string const& s) {
    auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::vector<std::string> split_comma(std::string const& s) {
    std::vector<std::string> result;
    std::istringstream iss(s);
    std::string token;
    while (std::getline(iss, token, ',')) {
        auto t = trim(token);
        if (!t.empty()) result.push_back(t);
    }
    return result;
}

std::vector<std::string> split_statements(std::string const& content) {
    std::vector<std::string> stmts;
    std::string current;
    for (char c : content) {
        if (c == ';') {
            auto t = trim(current);
            if (!t.empty()) stmts.push_back(t);
            current.clear();
        } else {
            current += c;
        }
    }
    auto t = trim(current);
    if (!t.empty()) stmts.push_back(t);
    return stmts;
}

struct ParsedSignal {
    std::string name;
    bool complemented{false};
};

ParsedSignal parse_signal_token(std::string const& tok) {
    ParsedSignal ps;
    if (!tok.empty() && tok[0] == '~') {
        ps.complemented = true;
        ps.name = tok.substr(1);
    } else {
        ps.name = tok;
    }
    return ps;
}

Signal lookup_signal(std::unordered_map<std::string, Signal> const& map,
                     ParsedSignal const& ps,
                     XagNetwork const& xag) {
    if (ps.name == "1'b0" || ps.name == "1\\'b0")
        return ps.complemented ? xag.get_constant(true) : xag.get_constant(false);
    if (ps.name == "1'b1" || ps.name == "1\\'b1")
        return ps.complemented ? xag.get_constant(false) : xag.get_constant(true);

    auto it = map.find(ps.name);
    if (it == map.end())
        throw std::runtime_error("Undefined signal: " + ps.name);
    return ps.complemented ? !it->second : it->second;
}

XagNetwork read_verilog(std::string const& path) {
    std::string content = read_file(path);
    auto stmts = split_statements(content);

    XagNetwork xag;
    std::unordered_map<std::string, Signal> sig;
    std::vector<std::string> output_names;

    for (auto const& stmt : stmts) {
        std::istringstream iss(stmt);
        std::string keyword;
        iss >> keyword;

        if (keyword == "module" || keyword == "endmodule") {
            continue;
        }

        if (keyword == "input") {
            std::string rest;
            std::getline(iss, rest);
            for (auto const& name : split_comma(rest)) {
                Signal pi = xag.create_pi();
                sig[name] = pi;
            }
            continue;
        }

        if (keyword == "output") {
            std::string rest;
            std::getline(iss, rest);
            for (auto const& name : split_comma(rest))
                output_names.push_back(name);
            continue;
        }

        if (keyword == "wire") {
            continue;
        }

        if (keyword == "assign") {
            std::string lhs, eq, tok1, op, tok2;
            iss >> lhs >> eq >> tok1;

            if (!(iss >> op)) {
                auto ps = parse_signal_token(tok1);
                sig[lhs] = lookup_signal(sig, ps, xag);
                continue;
            }

            iss >> tok2;
            auto ps0 = parse_signal_token(tok1);
            auto ps1 = parse_signal_token(tok2);
            Signal s0 = lookup_signal(sig, ps0, xag);
            Signal s1 = lookup_signal(sig, ps1, xag);

            Signal result(0, false);
            if (op == "&") {
                result = xag.create_and(s0, s1);
            } else if (op == "^") {
                result = xag.create_xor(s0, s1);
            } else if (op == "|") {
                result = !xag.create_and(!s0, !s1);
            } else {
                throw std::runtime_error("Unknown operator: " + op);
            }
            sig[lhs] = result;
            continue;
        }
    }

    for (auto const& name : output_names) {
        auto it = sig.find(name);
        if (it == sig.end())
            throw std::runtime_error("Output not defined: " + name);
        xag.create_po(it->second);
    }

    return xag;
}

uint32_t decode_aiger_delta(std::istream& in) {
    uint32_t x = 0;
    uint32_t shift = 0;
    int byte;
    do {
        byte = in.get();
        if (byte == EOF)
            throw std::runtime_error("Unexpected EOF in AIGER binary");
        x |= static_cast<uint32_t>(byte & 0x7F) << shift;
        shift += 7;
    } while (byte & 0x80);
    return x;
}

XagNetwork read_aiger_common(std::string const& path, bool binary) {
    std::ifstream ifs(path, binary ? std::ios::binary : std::ios::in);
    if (!ifs) throw std::runtime_error("Cannot open file: " + path);

    std::string header_tag;
    uint32_t M, I, L, O, A;
    ifs >> header_tag >> M >> I >> L >> O >> A;

    if (L != 0)
        throw std::runtime_error("Latches not supported in AIGER reader");

    std::string rest_of_header;
    std::getline(ifs, rest_of_header);

    std::vector<Signal> lit_to_signal(2 * (M + 1) + 1);
    lit_to_signal[0] = Signal(0, false);
    lit_to_signal[1] = Signal(0, true);

    XagNetwork xag;

    for (uint32_t i = 0; i < I; ++i) {
        Signal pi = xag.create_pi();
        uint32_t var = i + 1;
        lit_to_signal[2 * var] = pi;
        lit_to_signal[2 * var + 1] = !pi;
    }

    std::vector<uint32_t> output_lits(O);
    if (binary) {
        for (uint32_t i = 0; i < O; ++i) {
            ifs >> output_lits[i];
        }
        std::string nl;
        std::getline(ifs, nl);
    } else {
        std::vector<uint32_t> input_lits(I);
        for (uint32_t i = 0; i < I; ++i)
            ifs >> input_lits[i];
        for (uint32_t i = 0; i < O; ++i)
            ifs >> output_lits[i];
    }

    for (uint32_t i = 0; i < A; ++i) {
        uint32_t var = I + 1 + i;
        uint32_t lhs_lit = 2 * var;
        uint32_t rhs0_lit, rhs1_lit;

        if (binary) {
            uint32_t d0 = decode_aiger_delta(ifs);
            uint32_t d1 = decode_aiger_delta(ifs);
            rhs0_lit = lhs_lit - d0;
            rhs1_lit = rhs0_lit - d1;
        } else {
            uint32_t lhs_read;
            ifs >> lhs_read >> rhs0_lit >> rhs1_lit;
            assert(lhs_read == lhs_lit);
        }

        Signal s0 = lit_to_signal[rhs0_lit];
        Signal s1 = lit_to_signal[rhs1_lit];
        Signal gate = xag.create_and(s0, s1);
        lit_to_signal[lhs_lit] = gate;
        lit_to_signal[lhs_lit + 1] = !gate;
    }

    for (uint32_t i = 0; i < O; ++i)
        xag.create_po(lit_to_signal[output_lits[i]]);

    return xag;
}

}  // namespace

XagNetwork read_network(std::string const& path) {
    std::string ext = get_extension(path);

    if (ext == ".v" || ext == ".verilog")
        return read_verilog(path);

    if (ext == ".aig")
        return read_aiger_common(path, true);

    if (ext == ".aag")
        return read_aiger_common(path, false);

    throw std::runtime_error("Unsupported file format: " + ext +
                             " (supported: .v, .verilog, .aig, .aag)");
}

}  // namespace exact_t
