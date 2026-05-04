#include "io/trace.hpp"

#include <cassert>

namespace exact_t {

JsonlFileSink::JsonlFileSink(std::string const& path) : out_(path) {
    assert(out_.is_open());
}

void JsonlFileSink::emit(std::string const& line) {
    out_ << line << '\n';
    out_.flush();
}

std::unique_ptr<TraceSink> make_trace_sink(std::string const& path) {
    if (path.empty()) return std::make_unique<NullSink>();
    return std::make_unique<JsonlFileSink>(path);
}

}  // namespace exact_t
