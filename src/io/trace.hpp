/*! \file trace.hpp
 *  \brief JSONL trace sink for per-node / per-round mapping events.
 */

#pragma once

#include <fstream>
#include <memory>
#include <string>

namespace exact_t {

/*! \brief Abstract sink for a line of JSON output. */
class TraceSink {
  public:
    virtual ~TraceSink() = default;
    virtual void emit(std::string const& line) = 0;
};

/*! \brief No-op sink. All emit() calls discard the input. */
class NullSink : public TraceSink {
  public:
    void emit(std::string const&) override {}
};

/*! \brief Appends each line to a file, flushing after every emit.
 *
 *  On crash the file still contains well-formed JSON lines up to the
 *  last flush boundary.
 */
class JsonlFileSink : public TraceSink {
  public:
    explicit JsonlFileSink(std::string const& path);
    void emit(std::string const& line) override;

  private:
    std::ofstream out_;
};

/*! \brief Factory: empty path -> NullSink, else -> JsonlFileSink. */
std::unique_ptr<TraceSink> make_trace_sink(std::string const& path);

}  // namespace exact_t
