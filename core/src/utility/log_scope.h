#pragma once

#define LOG_SCOPE_MARKER_START "Scope Start "
#define LOG_SCOPE_MARKER_END "Scope End "

#define WHERE std::string(__FILE__) + std::string(":") + std::to_string(__LINE__) + std::string(" ")  + std::string(__PRETTY_FUNCTION__)
#define LOG_SCOPE_FUNCTION LogScope log_scope_function_guard (WHERE);
/*
use as
  LogScope s (WHERE);
or
  LOG_SCOPE_FUNCTION
at start of function.
*/
struct LogScope {
    LogScope(std::string n) : name(n) {
#ifndef DISABLE_LOG_SCOPE
      std::cout << LOG_SCOPE_MARKER_START << name << std::endl;
#endif
    }
    ~LogScope() {
#ifndef DISABLE_LOG_SCOPE
      std::cout << LOG_SCOPE_MARKER_END << name << std::endl;
#endif
    }
private:
    std::string name;
};
