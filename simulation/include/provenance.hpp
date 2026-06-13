#pragma once
#include <cstdio>
#include <ctime>
#include <ostream>
#include <sstream>
#include <string>
#include <yaml-cpp/yaml.h>

namespace provenance {

inline std::string run_cmd(const char *cmd) {
  std::string out;
  if (FILE *p = popen(cmd, "r")) {
    char buf[256];
    if (fgets(buf, sizeof buf, p)) out.assign(buf);
    pclose(p);
  }
  while (!out.empty() && (out.back() == '\n' || out.back() == '\r'))
    out.pop_back();
  return out;
}

inline std::string git_describe() {
  std::string hash = run_cmd("git rev-parse --short HEAD 2>/dev/null");
  if (hash.empty()) return "unknown";
  if (!run_cmd("git status --porcelain 2>/dev/null").empty()) hash += "-dirty";
  return hash;
}

inline std::string timestamp_utc() {
  std::time_t now = std::time(nullptr);
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
  return buf;
}

inline void write_header(std::ostream &out, const std::string &tool,
                         const std::string &config_path,
                         const YAML::Node &config) {
  out << "# tool: " << tool << "\n"
      << "# generated: " << timestamp_utc() << "\n"
      << "# git: " << git_describe() << "\n"
      << "# config: " << config_path << "\n";
  std::istringstream dump(YAML::Dump(config));
  for (std::string line; std::getline(dump, line);)
    out << "#   " << line << "\n";
}

} // namespace provenance
