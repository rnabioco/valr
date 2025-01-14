#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <set>

using namespace Rcpp;

//[[Rcpp::export]]
DataFrame read_gtf_impl(std::string filename) {
  // Cols for the GTF file
  std::vector<std::string> seqname, source, feature;
  std::vector<int> start, end;
  std::vector<double> score;
  std::vector<char> strand, frame;

  // For parsing attrs
  std::set<std::string> all_keys;
  std::vector<std::unordered_map<std::string, std::string>> row_attrs;

  std::ifstream file(filename);
  if (!file.is_open()) {
    stop("Unable to open file: " + filename);
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line[0] == '#') continue; // Skip comments

    std::istringstream linestream(line);
    std::string seq, src, feat, attr;
    int st, en;
    double sc;
    char str, frm;

    linestream >> seq >> src >> feat >> st >> en >> sc >> str >> frm;
    std::getline(linestream, attr);

    // Trim leading spaces from attrs
    attr.erase(0, attr.find_first_not_of(" \t"));

    // Populate basic cols
    seqname.push_back(seq);
    source.push_back(src);
    feature.push_back(feat);
    start.push_back(st);
    end.push_back(en);
    score.push_back(sc);
    strand.push_back(str);
    frame.push_back(frm);

    // Parse attrs into a map for this row
    std::unordered_map<std::string, std::string> attr_map;
    std::istringstream attr_stream(attr);
    std::string key_value_pair;
    while (std::getline(attr_stream, key_value_pair, ';')) {
      size_t eq_pos = key_value_pair.find(' ');
      if (eq_pos == std::string::npos) continue;

      std::string key = key_value_pair.substr(0, eq_pos);
      std::string value = key_value_pair.substr(eq_pos + 1);

      // Trim spaces and quotes
      key.erase(0, key.find_first_not_of(" \t"));
      key.erase(key.find_last_not_of(" \t") + 1);

      value.erase(0, value.find_first_not_of(" \t\""));
      value.erase(value.find_last_not_of(" \t\"") + 1);

      attr_map[key] = value;
      all_keys.insert(key);
    }
    row_attrs.push_back(attr_map);
  }
  file.close();

  // Initialize attr cols
  std::unordered_map<std::string, std::vector<std::string>> attr_cols;
  for (const auto &key : all_keys) {
    attr_cols[key] = std::vector<std::string>(seqname.size(), std::string());
  }

  // Populate attr cols
  for (size_t i = 0; i < row_attrs.size(); ++i) {
    for (const auto &pair : row_attrs[i]) {
      attr_cols[pair.first][i] = pair.second;
    }
  }

  // Create the DataFrame
  List df = List::create(
    Named("seqname") = seqname,
    Named("source") = source,
    Named("feature") = feature,
    Named("start") = start,
    Named("end") = end,
    Named("score") = score,
    Named("strand") = strand,
    Named("frame") = frame
  );

  // Add attr cols as named entries
  for (const auto &col : attr_cols) {
    df[col.first] = col.second;
  }

  return DataFrame(df);
}
