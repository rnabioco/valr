#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <set>
#include <zlib.h>

using namespace Rcpp;

std::istream* open_file(std::string filename, std::ifstream& file_stream, gzFile& gz_file, bool& is_gzip) {
  if (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz") {
    is_gzip = true;
    gz_file = gzopen(filename.c_str(), "r");
    if (!gz_file) stop("Unable to open gzip file: " + filename);
    return nullptr; // Placeholder for gzip stream
  } else {
    is_gzip = false;
    file_stream.open(filename);
    if (!file_stream.is_open()) stop("Unable to open file: " + filename);
    return &file_stream;
  }
}

std::string gz_getline(gzFile gz_file) {
  char buffer[4096];
  std::string line;
  while (gzgets(gz_file, buffer, sizeof(buffer))) {
    line += buffer;
    if (!line.empty() && line.back() == '\n') {
      line.pop_back();
      break;
    }
  }
  return line;
}

// [[Rcpp::export]]
DataFrame read_gtf_impl(std::string filename) {
  // Cols for the GTF file
  std::vector<std::string> seqname, source, feature;
  std::vector<int> start, end;
  std::vector<double> score;
  std::vector<char> strand, frame;

  // For parsing attrs
  std::set<std::string> all_keys;
  std::vector<std::unordered_map<std::string, std::string>> row_attrs;

  std::ifstream file_stream;
  gzFile gz_file = nullptr;
  bool is_gzip;
  std::istream* input_stream = open_file(filename, file_stream, gz_file, is_gzip);

  std::string line;
  while (true) {
    if (is_gzip) {
      line = gz_getline(gz_file);
      if (line.empty() && gzeof(gz_file)) break;
    } else {
      if (!std::getline(*input_stream, line)) break;
    }

    if (line.empty() || line[0] == '#') continue; // Skip comments

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

  if (is_gzip) gzclose(gz_file);
  if (file_stream.is_open()) file_stream.close();

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
