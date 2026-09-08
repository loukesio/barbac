#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

// =============================================================================
// Build marker
// =============================================================================
const char* BUILD_ID = "barbac-2026-09-08-exact-partitions-v12";

// =============================================================================
// Distance routines
// =============================================================================
int levenshtein_banded(const char* s1, int len1,
                       const char* s2, int len2,
                       int max_dist) {
  if (std::abs(len1 - len2) > max_dist) return max_dist + 1;
  if (len1 == 0) return len2;
  if (len2 == 0) return len1;
  
  if (len1 > len2) {
    std::swap(s1, s2);
    std::swap(len1, len2);
  }
  
  const int INF = max_dist + 1;
  const int MAX_STACK = 96;
  int prev_stack[MAX_STACK + 1];
  int curr_stack[MAX_STACK + 1];
  std::vector<int> prev_heap, curr_heap;
  int* prev;
  int* curr;
  
  if (len2 <= MAX_STACK) {
    prev = prev_stack;
    curr = curr_stack;
  } else {
    prev_heap.resize(len2 + 1);
    curr_heap.resize(len2 + 1);
    prev = prev_heap.data();
    curr = curr_heap.data();
  }
  
  for (int j = 0; j <= len2; ++j) prev[j] = (j <= max_dist) ? j : INF;
  
  for (int i = 1; i <= len1; ++i) {
    int j_min = std::max(1, i - max_dist);
    int j_max = std::min(len2, i + max_dist);
    
    curr[0] = (i <= max_dist) ? i : INF;
    if (j_min > 1) curr[j_min - 1] = INF;
    if (j_max < len2) curr[j_max + 1] = INF;
    
    int row_min = INF;
    const char c1 = s1[i - 1];
    for (int j = j_min; j <= j_max; ++j) {
      const int cost = (c1 == s2[j - 1]) ? 0 : 1;
      const int del = prev[j] + 1;
      const int ins = curr[j - 1] + 1;
      const int sub = prev[j - 1] + cost;
      int val = del < ins ? del : ins;
      val = val < sub ? val : sub;
      curr[j] = val;
      if (val < row_min) row_min = val;
    }
    
    if (row_min > max_dist) return INF;
    std::swap(prev, curr);
  }
  
  return prev[len2];
}


int levenshtein_myers_64(const char* pattern, int m,
                         const char* text, int n,
                         int max_dist) {
  if (std::abs(m - n) > max_dist) return max_dist + 1;
  if (m == 0) return n;
  if (n == 0) return m;
  if (m > 63) return levenshtein_banded(pattern, m, text, n, max_dist);
  
  // Use the shorter string as the bit-vector pattern. This keeps the bit mask
  // compact and makes the routine symmetric for our purposes.
  if (m > n) {
    std::swap(pattern, text);
    std::swap(m, n);
  }
  
  uint64_t peq[256];
  for (int i = 0; i < 256; ++i) peq[i] = 0;
  for (int i = 0; i < m; ++i) {
    peq[static_cast<unsigned char>(pattern[i])] |= (1ULL << i);
  }
  
  uint64_t pv = ~0ULL;
  uint64_t mv = 0;
  int score = m;
  const uint64_t high_bit = 1ULL << (m - 1);
  
  for (int i = 0; i < n; ++i) {
    const uint64_t eq = peq[static_cast<unsigned char>(text[i])];
    const uint64_t xv = eq | mv;
    const uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;
    uint64_t ph = mv | ~(xh | pv);
    uint64_t mh = pv & xh;
    
    if (ph & high_bit) {
      ++score;
    } else if (mh & high_bit) {
      --score;
    }
    
    ph = (ph << 1) | 1ULL;
    mh <<= 1;
    pv = mh | ~(xv | ph);
    mv = ph & xv;
  }
  
  return score <= max_dist ? score : max_dist + 1;
}

inline int levenshtein_fast(const char* s1, int len1,
                            const char* s2, int len2,
                            int max_dist) {
  if (std::max(len1, len2) <= 63) {
    return levenshtein_myers_64(s1, len1, s2, len2, max_dist);
  }
  return levenshtein_banded(s1, len1, s2, len2, max_dist);
}

// =============================================================================
// 2-bit DNA packing
// =============================================================================
inline int dna_base(char c) {
  switch (c) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

inline bool pack_seq(const char* s, int len, uint64_t& out) {
  if (len < 0 || len > 32) return false;
  uint64_t v = 0;
  for (int i = 0; i < len; ++i) {
    const int b = dna_base(s[i]);
    if (b < 0) return false;
    v = (v << 2) | static_cast<uint64_t>(b);
  }
  out = v;
  return true;
}

inline int hamming_packed(uint64_t a, uint64_t b) {
  uint64_t x = a ^ b;
  x = (x | (x >> 1)) & 0x5555555555555555ULL;
  return __builtin_popcountll(x);
}

inline bool hamming_bit_prefilter(uint64_t a, uint64_t b, int max_dist) {
  // Safe only for Hamming rejection. If more than 2D packed bits differ,
  // more than D bases must differ.
  return __builtin_popcountll(a ^ b) <= 2 * max_dist;
}

// Every edit changes the A/C/G/T composition vector by an L1 distance of at
// most two: two for a substitution and one for an insertion or deletion.
// Therefore composition_l1 > 2D proves that Levenshtein distance is > D.
// Counting only A/C/G/T remains safe for sequences containing other symbols;
// it can make the bound weaker, never reject a reachable sequence.
using BaseComposition = std::array<int, 4>;

inline BaseComposition base_composition(const char* seq, int len) {
  BaseComposition out = {{0, 0, 0, 0}};
  for (int i = 0; i < len; ++i) {
    const int base = dna_base(seq[i]);
    if (base >= 0) ++out[base];
  }
  return out;
}

inline int composition_l1(const BaseComposition& a, const BaseComposition& b) {
  return std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]) +
    std::abs(a[2] - b[2]) + std::abs(a[3] - b[3]);
}

// =============================================================================
// Method parsing
// =============================================================================
bool parse_lv_method(const std::string& method) {
  if (method == "lv" || method == "levenshtein" || method == "edit") return true;
  if (method == "hamming" || method == "ham") return false;
  stop("Unknown method '%s'. Use 'lv'/'levenshtein' or 'hamming'/'ham'.", method);
  return true;
}

// =============================================================================
// Candidate indexes
// =============================================================================
class CandidateAccumulator {
  mutable std::vector<int> hit_counts;
  mutable std::vector<int> touched;
  
public:
  void add_hit(int idx, int n_centroids) const {
    if (idx < 0 || idx >= n_centroids) return;
    if (idx >= static_cast<int>(hit_counts.size())) hit_counts.resize(idx + 1, 0);
    if (hit_counts[idx] == 0) touched.push_back(idx);
    ++hit_counts[idx];
  }
  
  void flush(int min_hits, std::vector<int>& out) const {
    out.clear();
    out.reserve(touched.size());
    for (int idx : touched) {
      if (hit_counts[idx] >= min_hits) out.push_back(idx);
      hit_counts[idx] = 0;
    }
    touched.clear();
  }
};

// A fixed partition of the query into D+1 disjoint blocks guarantees one
// untouched block under D edits. For LV its position can shift by at most D.
// Partition boundaries affect search cost only, never the candidate guarantee.
// Learn collision probabilities from observed sequences, without a template or
// truth labels, and avoid putting an entire block inside a constant anchor.
class ExactPartitionIndex {
  struct Block {
    std::vector<int> positions;
    std::unordered_map<uint64_t, std::vector<int> > buckets;
  };
  struct Plan { std::vector<Block> blocks; };
  std::vector<Plan> plans;
  std::vector<int> fallback;
  int D;
  bool lv;
  mutable CandidateAccumulator acc;

  bool supported(const char* seq, int len) const {
    if (len <= D || len > (lv ? 64 : 32)) return false;
    for (int p = 0; p < len; ++p) {
      if (dna_base(seq[p]) < 0) return false;
      // Packed Hamming is case insensitive; edit distance compares characters.
      if (lv && seq[p] != "ACGT"[dna_base(seq[p])]) return false;
    }
    return true;
  }

  uint64_t code_at(const char* seq, const Block& block, int shift = 0) const {
    uint64_t code = 0;
    for (int p : block.positions)
      code = (code << 2) | static_cast<uint64_t>(dna_base(seq[p + shift]));
    return code;
  }

public:
  ExactPartitionIndex(int distance, bool is_lv,
                      const std::vector<const char*>& seqs,
                      const std::vector<int>& lengths,
                      const IntegerVector& counts, bool enabled = true) : plans(65), D(distance), lv(is_lv) {
    if (!enabled) return;
    std::vector<std::vector<std::array<double, 4> > > hist(65);
    for (size_t i = 0; i < seqs.size(); ++i) {
      int len = lengths[i];
      if (!supported(seqs[i], len)) continue;
      if (hist[len].empty()) hist[len].resize(len, {{0, 0, 0, 0}});
      for (int p = 0; p < len; ++p)
        hist[len][p][dna_base(seqs[i][p])] += counts[i];
    }
    for (int len = 1; len < 65; ++len) {
      if (hist[len].empty()) continue;
      const int B = D + 1;
      std::vector<double> info(len);
      for (int p = 0; p < len; ++p) {
        double total = 0, squares = 0;
        for (double x : hist[len][p]) { total += x; squares += x*x; }
        info[p] = total > 0 ? -std::log(std::max(1e-12, squares/(total*total))) : 0;
      }
      Plan& plan = plans[len];
      plan.blocks.resize(B);
      if (!lv) {
        // Hamming blocks need not be contiguous. Spread informative positions
        // over the blocks; fixed anchor positions then cost no extra postings.
        std::vector<int> order(len);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return info[a] > info[b]; });
        std::vector<double> load(B, 0);
        for (int p : order) {
          int best = 0;
          for (int b = 1; b < B; ++b)
            if (load[b] < load[best] ||
                (load[b] == load[best] && plan.blocks[b].positions.size() < plan.blocks[best].positions.size())) best = b;
          plan.blocks[best].positions.push_back(p);
          load[best] += info[p];
        }
      } else {
        // Minimise expected total posting-list size. LV blocks are contiguous
        // so an unchanged substring survives insertion/deletion shifts.
        std::vector<double> prefix(len + 1, 0);
        for (int p = 0; p < len; ++p) prefix[p+1] = prefix[p] + info[p];
        std::vector<std::vector<double> > cost(B+1, std::vector<double>(len+1, std::numeric_limits<double>::infinity()));
        std::vector<std::vector<int> > cut(B+1, std::vector<int>(len+1, -1));
        cost[0][0] = 0;
        for (int b = 1; b <= B; ++b)
          for (int end = b; end <= len; ++end)
            for (int start = std::max(b-1, end-32); start < end; ++start) {
              double value = cost[b-1][start] + std::exp(-(prefix[end]-prefix[start]));
              if (value < cost[b][end]) { cost[b][end] = value; cut[b][end] = start; }
            }
        if (cut[B][len] < 0) { plan.blocks.clear(); continue; }
        int end = len;
        for (int b = B; b > 0; --b) {
          int start = cut[b][end];
          for (int p = start; p < end; ++p) plan.blocks[b-1].positions.push_back(p);
          end = start;
        }
      }
    }
  }

  void add(const char* seq, int len, int id) {
    if (!supported(seq, len)) { fallback.push_back(id); return; }
    int first = lv ? std::max(1, len-D) : len;
    int last = lv ? std::min(64, len+D) : len;
    for (int qlen = first; qlen <= last; ++qlen) {
      for (Block& block : plans[qlen].blocks) {
        int lo = lv ? std::max(-D, -block.positions.front()) : 0;
        int hi = lv ? std::min(D, len-1-block.positions.back()) : 0;
        for (int shift = lo; shift <= hi; ++shift) {
          std::vector<int>& bucket = block.buckets[code_at(seq, block, shift)];
          // Multiple shifts of a repeated substring must count only once.
          if (bucket.empty() || bucket.back() != id) bucket.push_back(id);
        }
      }
    }
  }

  void query(const char* seq, int len, int n_centroids, std::vector<int>& out) const {
    if (!supported(seq, len) || plans[len].blocks.empty()) {
      out.resize(n_centroids);
      std::iota(out.begin(), out.end(), 0);
      return;
    }
    for (const Block& block : plans[len].blocks) {
      auto it = block.buckets.find(code_at(seq, block));
      if (it != block.buckets.end())
        for (int id : it->second) acc.add_hit(id, n_centroids);
    }
    for (int id : fallback) acc.add_hit(id, n_centroids);
    acc.flush(1, out);

  }
};

// =============================================================================
// Clustering model
// =============================================================================
struct Cluster {
  int centroid_idx;
  int centroid_len;
  uint64_t centroid_pack;
  bool centroid_packable;
  int centroid_count;
  int sum_counts;
  std::vector<int> members;
  std::vector<int> member_dists;
  
  Cluster() : centroid_idx(-1), centroid_len(0), centroid_pack(0),
  centroid_packable(false), centroid_count(0), sum_counts(0) {}
};

struct CandidateChoice {
  int cluster_id;
  int dist;
  double score;
  
  CandidateChoice() : cluster_id(-1), dist(0), score(-std::numeric_limits<double>::infinity()) {}
};

inline double clamp_error_rate(double error_rate) {
  if (!std::isfinite(error_rate) || error_rate <= 0.0) return 0.005;
  if (error_rate >= 0.25) return 0.25;
  return error_rate;
}

inline double distance_log_likelihood_score(int dist, int parent_count,
                                            int child_count, int len,
                                            double error_rate,
                                            bool is_lv) {
  // Score for choosing the most plausible absorbing parent. This is deliberately
  // count-aware so a massive d=2 parent can beat a tiny d=1 parent when warranted.
  const double e = clamp_error_rate(error_rate);
  const double alphabet = is_lv ? 5.0 : 3.0; // LV has substitutions plus indel paths.
  const double edit_penalty = -std::log(e / alphabet);
  const double noedit_bonus = std::log(std::max(1e-12, 1.0 - e));
  
  return std::log1p(static_cast<double>(parent_count)) -
    static_cast<double>(dist) * edit_penalty +
    static_cast<double>(std::max(0, len - dist)) * noedit_bonus -
    0.15 * std::log1p(static_cast<double>(child_count));
}

// =============================================================================
// Barcode design
// =============================================================================
// A designed library fixes some positions and randomises others. A read that
// differs from a centroid only at a fixed position cannot be a different
// barcode -- no barcode varies there -- so the difference is a sequencing
// error and the read belongs to that centroid. A read differing at a random
// position may genuinely be another barcode.
//
// Edit distance cannot express that: it counts both mismatches as one. Knowing
// which positions carry identity is information the distance metric does not
// have, and it is recoverable from the library itself, because a fixed
// position shows one base in nearly every read while a random one shows four.
struct DesignMask {
  uint64_t variable_bits;   // packed layout, low bit of each 2-bit lane
  int len;
  int n_variable;
  bool usable;
};

// A position counts as fixed when one base covers at least this share of reads.
// Sequencing error puts a designed anchor near 1 - error_rate; a randomised
// position sits near 0.25, so the two are far apart and the threshold is not
// delicate.
const double DESIGN_FIXED_SHARE = 0.90;

DesignMask learn_design(const std::vector<const char*>& ptr,
                        const std::vector<int>& slen,
                        const IntegerVector& counts,
                        int n, int modal_len) {
  DesignMask m;
  m.variable_bits = 0;
  m.len = modal_len;
  m.n_variable = 0;
  m.usable = false;
  if (modal_len <= 0 || modal_len > 32) return m;

  std::vector<std::array<double, 4> > freq(modal_len, {{0.0, 0.0, 0.0, 0.0}});
  std::vector<double> total(modal_len, 0.0);
  for (int i = 0; i < n; ++i) {
    if (slen[i] != modal_len) continue;
    const double w = static_cast<double>(std::max(1, counts[i]));
    for (int p = 0; p < modal_len; ++p) {
      const int b = dna_base(ptr[i][p]);
      if (b < 0) continue;
      freq[p][b] += w;
      total[p] += w;
    }
  }

  for (int p = 0; p < modal_len; ++p) {
    if (total[p] <= 0.0) continue;
    double best = 0.0;
    for (int b = 0; b < 4; ++b) best = std::max(best, freq[p][b]);
    if (best / total[p] < DESIGN_FIXED_SHARE) {
      // pack_seq writes position 0 into the highest lane
      m.variable_bits |= (1ULL << (2 * (modal_len - 1 - p)));
      ++m.n_variable;
    }
  }

  // Only meaningful when the design actually fixes something. An all-random
  // library has nothing to exploit, and treating it as designed would be a
  // licence to merge on noise.
  m.usable = (m.n_variable > 0 && m.n_variable < modal_len);
  return m;
}

// Mismatches that fall on identity-carrying positions. The rest are errors by
// construction, so they should not count towards "these may be two barcodes".
inline int variable_mismatches(const DesignMask& m, uint64_t a, uint64_t b) {
  uint64_t x = a ^ b;
  x = (x | (x >> 1)) & 0x5555555555555555ULL;   // one low bit per differing base
  return __builtin_popcountll(x & m.variable_bits);
}

inline double effective_merge_ratio(int dist, double base_ratio, bool is_lv) {
  if (dist <= 0) return 0.0;
  
  if (!is_lv) {
    // Hamming mode was creating many low-count FP/WS centroids when the d=2/d=3
    // guard was too aggressive. Keep Hamming close to the proven ratio rule:
    // protect only reasonably abundant children, and use the same ratio at all
    // distances so count-1..9 error reads still get absorbed.
    return base_ratio;
  }
  
  // LV mode benefits from a distance-aware guard because a d=3 edit-distance
  // child is less likely to be an error than a d=1 child, especially when indels
  // can separate nearby true barcodes.
  if (dist == 1) return base_ratio;
  if (dist == 2) return base_ratio * 3.0;
  if (dist == 3) return base_ratio * 5.0;
  return base_ratio * 8.0;
}

inline int effective_count_floor(int dist, bool is_lv) {
  if (dist <= 0) return 0;
  if (!is_lv) return 10;
  if (dist == 1) return 10;
  if (dist == 2) return 5;
  return 2;
}

inline bool merge_blocked_by_ratio(int dist, int child_count, int parent_count,
                                   double base_ratio, bool is_lv) {
  const int floor = effective_count_floor(dist, is_lv);
  if (child_count < floor) return false;
  const double required = effective_merge_ratio(dist, base_ratio, is_lv) *
    static_cast<double>(child_count);
  return static_cast<double>(parent_count) < required;
}

inline double edit_sequence_probability(int dist, int len, double error_rate) {
  if (dist <= 0) return 1.0;
  const double e = clamp_error_rate(error_rate);
  return std::pow(e / 3.0, static_cast<double>(dist)) *
    std::pow(std::max(1e-12, 1.0 - e), static_cast<double>(std::max(0, len - dist)));
}

inline bool local_error_expectation_promotes(int dist, int child_count,
                                             int parent_count, int len,
                                             double error_rate) {
  if (dist <= 1) return false;
  
  const int floor = (dist == 2) ? 5 : 2;
  if (child_count < floor) return false;
  
  const double p_err = edit_sequence_probability(dist, len, error_rate);
  const double expected = std::max(1e-12, static_cast<double>(parent_count) * p_err);
  
  // Conservative local Poisson-style rule, not Shepherd's exact binomial Bayes
  // factor: promote only when the observed child count is far above the number
  // expected from this exact edit path. This catches nearby true barcodes at
  // d=3 without promoting the many singleton errors.
  if (static_cast<double>(child_count) >= std::max(2.0, 25.0 * expected + 2.0)) {
    return true;
  }
  
  return false;
}

inline double log_binomial_pmf(int observed, int trials, double probability) {
  if (observed < 0 || trials < observed || probability <= 0.0 || probability >= 1.0) {
    return -std::numeric_limits<double>::infinity();
  }
  return std::lgamma(static_cast<double>(trials) + 1.0) -
    std::lgamma(static_cast<double>(observed) + 1.0) -
    std::lgamma(static_cast<double>(trials - observed) + 1.0) +
    static_cast<double>(observed) * std::log(probability) +
    static_cast<double>(trials - observed) * std::log1p(-probability);
}

inline bool shepherd_hamming_bayes_promotes(int dist, int child_count,
                                            int parent_count, int len,
                                            int highest_count,
                                            double error_rate) {
  // Shepherd merges every distance-1 neighbor. At d>=2 it merges only when
  // its single-time-point Bayes score exceeds -4; the complementary decision
  // here promotes an already-absorbed sequence back to a centroid. Keeping the
  // exact score makes Hamming refinement faithful to the published algorithm,
  // while barbac's later cleanup can still improve on its one-pass ordering.
  if (dist <= 1) return false;

  const double e = clamp_error_rate(error_rate);
  const double p_no_error = std::pow(1.0 - e, static_cast<double>(len));
  const int inferred_reads = static_cast<int>(static_cast<double>(parent_count) / p_no_error);
  const int trials = std::max(inferred_reads, parent_count + child_count);
  const double p_exact = std::pow(e / 3.0, static_cast<double>(dist)) *
    std::pow(1.0 - e, static_cast<double>(std::max(0, len - dist)));
  const double log_denom = static_cast<double>(len) * std::log(4.0) +
    std::log(static_cast<double>(std::max(1, highest_count)));
  const double log_k = log_binomial_pmf(child_count, trials, p_exact) +
    std::log(p_exact) + log_denom;
  return log_k <= -4.0;
}

void add_cluster(std::vector<Cluster>& clusters,
                 int seq_idx, int seq_len, int count,
                 uint64_t seq_pack, bool seq_packable) {
  Cluster c;
  c.centroid_idx = seq_idx;
  c.centroid_len = seq_len;
  c.centroid_pack = seq_pack;
  c.centroid_packable = seq_packable;
  c.centroid_count = count;
  c.sum_counts = count;
  c.members.push_back(seq_idx);
  c.member_dists.push_back(0);
  clusters.push_back(std::move(c));
}

} // namespace

// =============================================================================
// R exports
// =============================================================================

// [[Rcpp::export]]
std::string barbac_build_id() {
  return std::string(BUILD_ID);
}

// Deterministic salted FNV-1a hash of each sequence, as a non-negative integer.
//
// super_cluster2() uses this to place count-tied barcodes in an order that does
// not depend on their bases. Any tie order is arbitrary, so varying `salt`
// resamples that arbitrary choice and lets a caller measure how much of a
// result rests on it. The value is a pure function of the sequence and the
// salt, so each salt still gives a fully reproducible clustering.
// [[Rcpp::export]]
IntegerVector barbac_seq_order_key(CharacterVector seqs, int salt) {
  const R_xlen_t n = seqs.size();
  IntegerVector out(n);
  const uint64_t basis = 1469598103934665603ULL ^
    (static_cast<uint64_t>(static_cast<uint32_t>(salt)) * 1099511628211ULL);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (seqs[i] == NA_STRING) { out[i] = NA_INTEGER; continue; }
    const char* p = CHAR(STRING_ELT(seqs, i));
    uint64_t h = basis;
    while (*p) {
      h ^= static_cast<unsigned char>(*p++);
      h *= 1099511628211ULL;
    }
    // Fold to 31 bits so the result is a valid non-negative R integer.
    out[i] = static_cast<int>((h ^ (h >> 32)) & 0x7FFFFFFFULL);
  }
  return out;
}

// Order equal-abundance observations using their one-edit error cloud. More
// abundant neighbours are excluded: proximity to a large unrelated lineage
// is not independent evidence that a low-count observation is real.
// [[Rcpp::export]]
NumericVector barbac_support_order_key(CharacterVector seqs, IntegerVector counts,
                                      std::string method) {
  const int n = seqs.size();
  if (counts.size() != n) stop("Sequences and counts must have equal length.");
  const bool lv = parse_lv_method(method);
  std::vector<const char*> ptr(n);
  std::vector<int> len(n);
  std::vector<uint64_t> packed(n);
  std::vector<bool> packable(n);
  std::unordered_map<int, int> count_frequency;
  for (int i = 0; i < n; ++i) {
    if (seqs[i] == NA_STRING) stop("Sequences must not contain NA.");
    ptr[i] = CHAR(STRING_ELT(seqs, i));
    len[i] = std::strlen(ptr[i]);
    packable[i] = pack_seq(ptr[i], len[i], packed[i]);
    ++count_frequency[counts[i]];
  }
  ExactPartitionIndex index(1, lv, ptr, len, counts);
  for (int i = 0; i < n; ++i) index.add(ptr[i], len[i], i);
  NumericVector support(n);
  std::vector<int> candidates;
  for (int i = 0; i < n; ++i) {
    if (count_frequency[counts[i]] < 2) continue;
    if (i % 10000 == 0) checkUserInterrupt();
    index.query(ptr[i], len[i], n, candidates);
    for (int j : candidates) {
      if (j == i || counts[j] > counts[i]) continue;
      int dist = 2;
      if (lv) {
        if (std::abs(len[i]-len[j]) > 1) continue;
        dist = levenshtein_fast(ptr[i], len[i], ptr[j], len[j], 1);
      } else if (len[i] == len[j] && packable[i] && packable[j]) {
        dist = hamming_packed(packed[i], packed[j]);
      }
      if (dist == 1) support[i] += counts[j];
    }
  }
  return support;
}

// Fast abundance-ranked barcode centroid clustering.
 //
 // Internal Rcpp export -- not part of the user-facing R API.
 // Wrapped from R by super_cluster2() in R/11_super_cluster2.R.
 //
 // barcodes:          Character vector of barcode sequences, already sorted
 //                    by descending abundance for best performance/accuracy.
 // counts:            Integer read counts in the same order as `barcodes`.
 // max_distance:      Maximum Hamming or Levenshtein distance.
 // method:            One of "lv", "levenshtein", "hamming", or "ham".
 // kmer_size:         Retained for API compatibility; partitions are learned.
 // min_shared_kmers:  Retained for API compatibility; the Hamming partition
 //                    index computes its own lossless threshold.
 // use_kmer_filter:   If FALSE, scan all centroids; useful for debugging.
 // merge_ratio:       Base abundance ratio for the distance-aware merge guard.
 // error_rate:        Approximate per-base error rate used for parent scoring.
 // verbose:           Print progress and merge-rule diagnostics.
 // returns:           A list of clusters and diagnostic counters.
 // [[Rcpp::export]]
 List barbac_cpp_centroid_cluster_optimized(
     CharacterVector barcodes,
     IntegerVector counts,
     double max_distance,
     std::string method,
     int kmer_size = 5,
     int min_shared_kmers = 2,
     bool use_kmer_filter = true,
     double merge_ratio = 20.0,
     double error_rate = 0.005,
     bool verbose = true,
     bool use_design = false) {
   
   const int n = barcodes.size();
   if (counts.size() != n) stop("`barcodes` and `counts` must have the same length.");
   
   const int D = static_cast<int>(std::floor(max_distance + 1e-9));
   if (D < 0) stop("`max_distance` must be non-negative.");
   const bool is_lv = parse_lv_method(method);
   const double err = clamp_error_rate(error_rate);
   
   if (n == 0) {
     return List::create(
       Named("cluster_id") = CharacterVector(0),
       Named("central_barcode") = CharacterVector(0),
       Named("all_barcodes") = List(0),
       Named("all_counts") = List(0),
       Named("sum_counts") = IntegerVector(0),
       Named("blocked_by_dist") = IntegerVector(D + 1),
       Named("candidate_count") = IntegerVector(0),
       Named("build_id") = barbac_build_id());
   }
   
   std::vector<const char*> ptr(n);
   std::vector<int> slen(n);
   std::vector<uint64_t> packed(n, 0);
   std::vector<bool> packable(n, false);
   std::vector<BaseComposition> composition(n);
   
   for (int i = 0; i < n; ++i) {
     if (barcodes[i] == NA_STRING) stop("`barcodes` contains NA at position %d.", i + 1);
     ptr[i] = CHAR(STRING_ELT(barcodes, i));
     slen[i] = static_cast<int>(std::strlen(ptr[i]));
     uint64_t p = 0;
     packable[i] = pack_seq(ptr[i], slen[i], p);
     if (is_lv && packable[i]) {
       for (int pos = 0; pos < slen[i]; ++pos)
         if (ptr[i][pos] != "ACGT"[dna_base(ptr[i][pos])]) packable[i] = false;
     }
     packed[i] = packable[i] ? p : 0;
     composition[i] = base_composition(ptr[i], slen[i]);
   }

   // Hamming mode assumes one barcode length. Reads carrying an indel break
   // that assumption, and the partition index -- keyed by length -- will never
   // offer them their parent, so each founds a cluster of its own. Rescuing
   // them by edit distance is worth it while they are a trace contaminant, and
   // is the wrong thing to do once they are common: at that point the data
   // wants Levenshtein, and scanning every length from every read would be both
   // slower and less accurate than simply asking for it.
   double off_length_fraction = 0.0;
   int modal_len = 0;
   if (n > 0) {
     std::unordered_map<int, int> len_hist;
     for (int i = 0; i < n; ++i) ++len_hist[slen[i]];
     int modal_count = 0;
     for (std::unordered_map<int, int>::const_iterator it = len_hist.begin();
          it != len_hist.end(); ++it) {
       if (it->second > modal_count) { modal_count = it->second; modal_len = it->first; }
     }
     off_length_fraction = 1.0 - static_cast<double>(modal_count) / static_cast<double>(n);
   }

   // Which positions carry barcode identity, read off the library itself.
   DesignMask design;
   design.usable = false;
   if (use_design) {
     design = learn_design(ptr, slen, counts, n, modal_len);
   }
   if (verbose && use_design) {
     if (design.usable) {
       Rcout << "  Design            : " << design.n_variable << " of "
             << design.len << " positions carry identity ("
             << (design.len - design.n_variable) << " fixed)\n";
     } else {
       Rcout << "  Design            : no fixed positions found; scoring unchanged\n";
     }
   }
   const double off_length_limit = 0.02;
   const bool hamming_rescue_indels = !is_lv && off_length_fraction > 0.0 &&
     off_length_fraction <= off_length_limit;
   if (!is_lv && off_length_fraction > off_length_limit) {
     Rf_warning("%.1f%% of sequences differ from the modal barcode length. Hamming "
                "distance is undefined between sequences of different lengths, so "
                "these cannot be compared and will each form their own cluster. "
                "Use method = \"lv\" for data containing indels.",
                100.0 * off_length_fraction);
   }

   if (verbose) {
     Rcout << "  Build ID          : " << BUILD_ID << "\n";
     Rcout << "  Assignment        : likelihood best-parent + distance-aware merge guard\n";
     Rcout << "  Method            : " << (is_lv ? "Levenshtein" : "Hamming") << "\n";
     if (!is_lv) {
       Rcout << "  Note              : Hamming mode is substitution-only; use LV for indel/shift-sensitive benchmarks\n";
     }
     Rcout << "  Index             : "
           << (use_kmer_filter ? (is_lv ? "exact information-balanced LV partitions" : "exact information-balanced Hamming partitions") : "OFF/full scan")
           << "\n";
     Rcout << "  Error rate        : " << err << "\n";
     Rcout << "  Merge rules (base_ratio=" << merge_ratio << "):\n";
     for (int d = 1; d <= D; ++d) {
       Rcout << "    d=" << d
             << " : floor=" << effective_count_floor(d, is_lv)
             << ", effective_ratio=" << effective_merge_ratio(d, merge_ratio, is_lv)
             << "\n";
     }
     R_FlushConsole();
   }
   
   std::vector<Cluster> clusters;
   clusters.reserve(std::max(16, n / 8));

   // Centroid ids grouped by length, used only by Hamming mode's cross-length
   // rescue below. The Hamming partition index is keyed by sequence length, so
   // it can only ever propose same-length candidates; this is what lets a read
   // carrying an indel still find the parent it came from.
   std::unordered_map<int, std::vector<int> > centroids_by_len;
   // A cross-length comparison is only worth making while they stay rare. On
   // predominantly fixed-length data (what Hamming mode is for) the indel-
   // bearing reads are a fraction of a percent and this costs almost nothing;
   // the budget stops it degenerating into a full scan on data that is really
   // Levenshtein's job, where the user should be in LV mode anyway.
   const long long cross_len_budget = 400LL * 1000LL * 1000LL;
   long long cross_len_comparisons = 0;
   bool cross_len_budget_hit = false;
   
   ExactPartitionIndex index(D, is_lv, ptr, slen, counts, use_kmer_filter && D > 0);
   ExactPartitionIndex lv_near_index(1, true, ptr, slen, counts, use_kmer_filter && is_lv && D > 1);
   
   std::vector<int> candidates;
   std::vector<int> cross_len_candidates;
   candidates.reserve(4096);
   std::vector<int> blocked_by_dist(D + 1, 0);
   std::vector<int> best_by_dist(D + 1, 0);
   long long total_candidates_seen = 0;
   long long lv_verifications = 0;
   long long lv_composition_rejects = 0;
   long long lv_fast_accepts = 0;
   long long hamming_prefilter_rejects = 0;
   long long lv_seed_queries = 0;
   long long cross_len_queries = 0;
   long long design_anchor_only_absorbs = 0;
   long long lv_seed_candidates = 0;
   long long lv_long_seed_queries = 0;
   long long lv_long_seed_candidates = 0;
   long long lv_hamming_stage_assignments = 0;
   long long shepherd_promoted = 0;
   long long shepherd_reassigned = 0;
   long long post_promotion_absorbed = 0;
   int no_candidate_count = 0;
   
   if (verbose) {
     Rcout << "Pass 1: abundance-ranked centroid assignment"
           << " | n=" << n
           << " | d=" << D
           << " | filter=" << (use_kmer_filter ? "ON" : "OFF")
           << std::endl;
     R_FlushConsole();
   }
   
   const auto t0 = std::chrono::steady_clock::now();
   
   for (int i = 0; i < n; ++i) {
     if (i % 10000 == 0) checkUserInterrupt();
     if (verbose && i > 0 && i % 100000 == 0) {
       const auto tn = std::chrono::steady_clock::now();
       const double elapsed = std::chrono::duration<double>(tn - t0).count();
       const double rate = elapsed > 0 ? static_cast<double>(i) / elapsed : 0.0;
       const double eta = rate > 0 ? static_cast<double>(n - i) / rate : 0.0;
       Rcout << "  [" << static_cast<int>(100.0 * i / n) << "%] "
             << i << "/" << n
             << " clusters=" << clusters.size()
             << " no_cand=" << no_candidate_count
             << " lv_dp=" << lv_verifications
             << " lv_fast=" << lv_fast_accepts
             << " lv_seed_q=" << lv_seed_queries
             << " lv_long_q=" << lv_long_seed_queries
             << " blocked[";
       for (int d = 1; d <= D; ++d) {
         if (d > 1) Rcout << "/";
         Rcout << "d" << d << "=" << blocked_by_dist[d];
       }
       Rcout << "] " << static_cast<int>(rate) << " seq/s"
             << " ETA=" << static_cast<int>(eta / 60) << "m" << (static_cast<int>(eta) % 60) << "s"
             << std::endl;
       R_FlushConsole();
     }
     
     const char* s = ptr[i];
     const int sl = slen[i];
     const int cnt = counts[i];
     const bool pk = packable[i];
     const uint64_t pv = packed[i];
     
     if (clusters.empty()) {
       add_cluster(clusters, i, sl, cnt, pv, pk);
       centroids_by_len[sl].push_back(0);
       if (use_kmer_filter && D > 0) {
         index.add(s, sl, 0);
         if (is_lv && D > 1) lv_near_index.add(s, sl, 0);
       }
       continue;
     }
     
     CandidateChoice best_absorb;
     CandidateChoice best_blocked;
     
     auto scan_candidates = [&](const std::vector<int>& scan_set, bool allow_lv_verify) {
       for (int j : scan_set) {
         Cluster& cl = clusters[j];
         if (std::abs(sl - cl.centroid_len) > D) continue;
         
         int dist = D + 1;
         
         if (sl == cl.centroid_len && pk && cl.centroid_packable) {
           if (!is_lv) {
             if (!hamming_bit_prefilter(pv, cl.centroid_pack, D)) {
               ++hamming_prefilter_rejects;
               continue;
             }
             dist = hamming_packed(pv, cl.centroid_pack);
           } else {
             const int ham = hamming_packed(pv, cl.centroid_pack);
             if (ham <= std::min(D, 2)) {
               // This is the common case on fixed-length barcode data and is as
               // cheap as Hamming mode while still running under the LV method.
               dist = ham;
               ++lv_fast_accepts;
             } else if (allow_lv_verify) {
               if (composition_l1(composition[i], composition[cl.centroid_idx]) > 2 * D) {
                 ++lv_composition_rejects;
                 continue;
               }
               ++lv_verifications;
               dist = levenshtein_fast(s, sl, ptr[cl.centroid_idx], cl.centroid_len, D);
             } else {
               continue;
             }
           }
         } else if (allow_lv_verify) {
           // Sequences of different lengths have no Hamming distance, so both
           // modes measure this pair by edit distance. In Hamming mode this is
           // reached only from the cross-length rescue below, where the pair
           // differs in length precisely because one of them carries an indel.
           if (composition_l1(composition[i], composition[cl.centroid_idx]) > 2 * D) {
             ++lv_composition_rejects;
             continue;
           }
           ++lv_verifications;
           dist = levenshtein_fast(s, sl, ptr[cl.centroid_idx], cl.centroid_len, D);
         } else {
           continue;
         }
         
         if (dist <= 0 || dist > D) continue;
         
         const double score = distance_log_likelihood_score(
           dist, cl.centroid_count, cnt, std::max(sl, cl.centroid_len), err, is_lv);
         
         // With a known design, a mismatch on a fixed position is a sequencing
         // error rather than evidence of a second barcode, so only mismatches
         // on identity-carrying positions argue against absorbing. When none
         // of them do, the pair cannot be two different barcodes and the
         // abundance guard has nothing to protect.
         int guard_dist = dist;
         if (design.usable && pk && cl.centroid_packable &&
             sl == design.len && cl.centroid_len == design.len) {
           guard_dist = variable_mismatches(design, pv, cl.centroid_pack);
           if (guard_dist == 0) ++design_anchor_only_absorbs;
         }

         CandidateChoice& target =
           (guard_dist > 0 &&
            merge_blocked_by_ratio(guard_dist, cnt, cl.centroid_count, merge_ratio, is_lv))
           ? best_blocked
         : best_absorb;
         
         if (target.cluster_id < 0 || score > target.score ||
             (score == target.score &&
              (cl.centroid_count > clusters[target.cluster_id].centroid_count ||
               (cl.centroid_count == clusters[target.cluster_id].centroid_count && j < target.cluster_id)))) {
           target.cluster_id = j;
           target.dist = dist;
           target.score = score;
         }
       }
     };
     
     bool used_lv_seed_this_query = false;
     
     if (use_kmer_filter && D > 0) {
       // The Hamming partition is lossless within the configured radius.
       // LV uses its own shift-aware partitions and the score bound below.
       if (!is_lv) {
         index.query(s, sl, static_cast<int>(clusters.size()), candidates);
         total_candidates_seen += static_cast<long long>(candidates.size());
         if (candidates.empty()) ++no_candidate_count;
         scan_candidates(candidates, false);
       }

       // Hamming mode: rescue reads carrying an indel.
       //
       // The partition index is keyed by sequence length, so a read one base
       // shorter than its parent is never even offered as a candidate and ends
       // up founding a cluster of its own -- a false positive for every indel
       // in the data. Same-length candidates were already settled losslessly
       // above, so only the other lengths need looking at, and on the
       // fixed-length data this mode is meant for those centroid lists are
       // nearly empty. Cost is therefore paid per indel-bearing read rather
       // than per read.
       if (hamming_rescue_indels && D > 0 && best_absorb.cluster_id < 0 &&
           !cross_len_budget_hit) {
         cross_len_candidates.clear();
         for (int cand_len = sl - D; cand_len <= sl + D; ++cand_len) {
           if (cand_len == sl || cand_len <= 0) continue;
           std::unordered_map<int, std::vector<int> >::const_iterator it =
             centroids_by_len.find(cand_len);
           if (it == centroids_by_len.end()) continue;
           cross_len_candidates.insert(cross_len_candidates.end(),
                                       it->second.begin(), it->second.end());
         }
         if (!cross_len_candidates.empty()) {
           cross_len_comparisons += static_cast<long long>(cross_len_candidates.size());
           if (cross_len_comparisons > cross_len_budget) {
             cross_len_budget_hit = true;
           } else {
             ++cross_len_queries;
             total_candidates_seen += static_cast<long long>(cross_len_candidates.size());
             scan_candidates(cross_len_candidates, true);
           }
         }
       }

       if (is_lv) {
         // First find every distance-1 parent. Skip the wider search only if
         // even the most abundant possible distance-2 parent cannot beat it.
         // This is a likelihood upper bound, not the old first-hit heuristic.
         if (D > 1) {
           ++lv_long_seed_queries;
           lv_near_index.query(s, sl, static_cast<int>(clusters.size()), candidates);
           lv_long_seed_candidates += static_cast<long long>(candidates.size());
           total_candidates_seen += static_cast<long long>(candidates.size());
           scan_candidates(candidates, true);
         }
         const double unseen_upper = distance_log_likelihood_score(
           2, counts[0], cnt, sl, err, true);
         if (D <= 1 || best_absorb.cluster_id < 0 || best_absorb.score <= unseen_upper) {
           used_lv_seed_this_query = true;
           ++lv_seed_queries;
           index.query(s, sl, static_cast<int>(clusters.size()), candidates);
           lv_seed_candidates += static_cast<long long>(candidates.size());
           total_candidates_seen += static_cast<long long>(candidates.size());
           scan_candidates(candidates, true);
         }
       }
     } else {
       candidates.resize(clusters.size());
       std::iota(candidates.begin(), candidates.end(), 0);
       total_candidates_seen += static_cast<long long>(candidates.size());
       if (candidates.empty()) ++no_candidate_count;
       scan_candidates(candidates, true);
     }
     
     bool assigned = false;
     if (best_absorb.cluster_id >= 0) {
       if (best_absorb.dist >= 0 && best_absorb.dist <= D) ++best_by_dist[best_absorb.dist];
       
       Cluster& parent = clusters[best_absorb.cluster_id];
       parent.members.push_back(i);
       parent.member_dists.push_back(best_absorb.dist);
       parent.sum_counts += cnt;
       assigned = true;
       if (is_lv && !used_lv_seed_this_query) ++lv_hamming_stage_assignments;
     } else if (best_blocked.cluster_id >= 0) {
       if (best_blocked.dist >= 0 && best_blocked.dist <= D) {
         ++best_by_dist[best_blocked.dist];
         ++blocked_by_dist[best_blocked.dist];
       }
     }
     
     if (!assigned) {
       const int new_id = static_cast<int>(clusters.size());
       add_cluster(clusters, i, sl, cnt, pv, pk);
       centroids_by_len[sl].push_back(new_id);
       if (use_kmer_filter && D > 0) {
         index.add(s, sl, new_id);
         if (is_lv && D > 1) lv_near_index.add(s, sl, new_id);
       }
     }
   }
   
   
   // Local statistical refinement: greedy centroid assignment is good at
   // absorbing errors, but it can hide real nearby barcodes inside a larger
   // parent. Conservatively promote only d>=2 children whose counts are far
   // above the exact-edit error expectation, then let those promoted barcodes
   // reclaim better-explained later members from the same local cluster.
   if (D > 1) {
     std::vector<Cluster> refined;
     refined.reserve(clusters.size() + 256);
     std::vector<bool> is_promoted;
     is_promoted.reserve(clusters.size() + 256);
     
     for (const Cluster& cl : clusters) {
       if (cl.members.empty()) continue;
       
       const int root_seq = cl.centroid_idx;
       const int root_cluster = static_cast<int>(refined.size());
       add_cluster(refined, root_seq, cl.centroid_len, cl.centroid_count,
                   cl.centroid_pack, cl.centroid_packable);
       is_promoted.push_back(false);
       
       std::vector<int> promoted_cluster_for_pos(cl.members.size(), -1);
       std::vector<int> promoted_positions;
       for (int pos = 1; pos < static_cast<int>(cl.members.size()); ++pos) {
         const int seq_id = cl.members[pos];
         const int dist = pos < static_cast<int>(cl.member_dists.size()) ? cl.member_dists[pos] : D + 1;
         const int child_count = counts[seq_id];
         const int comparison_len = std::max(slen[seq_id], cl.centroid_len);
         const bool promote = is_lv
           ? local_error_expectation_promotes(dist, child_count, cl.centroid_count,
                                              comparison_len, err)
           : shepherd_hamming_bayes_promotes(dist, child_count, cl.centroid_count,
                                             comparison_len, counts[0], err);
         if (promote) {
           const int new_cluster = static_cast<int>(refined.size());
           add_cluster(refined, seq_id, slen[seq_id], counts[seq_id], packed[seq_id], packable[seq_id]);
           is_promoted.push_back(true);
           promoted_cluster_for_pos[pos] = new_cluster;
           promoted_positions.push_back(pos);
           ++shepherd_promoted;
         }
       }
       
       for (int pos = 1; pos < static_cast<int>(cl.members.size()); ++pos) {
         const int seq_id = cl.members[pos];
         if (promoted_cluster_for_pos[pos] >= 0) continue;
         
         int best_cluster = root_cluster;
         int best_dist = pos < static_cast<int>(cl.member_dists.size()) ? cl.member_dists[pos] : D + 1;
         double best_score = (best_dist > 0 && best_dist <= D)
           ? distance_log_likelihood_score(
               best_dist, refined[root_cluster].centroid_count, counts[seq_id],
                                                                      std::max(slen[seq_id], refined[root_cluster].centroid_len), err, is_lv)
             : -std::numeric_limits<double>::infinity();
         
         for (int ppos : promoted_positions) {
           const int promoted_cluster = promoted_cluster_for_pos[ppos];
           if (promoted_cluster < 0) continue;
           const int promoted_seq = cl.members[ppos];
           if (seq_id == promoted_seq) continue;
           if (std::abs(slen[seq_id] - slen[promoted_seq]) > D) continue;
           
           int d_new = D + 1;
           if (slen[seq_id] == slen[promoted_seq] && packable[seq_id] && packable[promoted_seq]) {
             const int ham = hamming_packed(packed[seq_id], packed[promoted_seq]);
             if (!is_lv) {
               d_new = (ham <= D) ? ham : D + 1;
             } else if (ham <= std::min(D, 2)) {
               d_new = ham;
             } else {
               d_new = levenshtein_fast(ptr[seq_id], slen[seq_id], ptr[promoted_seq], slen[promoted_seq], D);
             }
           } else if (is_lv) {
             d_new = levenshtein_fast(ptr[seq_id], slen[seq_id], ptr[promoted_seq], slen[promoted_seq], D);
           }
           
           if (d_new <= 0 || d_new > D) continue;
           const double new_score = distance_log_likelihood_score(
             d_new, refined[promoted_cluster].centroid_count, counts[seq_id],
                                                                    std::max(slen[seq_id], refined[promoted_cluster].centroid_len), err, is_lv);
           if (new_score > best_score) {
             best_score = new_score;
             best_dist = d_new;
             best_cluster = promoted_cluster;
           }
         }
         
         refined[best_cluster].members.push_back(seq_id);
         refined[best_cluster].member_dists.push_back(best_dist);
         refined[best_cluster].sum_counts += counts[seq_id];
         if (best_cluster != root_cluster) ++shepherd_reassigned;
       }
     }

     // A promoted barcode was not available as a parent during pass 1. Its
     // low-count error variants may therefore already have founded separate
     // clusters. Revisit only those pre-existing roots and only against newly
     // promoted Hamming centroids. Shepherd's conservative unconditional cases
     // are used here: distance 1, or a singleton within the requested radius.
     // This closes the ordering hole without re-merging protected d>=2
     // multi-read barcodes or changing LV behavior.
     if (!is_lv && shepherd_promoted > 0) {
       std::vector<int> promoted_ids;
       promoted_ids.reserve(static_cast<size_t>(shepherd_promoted));
       for (int cid = 0; cid < static_cast<int>(refined.size()); ++cid) {
         if (is_promoted[cid]) promoted_ids.push_back(cid);
       }

       std::vector<int> absorb_into(refined.size(), -1);
       for (int cid = 0; cid < static_cast<int>(refined.size()); ++cid) {
         if (is_promoted[cid]) continue;
         const Cluster& child = refined[cid];
         int best_parent = -1;
         int best_dist = D + 1;
         int best_parent_count = -1;
         for (int pid : promoted_ids) {
           const Cluster& parent = refined[pid];
           if (child.centroid_len != parent.centroid_len ||
               !child.centroid_packable || !parent.centroid_packable) continue;
           if (!hamming_bit_prefilter(child.centroid_pack, parent.centroid_pack, D)) continue;
           const int dist = hamming_packed(child.centroid_pack, parent.centroid_pack);
           if (dist <= 0 || dist > D) continue;
           if (dist != 1 && child.centroid_count != 1) continue;
           if (dist < best_dist ||
               (dist == best_dist && parent.centroid_count > best_parent_count)) {
             best_parent = pid;
             best_dist = dist;
             best_parent_count = parent.centroid_count;
           }
         }
         absorb_into[cid] = best_parent;
       }

       for (int cid = 0; cid < static_cast<int>(refined.size()); ++cid) {
         const int pid = absorb_into[cid];
         if (pid < 0) continue;
         Cluster& parent = refined[pid];
         Cluster& child = refined[cid];
         parent.members.insert(parent.members.end(), child.members.begin(), child.members.end());
         parent.member_dists.insert(parent.member_dists.end(),
                                    child.member_dists.begin(), child.member_dists.end());
         parent.sum_counts += child.sum_counts;
         ++post_promotion_absorbed;
       }

       if (post_promotion_absorbed > 0) {
         std::vector<Cluster> cleaned;
         cleaned.reserve(refined.size() - static_cast<size_t>(post_promotion_absorbed));
         for (int cid = 0; cid < static_cast<int>(refined.size()); ++cid) {
           if (absorb_into[cid] < 0) cleaned.push_back(std::move(refined[cid]));
         }
         refined.swap(cleaned);
       }
     }
     
     clusters.swap(refined);
   }
   
   const auto t1 = std::chrono::steady_clock::now();
   const double elapsed = std::chrono::duration<double>(t1 - t0).count();
   
   if (verbose) {
     int total_blocked = 0;
     Rcout << "  [100%] done in "
           << static_cast<int>(elapsed / 60) << "m" << (static_cast<int>(elapsed) % 60) << "s"
           << " clusters=" << clusters.size()
           << " no_cand=" << no_candidate_count
           << " avg_cand=" << (n > 0 ? static_cast<double>(total_candidates_seen) / n : 0.0)
           << " lv_dp=" << lv_verifications
           << " lv_comp_reject=" << lv_composition_rejects
           << " lv_fast=" << lv_fast_accepts
           << " lv_seed_q=" << lv_seed_queries
           << " lv_seed_cand=" << lv_seed_candidates
           << " lv_long_q=" << lv_long_seed_queries
           << " lv_long_cand=" << lv_long_seed_candidates
           << " promoted=" << shepherd_promoted
           << " reassigned=" << shepherd_reassigned
           << " post_promote_absorb=" << post_promotion_absorbed
           << " ham_reject=" << hamming_prefilter_rejects
           << " blocked[";
     for (int d = 1; d <= D; ++d) {
       if (d > 1) Rcout << "/";
       Rcout << "d" << d << "=" << blocked_by_dist[d];
       total_blocked += blocked_by_dist[d];
     }
     Rcout << "] total=" << total_blocked
           << " " << (elapsed > 0 ? static_cast<int>(n / elapsed) : 0) << " seq/s"
           << std::endl;
     R_FlushConsole();
   }
   
   const int nc = static_cast<int>(clusters.size());
   CharacterVector r_cluster_id(nc);
   CharacterVector r_central(nc);
   List r_all_bc(nc);
   List r_all_cnt(nc);
   IntegerVector r_sum(nc);
   
   for (int c = 0; c < nc; ++c) {
     const Cluster& cl = clusters[c];
     r_cluster_id[c] = "group" + std::to_string(c + 1);
     r_central[c] = barcodes[cl.centroid_idx];
     r_sum[c] = cl.sum_counts;
     
     const int m = static_cast<int>(cl.members.size());
     CharacterVector bc_vec(m);
     IntegerVector ct_vec(m);
     for (int j = 0; j < m; ++j) {
       const int seq_id = cl.members[j];
       bc_vec[j] = barcodes[seq_id];
       ct_vec[j] = counts[seq_id];
     }
     r_all_bc[c] = bc_vec;
     r_all_cnt[c] = ct_vec;
   }
   
   IntegerVector r_blocked(D + 1);
   IntegerVector r_best(D + 1);
   for (int d = 0; d <= D; ++d) {
     r_blocked[d] = blocked_by_dist[d];
     r_best[d] = best_by_dist[d];
   }
   
   return List::create(
     Named("cluster_id") = r_cluster_id,
     Named("central_barcode") = r_central,
     Named("all_barcodes") = r_all_bc,
     Named("all_counts") = r_all_cnt,
     Named("sum_counts") = r_sum,
     Named("blocked_by_dist") = r_blocked,
     Named("best_match_by_dist") = r_best,
     Named("candidate_count") = NumericVector::create(
       Named("total") = static_cast<double>(total_candidates_seen),
       Named("average") = n > 0 ? static_cast<double>(total_candidates_seen) / n : 0.0,
       Named("no_candidate") = no_candidate_count),
       Named("distance_count") = NumericVector::create(
         Named("lv_verifications") = static_cast<double>(lv_verifications),
         Named("lv_composition_rejects") = static_cast<double>(lv_composition_rejects),
         Named("lv_fast_accepts") = static_cast<double>(lv_fast_accepts),
         Named("lv_seed_queries") = static_cast<double>(lv_seed_queries),
         Named("lv_seed_candidates") = static_cast<double>(lv_seed_candidates),
         Named("lv_long_seed_queries") = static_cast<double>(lv_long_seed_queries),
         Named("lv_long_seed_candidates") = static_cast<double>(lv_long_seed_candidates),
         Named("lv_hamming_stage_assignments") = static_cast<double>(lv_hamming_stage_assignments),
         Named("hamming_prefilter_rejects") = static_cast<double>(hamming_prefilter_rejects),
         Named("cross_length_queries") = static_cast<double>(cross_len_queries)),
         Named("refinement_count") = NumericVector::create(
           Named("shepherd_promoted") = static_cast<double>(shepherd_promoted),
           Named("shepherd_reassigned") = static_cast<double>(shepherd_reassigned),
           Named("post_promotion_absorbed") = static_cast<double>(post_promotion_absorbed)),
           Named("method") = is_lv ? "levenshtein" : "hamming",
           Named("build_id") = barbac_build_id());
 }
