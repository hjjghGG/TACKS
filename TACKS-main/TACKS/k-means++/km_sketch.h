#ifndef KM_SKETCH_H
#define KM_SKETCH_H

#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <hash.h>
#include <Parameter.h>

struct OutputResult5 {
    uint32_t total_pkts;
    uint32_t total_matching;

    struct BucketDetail {
        int bucket_id;
        uint32_t b_total;
        uint32_t b_matching;
        float b_center;
        float b_sum_sq;
    };
    std::vector<BucketDetail> top_k_buckets;
};
struct KMBucket {
    uint32_t total_cnt = 0;
    uint32_t matching_cnt = 0;
    float center = 0.0f;
    float sum_sq_diff = 0.0f;

    void update(uint16_t latency, bool is_match) {
        total_cnt++;
        if (is_match) {
            matching_cnt++;
        }

        float old_center = center;
        center += (static_cast<float>(latency) - center) / total_cnt;
        sum_sq_diff += (static_cast<float>(latency) - old_center) * (static_cast<float>(latency) - center);
    }
};

class KMSketch {
public:
    std::vector<KMBucket> buckets;

    KMSketch() {
        buckets.reserve(S_NUM);
    }

    void update(uint16_t latency) {
        if (buckets.empty()) {
            add_bucket(latency);
            return;
        }

        int best_idx = -1;
        float min_diff = 1e9;
        for (int i = 0; i < (int)buckets.size(); ++i) {
            float diff = std::abs(static_cast<float>(latency) - buckets[i].center);
            if (diff < min_diff) {
                min_diff = diff;
                best_idx = i;
            }
        }

        bool is_match = (buckets[best_idx].center > 0 && (min_diff / buckets[best_idx].center) < alpha) || min_diff < beta;

        if (!is_match && buckets.size() < S_NUM) {
            add_bucket(latency);
            return;
        }

        buckets[best_idx].update(latency, is_match);
    }

private:
    void add_bucket(uint16_t latency) {
        KMBucket b;
        b.center = static_cast<float>(latency);
        b.total_cnt = 1;
        b.matching_cnt = 1;
        buckets.push_back(b);
    }
};

struct Cell_km {
    uint32_t id;
    uint32_t cnt;
    uint32_t t_x;
    KMSketch sketch;

    Cell_km() : id(0), cnt(0), t_x(0) {}

    void reset(uint32_t new_id, uint32_t timestamp) {
        id = new_id;
        cnt = 1;
        t_x = timestamp;
        sketch = KMSketch();
    }
};

struct Bucket4 {
    std::vector<Cell_km> cells;
    uint16_t cnt_n;
    explicit Bucket4(int n) : cells(n), cnt_n(0) {}
};

class L8 {
public:
    explicit L8(uint16_t sampling_rate, uint16_t m, uint8_t n, float theta)
        : sg_(sampling_rate), m_(m), n_(n), theta_(theta), cnt_g_(0) {
        buckets_.resize(m_, Bucket4(n_));
    }

    void insert(uint32_t id, uint32_t value) {
        if (++cnt_g_ % sg_ != 0) return;

        uint32_t idx = Hash::BOBHash32((uint8_t*)&id, sizeof(uint32_t), SEED) % m_;
        Bucket4& bkt = buckets_[idx];

        int min_idx = 0;
        uint16_t min_cnt = UINT16_MAX;
        for (int i = 0; i < n_; ++i) {
            if (bkt.cells[i].id == id && bkt.cells[i].cnt > 0) {
                bkt.cells[i].cnt++;
                uint16_t latency = std::round((double)(value - bkt.cells[i].t_x) / sg_);
                bkt.cells[i].t_x = value;
                bkt.cells[i].sketch.update(latency);
                return;
            }
            if (bkt.cells[i].cnt == 0) {
                bkt.cells[i].reset(id, value);
                return;
            }
            if (bkt.cells[i].cnt < min_cnt) {
                min_cnt = bkt.cells[i].cnt;
                min_idx = i;
            }
        }

        bkt.cnt_n++;
        if ((double)bkt.cnt_n * theta_ > (double)min_cnt) {
            bkt.cells[min_idx].reset(id, value);
            bkt.cnt_n = 0;
        }
    }

    std::unordered_map<uint32_t, OutputResult5> output_top_k(int k) {
        std::unordered_map<uint32_t, OutputResult5> results;
        for (auto& bkt : buckets_) {
            for (auto& cell : bkt.cells) {
                if (cell.cnt == 0) continue;
                OutputResult5 res = {0, 0, {}};
                for (int i = 0; i < (int)cell.sketch.buckets.size(); ++i) {
                    auto& b = cell.sketch.buckets[i];
                    res.total_pkts += b.total_cnt;
                    res.total_matching += b.matching_cnt;
                    res.top_k_buckets.push_back({i, b.total_cnt, b.matching_cnt, b.center, b.sum_sq_diff});
                }
                std::sort(res.top_k_buckets.begin(), res.top_k_buckets.end(),
                    [](const auto& a, const auto& b){ return a.b_total > b.b_total; });

                int limit = std::min((int)res.top_k_buckets.size(), k);
                for (int j = 0; j < limit; ++j) {
                    res.top_k_buckets.push_back(res.top_k_buckets[j]);
                }
                results[cell.id] = res;
            }
        }
        return results;
    }

private:
    uint16_t sg_, m_;
    uint8_t n_;
    float theta_;
    uint32_t cnt_g_;
    std::vector<Bucket4> buckets_;
};

#endif