#ifndef LOG_SKETCH_H
#define LOG_SKETCH_H

#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <hash.h>
#include <Parameter.h>
#include <cstdlib>

#define EMA
//#define SIMD
#ifdef SIMD
#include <immintrin.h>
#endif

struct LogBucket {
    uint32_t total_cnt = 0;
    uint32_t matching_cnt = 0;
    float center = 0.0f;
    float sum_sq_diff = 0.0f;

#ifdef Stb
    float stable = 0;
#endif


    inline void update(uint16_t latency) {
        total_cnt++;

        if (center == 0) {
            center = static_cast<float>(latency);
        }

        float diff = std::abs(static_cast<float>(latency) - center);
        float times = diff / center;
        if ((center > 0 && (times) < alpha) || diff < beta) {
            matching_cnt++;

#ifdef Stb
            stable += stb_para / (1 - stb_para);
#endif

        }else
        {
#ifdef Stb
            stable -= 1;
#endif
        }
#ifdef Stb
        if (total_cnt % 100 == 0) {
           // if (stable < 0) {
                float new_cnt = std::max(0.0f, (float)((matching_cnt) + alpha * stable));
                matching_cnt = uint32_t(new_cnt);
          //  }
            stable *= 0.1;
        }
#endif

        float d = static_cast<float>(latency) - center;
        center += (static_cast<float>(latency) - center ) *( alpha_ema/(alpha_ema + total_cnt));
        {
            d *=  static_cast<float>(latency) - center;
            sum_sq_diff += d;
        }
    }
};

class LogSketch {
private:
    uint32_t current_bins_ = 0;
    uint8_t bucket_ids_[MAX_BIN];
    LogBucket real_buckets_[MAX_BIN];

#ifdef HOT_CACHE
    // Single-entry last-hit cache
    // 0xFF = invalid; stores {bin_idx, bucket_idx} of the last exact match.
    uint8_t last_hit_bin_    = 0xFF;
    uint8_t last_hit_bucket_ = 0;
#endif


public:
    double base_a;

    struct LogBucketProxy {
        uint32_t total_cnt;
        uint32_t matching_cnt;
        float center;
        float sum_sq_diff;
        uint8_t bucket_id;
    };

    struct BucketProxy {
        LogSketch* self;

        size_t size() const { return self->current_bins_; }

        LogBucketProxy& operator[](size_t i) {
            thread_local LogBucketProxy p;
            p.total_cnt = self->real_buckets_[i].total_cnt;
            p.matching_cnt = self->real_buckets_[i].matching_cnt;
            p.center = self->real_buckets_[i].center;
            p.sum_sq_diff = self->real_buckets_[i].sum_sq_diff;
            p.bucket_id = self->bucket_ids_[i];
            return p;
        }

        const LogBucketProxy& operator[](size_t i) const {
            thread_local LogBucketProxy p;
            p.total_cnt = self->real_buckets_[i].total_cnt;
            p.matching_cnt = self->real_buckets_[i].matching_cnt;
            p.center = self->real_buckets_[i].center;
            p.sum_sq_diff = self->real_buckets_[i].sum_sq_diff;
            p.bucket_id = self->bucket_ids_[i];
            return p;
        }
    };

    BucketProxy buckets;

    explicit LogSketch(double a = DEFAULT_BASE) : current_bins_(0), base_a(a) {
        buckets.self = this;
    }

    LogSketch(const LogSketch& other) : current_bins_(other.current_bins_), base_a(other.base_a)
#ifdef HOT_CACHE
                                        ,last_hit_bin_(other.last_hit_bin_), last_hit_bucket_(other.last_hit_bucket_)
#endif
    {
        buckets.self = this;
        for (uint32_t i = 0; i < current_bins_; ++i) {
            real_buckets_[i] = other.real_buckets_[i];
            bucket_ids_[i] = other.bucket_ids_[i];
        }
    }
    LogSketch& operator=(const LogSketch& other) {
        if (this != &other) {
            base_a = other.base_a;
            current_bins_ = other.current_bins_;
#ifdef HOT_CACHE
            last_hit_bin_    = other.last_hit_bin_;
            last_hit_bucket_ = other.last_hit_bucket_;
#endif

            buckets.self = this;
            for (uint32_t i = 0; i < current_bins_; ++i) {
                real_buckets_[i] = other.real_buckets_[i];
                bucket_ids_[i] = other.bucket_ids_[i];
            }
        }
        return *this;
    }

    inline void update(uint16_t latency, int idx) {

#ifdef SIMD
        if (current_bins_ > 0) {
            // idx → mm256
            __m256i v_target = _mm256_set1_epi8(target_idx);
            //load index
            __m256i v_ids = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(bucket_ids_));
            //compare
            __m256i v_cmp = _mm256_cmpeq_epi8(v_ids, v_target);

            uint32_t mask = _mm256_movemask_epi8(v_cmp);

            uint32_t valid_mask = (1U << current_bins_) - 1;
            mask &= valid_mask;
            // found
            if (mask != 0) {
                match_idx = __builtin_ctz(mask);
                found = true;
            }
        }
        //not found
        if (!found && current_bins_ == MAX_BIN) {
            for (uint32_t i = 0; i < current_bins_; ++i) {
                int d = std::abs(static_cast<int>(bucket_ids_[i]) - (idx));
                if (d < min_diff) {
                    min_diff = d;
                    diff_index = i;
                }
            }
        }

#else

#ifdef  HOT_CACHE
        // hot_cache
        if (last_hit_bin_ == idx) {
            real_buckets_[last_hit_bucket_].update(latency);
            return;
        }
#endif

        int min_diff = INFINITY;
        uint16_t diff_index = 0;

        for (uint32_t i = 0; i < current_bins_; ++i) {
            int d = std::abs(static_cast<int>(bucket_ids_[i]) - (idx));
            if (d == 0) {
                real_buckets_[i].update(latency);
#ifdef HOT_CACHE
                last_hit_bin_    = static_cast<uint8_t>(idx);
                last_hit_bucket_ = static_cast<uint8_t>(i);
#endif

                return;
            }
            if (d<min_diff) {
                min_diff = d;
                diff_index = i;
            }
        }

        if (current_bins_ < MAX_BIN) {
            uint32_t i = current_bins_++;
            real_buckets_[i] = LogBucket();
            bucket_ids_[i] = static_cast<uint8_t>(idx);
            real_buckets_[i].update(latency);
            return;
        }

        uint32_t min_idx = 0;
        uint32_t min_cnt = real_buckets_[0].total_cnt;
        for (uint32_t i = 1; i < MAX_BIN; ++i) {
            if (real_buckets_[i].total_cnt < min_cnt) {
                min_cnt = real_buckets_[i].total_cnt;
                min_idx = i;
            }
        }

        if (rand() % (min_cnt*5 + 10) == 0) {
            real_buckets_[min_idx] = LogBucket();
            bucket_ids_[min_idx] = static_cast<uint8_t>(idx);
            real_buckets_[min_idx].update(latency);
        } else {
            real_buckets_[diff_index].total_cnt++;

            float d = static_cast<float>(latency) - real_buckets_[diff_index].center;
            real_buckets_[diff_index].sum_sq_diff += d * d;
        }

#endif
    }


private:
};

struct OutputResult2 {
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

struct Cell_sp {
    uint32_t id;
    uint32_t cnt;
    uint32_t t_x;
    LogSketch sketch;

    Cell_sp() : id(0), cnt(0), t_x(0), sketch(DEFAULT_BASE) {}

    void reset(uint32_t new_id, uint32_t timestamp) {
        id = new_id;
        cnt = 1;
        t_x = timestamp;
        sketch = LogSketch(DEFAULT_BASE);
    }
};

struct Bucket2 {
    std::vector<Cell_sp> Cell_sps;
    uint16_t cnt_n;
    explicit Bucket2(int n) : Cell_sps(n), cnt_n(0) {}
};

class L1 {
public:
    explicit L1(uint16_t sampling_rate, uint16_t m, uint8_t n, float theta)
        : sg_(sampling_rate), m_(m), n_(n), theta_(theta), cnt_g_(0)
    {
        buckets_.resize(m_, Bucket2(n_));
    }

    bool should_sample() {
        return (++cnt_g_ % sg_ == 0);
    }

    void insert(uint32_t id, uint32_t value) {
        if (!should_sample()) return;

        uint32_t idx = hash_func(id) % m_;
        Bucket2& bkt = buckets_[idx];

        int min_idx = 0;
        uint16_t min_cnt = UINT16_MAX;
        bool found = false;

        for (int i = 0; i < n_; ++i) {
            if (bkt.Cell_sps[i].id == id && bkt.Cell_sps[i].cnt > 0) {
                bkt.Cell_sps[i].cnt++;
                uint16_t latency = std::round((double)(value - bkt.Cell_sps[i].t_x) / sg_);
                bkt.Cell_sps[i].t_x = value;

                int index = 0;
                if (latency > 0) {
                    index = static_cast<int>(std::log(latency) / std::log(DEFAULT_BASE)) + 1;
                }
                bkt.Cell_sps[i].sketch.update(latency, index);
                found = true;
                break;
            }
            if (bkt.Cell_sps[i].cnt == 0) {
                bkt.Cell_sps[i].reset(id, value);
                return;
            }
            if (bkt.Cell_sps[i].cnt < min_cnt) {
                min_cnt = bkt.Cell_sps[i].cnt;
                min_idx = i;
            }
        }

        if (found) return;

        bkt.cnt_n++;
        if ((double)bkt.cnt_n * theta_ > (double)min_cnt) {
            bkt.Cell_sps[min_idx].reset(id, value);
            bkt.cnt_n = 0;
        }
    }

    std::unordered_map<uint32_t, OutputResult2> output_top_k(int k) {
        std::unordered_map<uint32_t, OutputResult2> results;
        for (int i = 0; i < m_; ++i) {
            for (auto& Cell_sp : buckets_[i].Cell_sps) {
                if (Cell_sp.cnt == 0) continue;

                OutputResult2 res;
                res.total_pkts = 0;
                res.total_matching = 0;

                std::vector<OutputResult2::BucketDetail> all_bkts;
                for (int b_idx = 0; b_idx < (int)Cell_sp.sketch.buckets.size(); ++b_idx) {
                    auto& b = Cell_sp.sketch.buckets[b_idx];
                    res.total_pkts += b.total_cnt;
                    res.total_matching += b.matching_cnt;
                    if (b.total_cnt > 0) {
                        all_bkts.push_back({b.bucket_id, b.total_cnt, b.matching_cnt, b.center, b.sum_sq_diff});
                    }
                }

                std::sort(all_bkts.begin(), all_bkts.end(),
                    [](const OutputResult2::BucketDetail& a, const OutputResult2::BucketDetail& b) {
                    return a.b_total > b.b_total;
                });

                int limit = std::min((int)all_bkts.size(), k);
                for (int j = 0; j < limit; ++j) {
                    res.top_k_buckets.push_back(all_bkts[j]);
                }
                results[Cell_sp.id] = res;
            }
        }
        return results;
    }

private:
    uint16_t sg_, m_;
    uint8_t n_;
    float theta_;
    uint32_t cnt_g_;
    std::vector<Bucket2> buckets_;
    static uint32_t hash_func(uint32_t id) {
        return Hash::BOBHash32((uint8_t*)&id, sizeof(uint32_t), SEED);
    }
};

#endif//LOG_SKETCH_H
