#ifndef LOG_SKETCH_H
#define LOG_SKETCH_H
#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <hash.h>
#include <Parameter.h>
#define EMA
#define Stb
//#define decay_EMA

struct LogBucket {
    float center = 0.0f;
#ifdef Stb
    float stable = 0;
#endif
    uint16_t total_cnt = 0;
    uint16_t matching_cnt = 0;

#ifndef EMA
    float accumulate = 0.0f;
    uint16_t acc_cnt = 0;
#else
    //float stable = 0.0f;
#endif

    void update(uint16_t latency) {
        total_cnt++;

        if (center == 0) {
            center = static_cast<float>(latency);
        }

        float diff = std::abs(static_cast<float>(latency) - center);
        float times = diff / center;
        if ((center > 0 && (times) < alpha) || diff < beta) {
            ++matching_cnt;
        }
#ifdef EMA

        if (times > thresh) return;
#ifdef Stb
        if (times <= alpha * 2 )
        {
            stable += alpha * 2/(1-alpha * 2);
        }else
        {
            stable -= 1;
        }
        if (total_cnt % 100 ==0)
        {
            //if (stable < 0)
            {
                float new_cnt = std::max(0.0,(float)(matching_cnt)+  stable * alpha);
                matching_cnt = uint32_t(new_cnt);
            }
#ifndef decay_EMA
            stable = 0;
#endif
        }
#ifdef decay_EMA
        stable *=0.98;
#endif

        if ( stable > 0)
        {
            return;
        }
#endif

        float epsilon =alpha_ema;


        center = (1.0f - epsilon) * center + epsilon * static_cast<float>(latency);
#else
        accumulate += (static_cast<float>(latency) - center);
        acc_cnt++;

        if (acc_cnt >= 32 || std::abs(accumulate) > 300000) {
            auto fix = static_cast<float>(accumulate / acc_cnt);
            center += fix;
            if (center < 1) center = 1;
            accumulate = 0.0f;
            acc_cnt = 0;
        } // better accuracy but memory inefficient
#endif

    }
};

class LogSketch {
public:
    std::array<LogBucket,MAX_BIN> buckets;

    explicit LogSketch()  {
    }

    void update(uint16_t latency) {
        int idx = 0;
        if (latency > 0) {
            idx = std:: min(MAX_BIN-1,
                static_cast<int>(std::log(latency) / std::log(DEFAULT_BASE)));
        }
        idx = std::min(idx, (int)buckets.size() - 1);
        buckets[idx].update(latency);
    }
};

struct OutputResult2 {
    uint32_t total_pkts;
    uint32_t total_matching;

    struct BucketDetail {
        int bucket_id;
        uint32_t b_total;
        uint32_t b_matching;
        float b_center;
    };
    std::vector<BucketDetail> top_k_buckets;
};

struct Cell_sp {
    uint32_t id;
    uint32_t cnt;
    uint32_t t_x;
    LogSketch sketch;

    Cell_sp() : id(0), cnt(0), t_x(0), sketch() {}

    void reset(uint32_t new_id, uint32_t timestamp) {
        id = new_id;
        cnt = 1;
        t_x = timestamp;
        sketch = LogSketch();
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
        : sg_(sampling_rate), m_(m), n_(n), theta_(theta), cnt_g_(0) {
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
                uint16_t latency = ((double)(value - bkt.Cell_sps[i].t_x) / sg_);
                bkt.Cell_sps[i].t_x = value;
                bkt.Cell_sps[i].sketch.update(latency);
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
                        all_bkts.push_back({b_idx, b.total_cnt, b.matching_cnt, b.center});
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
#endif //LOG_SKETCH_H