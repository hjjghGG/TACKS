#ifndef DSSKETCH_H
#define DSSKETCH_H

#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include "hash.h"
#include "Parameter.h"

// D-Stream
struct DSBin {
    uint32_t total_cnt = 0;
    uint32_t matching_cnt = 0;
    float sum_sq_diff = 0.0f;

    float D = 0.0;       // D
    float MP_1 = 0.0;    // D1
    uint32_t last_tg = 0; // T_x

    float get_center() const {
        return D > 0 ? static_cast<float>(MP_1 / D) : 0.0f;
    }

    // Decay
    void decay(uint32_t current_tg, float lambda) {
        if (current_tg > last_tg) {
            float factor = std::pow(lambda, static_cast<float>(current_tg - last_tg));
            D *= factor;
            MP_1 *= factor;
            last_tg = current_tg;
        }
    }
};

class DSSketch {
public:
    std::unordered_map<int, DSBin> grid;

    // D-Stream param
    int W_ = 100;               // W_
    float lambda_ = 0.99;    // Decay
    uint32_t T_gap_ = 1000;   // Time window size
    float D_L_ = 0.1;        // noise thresh
    uint32_t B_ = D_NUM;      // Max b

    DSSketch() {}

    void update(uint16_t latency, uint32_t timestamp) {
        uint32_t tg_curr = timestamp / T_gap_;

        int best_idx = -1;
        float min_diff = 1e9;
        for (auto& pair : grid) {
            float diff = std::abs(static_cast<float>(latency) - pair.second.get_center());
            if (diff < min_diff) {
                min_diff = diff;
                best_idx = pair.first;
            }
        }

        bool is_match = false;
        if (best_idx != -1) {
            float center = grid[best_idx].get_center();
            is_match = (center > 0 && (min_diff / center) < alpha) || min_diff < beta;
        }

        // map
        int target_idx = latency / W_;

        // check B_
        if (grid.find(target_idx) == grid.end()) {
            if (grid.size() >= B_) {
                cleanup_and_cluster(tg_curr);
                if (grid.size() >= B_) {
                    return; 
                }
            }
            grid[target_idx] = DSBin();
            grid[target_idx].last_tg = tg_curr;
        }

        grid[target_idx].decay(tg_curr, lambda_);
        grid[target_idx].D += 1.0;
        grid[target_idx].MP_1 += latency;
        grid[target_idx].total_cnt++;

        {
            float new_center = grid[target_idx].get_center();
            float d = static_cast<float>(latency) - new_center;
            grid[target_idx].sum_sq_diff += d * d;
        }
        
        if (is_match) {
            grid[target_idx].matching_cnt++;
        }
    }

private:
    void cleanup_and_cluster(uint32_t tg_curr) {
        std::vector<int> to_remove;
        for (auto& pair : grid) {
            pair.second.decay(tg_curr, lambda_);
            if (pair.second.D < D_L_) {
                to_remove.push_back(pair.first);
            }
        }
        for (int idx : to_remove) {
            grid.erase(idx);
        }

        while (grid.size() >= B_) {
            int best_k1 = -1, best_k2 = -1;
            float min_dist = 1e9;
            
            for (auto it1 = grid.begin(); it1 != grid.end(); ++it1) {
                auto it2 = it1; ++it2;
                for (; it2 != grid.end(); ++it2) {
                    float dist = std::abs(it1->second.get_center() - it2->second.get_center());
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_k1 = it1->first;
                        best_k2 = it2->first;
                    }
                }
            }
            
            if (best_k1 != -1 && best_k2 != -1) {
                grid[best_k1].D += grid[best_k2].D;
                grid[best_k1].MP_1 += grid[best_k2].MP_1;
                grid[best_k1].total_cnt += grid[best_k2].total_cnt;
                grid[best_k1].matching_cnt += grid[best_k2].matching_cnt;
                grid[best_k1].sum_sq_diff += grid[best_k2].sum_sq_diff;
                grid[best_k1].last_tg = std::max(grid[best_k1].last_tg, grid[best_k2].last_tg);
                
                grid.erase(best_k2);
            } else {
                break; 
            }
        }
    }
};

struct Cell_DS {
    uint32_t id;
    uint32_t cnt;
    uint32_t t_x;
    DSSketch sketch;

    Cell_DS() : id(0), cnt(0), t_x(0) {}

    void reset(uint32_t new_id, uint32_t timestamp) {
        id = new_id;
        cnt = 1;
        t_x = timestamp;
        sketch = DSSketch();
    }
};

struct Bucket3_DS {
    std::vector<Cell_DS> cells;
    uint16_t cnt_n;
    explicit Bucket3_DS(int n) : cells(n), cnt_n(0) {}
};

class L2_DS {
public:

    explicit L2_DS(uint16_t sampling_rate, uint16_t m, uint8_t n, float theta)
        : sg_(sampling_rate), m_(m), n_(n), theta_(theta), cnt_g_(0) {
        buckets_.resize(m_, Bucket3_DS(n_));
    }

    void insert(uint32_t id, uint32_t value) {
        if (++cnt_g_ % sg_ != 0) return;

        uint32_t idx = Hash::BOBHash32((uint8_t*)&id, sizeof(uint32_t), SEED) % m_;
        Bucket3_DS& bkt = buckets_[idx];

        int min_idx = 0;
        uint16_t min_cnt = UINT16_MAX;
        for (int i = 0; i < n_; ++i) {
            if (bkt.cells[i].id == id && bkt.cells[i].cnt > 0) {
                bkt.cells[i].cnt++;
                uint16_t latency = std::round((float)(value - bkt.cells[i].t_x) / sg_);
                bkt.cells[i].t_x = value;
                bkt.cells[i].sketch.update(latency, value);
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
        if ((float)bkt.cnt_n * theta_ > (float)min_cnt) {
            bkt.cells[min_idx].reset(id, value);
            bkt.cnt_n = 0;
        }
    }

    std::unordered_map<uint32_t, OutputResult2> output_top_k(int k) {
        std::unordered_map<uint32_t, OutputResult2> results;
        for (auto& bkt : buckets_) {
            for (auto& cell : bkt.cells) {
                if (cell.cnt == 0) continue;
                OutputResult2 res = {0, 0, {}};
                
                int b_idx = 0;
                for (auto& pair : cell.sketch.grid) {
                    auto& b = pair.second;
                    res.total_pkts += b.total_cnt;
                    res.total_matching += b.matching_cnt;
                    res.top_k_buckets.push_back({b_idx++, b.total_cnt, b.matching_cnt, b.get_center(), b.sum_sq_diff});
                }
                
                std::sort(res.top_k_buckets.begin(), res.top_k_buckets.end(),
                    [](const auto& a, const auto& b){ return a.b_total > b.b_total; });
                if ((int)res.top_k_buckets.size() > k) res.top_k_buckets.resize(k);
                
                results[cell.id] = res;
            }
        }
        return results;
    }

    double get_memory_mb() const {
        size_t total_bytes = sizeof(*this);

        total_bytes += m_ * sizeof(Bucket3_DS);

        for (const auto& bkt : buckets_) {
            total_bytes += n_ * sizeof(Cell_DS);

            for (const auto& cell : bkt.cells) {
                const auto& grid = cell.sketch.grid;

                total_bytes += grid.bucket_count() * sizeof(void*);

                size_t node_size = sizeof(std::pair<const int, DSBin>);

                total_bytes += grid.size() * node_size;
            }
        }

        return static_cast<double>(total_bytes) / (1024.0 * 1024.0);
    }
private:
    uint16_t sg_, m_;
    uint8_t n_;
    float theta_;
    uint32_t cnt_g_;
    std::vector<Bucket3_DS> buckets_;
};

#endif