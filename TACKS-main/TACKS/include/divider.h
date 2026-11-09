#ifndef PROJ_DIVIDER_H
#define PROJ_DIVIDER_H

#include <cstdint>
#include <vector>
#include <hash.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "p_square.h"
#include "Extractor.h"
//#define Test_Mode


struct OutputResult {
    std::vector<std::pair<uint16_t,uint16_t>> centers;
    std::vector<uint16_t> cnts;
    std::vector<float> sum_sq_diffs;
    uint32_t cnt;
    uint32_t cnt1;
};

struct Cell {
    uint32_t id;
    uint32_t cnt;
    uint32_t t_x;
    PSquare p_square;
    List3 l3;

    Cell() : id(0), cnt(0), t_x(0) {
        p_square = PSquare();
        l3 = List3();
    }

};

struct Bucket {
    std::vector<Cell> cells;
    uint16_t cnt_n;

    explicit Bucket(int n) : cells(n), cnt_n(0) {
        for (auto& cell : cells) {
            cell = Cell();
        }
    }
};

class divider {
public:
    explicit divider(uint16_t sampling_rate, uint16_t m, uint8_t n, float theta):
    sg_(sampling_rate), m(m), n(n), theta(theta), cnt_g_(0) {

        buckets.resize(m, Bucket(n));
    }


    bool should_sample() {
        cnt_g_++;
        if (cnt_g_ % sg_ == 0) {
            return true;
        }
        return false;
    }

    static void update_centers(Cell& cell) {

        std::vector<float> res = {};
        for(int j=1; j<Q_NUM+1; ++j) {
            res.emplace_back((cell.p_square.get_quantile(j)));
        }
        cell.l3.set_centers(fast_gap_cluster(res));

    }

    void insert(uint32_t id, uint32_t value) {
    if (!should_sample()) {
        return;
    }

    uint32_t idx = hash_func(id) % m;
    Bucket& bucket = buckets[idx];

    uint32_t hash_pos = hash_func(id) % n;
    int target_idx = static_cast<int>(hash_pos);


    if (bucket.cells[target_idx].id == id && bucket.cells[target_idx].cnt > 0) {
        bucket.cells[target_idx].cnt++;
        uint16_t latency= round((double)(value - bucket.cells[target_idx].t_x)/sg_);
        bucket.cells[target_idx].t_x = value;
        bucket.cells[target_idx].p_square.update(latency);
        if ( bucket.cells[target_idx].cnt > std::max(bucket.cnt_n,(uint16_t)15) &&
            bucket.cells[target_idx].p_square.get_max()>1
            && bucket.cells[target_idx].cnt % 15 == 1)
        {
            update_centers(bucket.cells[target_idx]);
        }
        if (bucket.cells[target_idx].l3.centers[0]>0)bucket.cells[target_idx].l3.update(latency);

        return;
    }

    int min_idx = 0;
    uint16_t min_cnt = UINT16_MAX;
    bool found = false;
    int empty_idx = -1;
    for (int i = 0; i < n; ++i) {
        if (bucket.cells[i].id == id && bucket.cells[i].cnt > 0) {
            bucket.cells[i].cnt++;
            uint16_t latency= round((double)(value - bucket.cells[i].t_x)/sg_);
            bucket.cells[i].t_x = value;
            bucket.cells[i].p_square.update(latency);
            if ( bucket.cells[i].cnt > std::min(bucket.cnt_n,(uint16_t)50) &&
            bucket.cells[i].p_square.get_m_count()>15
            && bucket.cells[i].cnt%15 == 1)
            {
                update_centers(bucket.cells[i]);
            }
            if (bucket.cells[i].l3.centers[0]>0)
            {
                bucket.cells[i].l3.update(latency);
            }

            found = true;
            break;
        }

        if (bucket.cells[i].cnt == 0) {
            empty_idx = i;
        }

        // min_count
        if (bucket.cells[i].cnt < min_cnt) {
            min_cnt = bucket.cells[i].cnt;
            min_idx = i;
        }
    }

        if (!found && empty_idx >= 0) {

            bucket.cells[empty_idx].id = id;
            bucket.cells[empty_idx].cnt = 1;
            bucket.cells[empty_idx].t_x = value;
            bucket.cells[empty_idx].p_square = PSquare(); // initial p-square l3
            bucket.cells[empty_idx].l3 = List3();
            return;
        }
    if (found) return;

    // (Collision)
    bucket.cnt_n++;

    // repalcement
    if ((double)bucket.cnt_n * theta > (double)min_cnt) {

        bucket.cells[min_idx].id = id;
        bucket.cells[min_idx].cnt = 1;
        bucket.cells[min_idx].t_x = value;
        bucket.cells[min_idx].p_square = PSquare(); // reset p-square l3
        bucket.cells[min_idx].l3 = List3();
        bucket.cnt_n = 0; // reset cnt_n
    } else {
    }
}


    uint16_t query(uint32_t id) {
        uint32_t idx = hash_func(id) % m;
        Bucket& bucket = buckets[idx];

        uint32_t hash_pos = hash_func(id) % n;
        int target_idx = static_cast<int>(hash_pos);

        if (bucket.cells[target_idx].id == id) {
            return bucket.cells[target_idx].cnt;
        }

        for (int i = 0; i < n; ++i) {
            if (bucket.cells[i].id == id) {
                return bucket.cells[i].cnt;
            }
        }
        return 0;
    }



    static std::vector<uint16_t> fast_gap_cluster(const std::vector<float>& quantiles) {
        int size = (int)quantiles.size();

        // bound check
        if (size <= 0 || C_NUM <= 0) return {};
        if (size <= C_NUM) {
            std::vector<uint16_t> centers;
            for (float v : quantiles) centers.push_back((uint16_t)std::max(1.0f, std::round(v)));
            return centers;
        }

        std::vector<std::pair<double, int>> gaps;
        for (int i = 0; i < size - 1; ++i) {
            double gap = (quantiles[i] == 0) ? quantiles[i+1] : (quantiles[i+1] / quantiles[i]);
            gaps.emplace_back(gap, i);
        }

        // 2. C_NUM - 1 biggest
        int num_splits = C_NUM - 1;
        std::partial_sort(gaps.begin(), gaps.begin() + num_splits, gaps.end(),
            [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                return a.first > b.first;
            });

        std::vector<int> split_indices;
        for (int i = 0; i < num_splits; ++i) {
            split_indices.push_back(gaps[i].second);
        }
        std::sort(split_indices.begin(), split_indices.end());

        std::vector<uint16_t> centers;
        int start_idx = 0;

        // mean with min value 1
        auto add_center = [&](int start, int end) {
            if (start > end) return;
            double sum = 0;
            for (int i = start; i <= end; ++i) sum += quantiles[i];
            uint16_t center = (uint16_t)std::max(1.0, round(sum / (end - start + 1)));
            centers.push_back(center);
        };

        for (int split_idx : split_indices) {
            add_center(start_idx, split_idx);
            start_idx = split_idx + 1;
        }

        add_center(start_idx, size - 1);

        // Ensure no two consecutive centers are equal
        for (size_t i = 1; i < centers.size(); ++i) {
            if (centers[i] <= centers[i-1]) {
                centers[i] = centers[i-1] + 1;
            }
        }

        return centers;
    }

    void output() {
        std::cout << "--- Divider Buckets Summary (Max|Min Count Cell per Bucket) ---" << std::endl;
        for (int i = 0; i < m; ++i) {
            Bucket& bucket = buckets[i];
            int max_idx = -1;
            int min_idx = -1;
            uint16_t max_cnt = 0;
            uint16_t min_cnt = UINT16_MAX;

            for (int j = 0; j < n; ++j) {
                if (bucket.cells[j].cnt > max_cnt) {
                    max_cnt = bucket.cells[j].cnt;
                    max_idx = j;
                }

                if (bucket.cells[j].cnt < min_cnt) {
                    min_cnt = bucket.cells[j].cnt;
                    min_idx = j;
                }
            }

            if (max_idx != -1) {
                Cell& best_cell = bucket.cells[max_idx];
                std::cout << "Bucket [" << i << "] -> Cell[" << max_idx << "] "
                          << "id: " << best_cell.id << ", Count: " << best_cell.cnt << "\n  Quantiles: ";

                for (int k = 1; k < Q_NUM+1; ++k) {
                    double q_val = best_cell.p_square.get_quantile(k);
                    std::cout << q_val << (k == Q_NUM ? "" : ", ");
                }
                std::cout << std::endl;
                best_cell.l3.output();
            } else {
                std::cout << "Bucket [" << i << "] is empty." << std::endl;
            }

            if (min_idx != -1) {
                Cell& best_cell = bucket.cells[min_idx];
                std::cout << "Bucket [" << i << "] -> Cell[" << min_idx << "] "
                          << "id: " << best_cell.id << ", Count: " << best_cell.cnt << "\n  Quantiles: ";

                for (int k = 1; k < 6; ++k) {
                    double q_val = best_cell.p_square.get_quantile(k);
                    std::cout << q_val << (k == 5 ? "" : ", ");
                }
                std::cout << std::endl;
                best_cell.l3.output();
            } else {
                std::cout << "Bucket [" << i << "] is empty." << std::endl;
            }
        }
        std::cout << "--------------------------------------------------------" << std::endl;
    }

    std::unordered_map<uint32_t, OutputResult> output2() {
        std::unordered_map<uint32_t, OutputResult> op;
        for (int i = 0; i < m; ++i) {
            Bucket& bucket = buckets[i];
            for (Cell c : bucket.cells) {
                OutputResult rs = OutputResult();
                rs.cnt = c.cnt;
                std::vector<std::pair<uint16_t,uint16_t>> centers(C_NUM);
                std::vector<uint16_t> counts(C_NUM);
                std::vector<float> sum_sq(C_NUM);
                for (int j = 0; j < C_NUM; ++j) {
                    centers[j].first = c.l3.centers[j];
                    centers[j].second = c.l3.counts[j];
                    counts[j] = (c.l3.cnt[j]);
                    sum_sq[j] = c.l3.sum_sq_diff[j];
                }
                rs.centers = centers;
                rs.cnts = counts;
                rs.sum_sq_diffs = sum_sq;
                op[c.id] = rs;
            }
        }
        return op;
    }


private:
    uint16_t sg_; //sample rate
    uint32_t cnt_g_; //num
    uint16_t m;
    uint8_t n;
    float theta;

    std::vector<Bucket> buckets;

    static uint32_t hash_func(uint32_t id) {
        return Hash::BOBHash32((uint8_t*)&id, sizeof(uint32_t), SEED);
    }

};

#endif //PROJ_DIVIDER_H