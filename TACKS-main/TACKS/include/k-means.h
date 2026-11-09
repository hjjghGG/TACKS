#ifndef PROJ_K_MEANS_H
#define PROJ_K_MEANS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <map>
#include <cstdint>

struct ClusterResult {
    std::vector<int> assignments;
    std::vector<double> centroids;
    std::vector<double> metric;
};

class KMeans {
public:
    static ClusterResult run(const std::vector<uint32_t>& data, int k = C_NUM/2, bool recursive = false) {
        if (data.empty()) return {};

        ClusterResult res = executeKMeans(data, k);

        if (recursive && !data.empty() && k > 1) {
            std::map<int, int> counts;
            for (int id : res.assignments) counts[id]++;

            int max_cluster_id = -1;
            int max_count = -1;
            for (auto const& [id, count] : counts) {
                if (count > max_count) {
                    max_count = count;
                    max_cluster_id = id;
                }
            }

            std::vector<uint32_t> sub_data;
            std::vector<int> sub_data_indices;
            for (size_t i = 0; i < data.size(); ++i) {
                if (res.assignments[i] == max_cluster_id) {
                    sub_data.push_back(data[i]);
                    sub_data_indices.push_back(i);
                }
            }

            if (sub_data.size() > static_cast<size_t>(k)) {
                ClusterResult sub_res = executeKMeans(sub_data, k);


                // reallocate
                int id_offset = k;
                for (size_t i = 0; i < sub_data_indices.size(); ++i) {
                    int original_idx = sub_data_indices[i];
                    res.assignments[original_idx] = id_offset + sub_res.assignments[i];
                }

                // update centroid
                std::vector<double> global_centroids;
                std::map<int, int> old_to_new_id_map;

                for (int i = 0; i < k; ++i) {
                    if (i != max_cluster_id) {
                        old_to_new_id_map[i] = global_centroids.size();
                        global_centroids.push_back(res.centroids[i]);
                    }
                }
                int sub_start_new_id = global_centroids.size();
                for (double c : sub_res.centroids) {
                    global_centroids.push_back(c);
                }

                // map assignments
                for (int& id : res.assignments) {
                    if (id >= id_offset) {
                        id = sub_start_new_id + (id - id_offset);
                    } else {
                        id = old_to_new_id_map[id];
                    }
                }

                res.centroids = global_centroids;
            }
        }

        res.metric = calculateStabilityMetric(data, res.assignments, res.centroids);

        return res;
    }

private:
    static ClusterResult executeKMeans(const std::vector<uint32_t>& data, int k) {
        std::vector<double> points;
        for (uint32_t v : data) points.push_back(static_cast<double>(v));

        int n = points.size();
        int actual_k = std::min(k, n);

        std::vector<double> centroids(actual_k);
        auto minmax = std::minmax_element(points.begin(), points.end());
        double min_val = *minmax.first;
        double range = *minmax.second - min_val;

        for (int i = 0; i < actual_k; ++i) {
            centroids[i] = min_val + (range * (i + 1) / (actual_k + 1));
        }

        std::vector<int> assignments(n, 0);
        bool changed = true;
        int iterations = 0;

        while (changed && iterations++ < 100) {
            changed = false;
            for (int i = 0; i < n; ++i) {
                int best_k = 0;
                double min_dist = std::abs(points[i] - centroids[0]);
                for (int j = 1; j < actual_k; ++j) {
                    double dist = std::abs(points[i] - centroids[j]);
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_k = j;
                    }
                }
                if (assignments[i] != best_k) {
                    assignments[i] = best_k;
                    changed = true;
                }
            }

            std::vector<double> sums(actual_k, 0);
            std::vector<int> counts(actual_k, 0);
            for (int i = 0; i < n; ++i) {
                sums[assignments[i]] += points[i];
                counts[assignments[i]]++;
            }
            for (int j = 0; j < actual_k; ++j) {
                if (counts[j] > 0) centroids[j] = sums[j] / counts[j];
            }
        }

        ClusterResult res;
        res.assignments = assignments;
        res.centroids = centroids;
        return res;
    }

    static std::vector<double> calculateStabilityMetric(
        const std::vector<uint32_t>& data,
        const std::vector<int>& assignments,
        const std::vector<double>& centroids) {

        int k = centroids.size();
        std::vector<double> metric(k+1, 0.0);
        std::vector<int> total_counts(k, 0);
        std::vector<int> stable_counts(k, 0);

        for (size_t i = 0; i < data.size(); ++i) {
            double point = static_cast<double>(data[i]);

            // find closed centroid
            int nearest_cluster = 0;
            double min_dist = std::abs(point - centroids[0]);
            for (int j = 1; j < k; ++j) {
                double dist = std::abs(point - centroids[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    nearest_cluster = j;
                }
            }

            total_counts[nearest_cluster]++;

            // if steady
            double centroid_val = centroids[nearest_cluster];
            double abs_diff = std::abs(point - centroid_val);
            double ratio = (centroid_val != 0.0) ? (abs_diff / std::abs(centroid_val)) :
                            ((point == 0.0) ? 0.0 : std::numeric_limits<double>::infinity());

            if (abs_diff < beta || ratio < alpha) {
                stable_counts[nearest_cluster]++;
            }
        }


        double sum_counts = 0;
        double sum_total_counts = 0;
        // cal
        for (int j = 0; j < k; ++j) {
            if (total_counts[j] > 0) {
                sum_total_counts += total_counts[j];
                sum_counts += stable_counts[j];
                metric[j] = static_cast<double>(stable_counts[j]) / total_counts[j];
            } else {
                metric[j] = 0.0;
            }
        }
        metric[k] = sum_counts / sum_total_counts;
        return metric;
    }
};

#endif