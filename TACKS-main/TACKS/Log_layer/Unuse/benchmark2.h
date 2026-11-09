#pragma once
#ifndef BENCHMARK2_H
#define BENCHMARK2_H

#include <chrono>
#include <vector>
#include <climits>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <algorithm>

#include "../Log_layer/log_sketch.h"
#include "BF.h"
#include "Parameter.h"
#include "Tracedata.h"


struct ComparisonResult2 {
    uint32_t flow_id;
    std::vector<double> extractor_metrics;
    bool extractor_available;
    std::string extractor_centers;
};

class BenchMark2 {
public:
    explicit BenchMark2(std::string &file) {
        traceData = readTraceFile(file);
    }

    explicit BenchMark2(TraceData trace_data) {
        traceData = trace_data;
    }

    void Error_test(uint16_t sampling_rate, float theta, uint16_t m, uint16_t n) {
        std::cout << "=== Initializing L1 (LogSketch) Divider ===" << std::endl;
        std::cout << "Parameters: sampling_rate=" << sampling_rate
                  << ", m=" << m << ", n=" << (int)n << ", theta=" << theta << std::endl;

        int bf_n = m * 3;
        DecayingBF bf = DecayingBF(bf_n, 3);
        L1 l1 = L1(sampling_rate, m, n, theta);

        TP start = now();
        long lines = 0;
        for (const auto &p : traceData.allPackets) {
            ++lines;
            if (lines % 5000 == 1 && lines > 5000) {
                bf.decay();
            }
            if (!bf.insert(p.first)) continue;
            l1.insert(p.first, p.second);
        }
        TP end = now();
        auto elapsed = durationms(end, start);

        std::cout << "Processed " << lines << " packets" << std::endl;
        std::cout << "Elapsed time: " << elapsed / 1000000.0 << " s" << std::endl;

        std::unordered_map<uint32_t, OutputResult2> output = l1.output_top_k(1000);
        std::vector<ComparisonResult2> all_results;

        int total_flows = 0;
        int extractor_available_count = 0;

        for (const auto &it : traceData.flowLatencies) {
            total_flows++;
            uint32_t id = it.first;

            if (output.find(id) == output.end()) continue;

            const std::vector<uint32_t> &latencies = it.second;
            std::vector<uint32_t> intervals;
            for (size_t i = 1; i < latencies.size(); ++i) {
                uint32_t diff = latencies[i] - latencies[i - 1];
                if (diff <= 65535) intervals.push_back(diff);
            }
            if (intervals.empty()) continue;

            ComparisonResult2 result;
            result.flow_id = id;
            result.extractor_available = true;
            extractor_available_count++;

            const OutputResult2 &ors = output[id];
            std::vector<double> centers;
            std::vector<uint16_t> b_cnts;
            std::vector<double> b_matchings;

            for (const auto &bkt : ors.top_k_buckets) {
                centers.push_back(bkt.b_center);
                b_cnts.push_back(static_cast<uint16_t>(bkt.b_total));
                b_matchings.push_back(static_cast<double>(bkt.b_matching));
            }

            result.extractor_metrics = calculate_center_metric(b_cnts, b_matchings, intervals, centers);

            std::stringstream ss;
            for (size_t i = 0; i < centers.size(); ++i) {
                ss << (i > 0 ? "|" : "") << static_cast<int>(std::round(centers[i]));
            }
            result.extractor_centers = ss.str();
            all_results.push_back(result);
        }

        std::cout << "Flows analyzed: " << total_flows << ", Available: " << extractor_available_count << std::endl;
        
        double throughput = (elapsed > 0) ? (1000000.0 * lines / elapsed) : 0;
        outputToCSV_Simple(all_results, static_cast<uint32_t>(throughput), sampling_rate, m, n, theta, bf_n);
    }

private:
    TraceData traceData;
    typedef std::chrono::high_resolution_clock::time_point TP;
    inline TP now() { return std::chrono::high_resolution_clock::now(); }
    inline double durationms(TP finish, TP start) {
        return std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1, 1000000>>>(finish - start).count();
    }

    std::vector<double> calculate_center_metric(
        const std::vector<uint16_t> &cnt,
        const std::vector<double> &counts,
        const std::vector<uint32_t> &data,
        const std::vector<double> &centers) {
        
        size_t k = centers.size();
        if (k == 0) return {};

        std::vector<double> result;
        double sum_counts = 0.0, sum_cnts = 0.0;

        for (size_t i = 0; i < k; ++i) {
            sum_counts += counts[i];
            sum_cnts += cnt[i];
            result.push_back(cnt[i] > 0 ? counts[i] / cnt[i] : 0.0);
        }
        result.push_back(sum_cnts == 0 ? 0 : sum_counts / sum_cnts); // ExtractorMetric[cite: 3]

        std::vector<int> total_assigned(k, 0), stable_assigned(k, 0);
        for (uint32_t point : data) {
            double dp = static_cast<double>(point);
            int best = 0;
            double min_dist = std::abs(dp - centers[0]);
            for (size_t j = 1; j < k; ++j) {
                double dist = std::abs(dp - centers[j]);
                if (dist < min_dist) { min_dist = dist; best = j; }
            }
            total_assigned[best]++;
            if (min_dist < beta || (centers[best] > 0 && (min_dist / centers[best]) < alpha))
                stable_assigned[best]++;
        }

        int sum_total = 0, sum_stable = 0;
        for (size_t i = 0; i < k; ++i) {
            sum_total += total_assigned[i];
            sum_stable += stable_assigned[i];
        }
        result.push_back(sum_total > 0 ? static_cast<double>(sum_stable) / sum_total : 0.0); // IdealMetric[cite: 3]

        return result;
    }

    double calculate_logsketch_memory_usage_mb(uint32_t bf_n, uint16_t m, uint8_t n) {
        size_t bf_size = (bf_n + 1) / 2;
        size_t log_bucket_struct_size = sizeof(LogBucket);
        int num_log_buckets = MAX_BIN;
        size_t sketch_total_size = sizeof(double) + sizeof(float) + sizeof(std::vector<LogBucket>) + (num_log_buckets * log_bucket_struct_size);
        size_t cell_size = (sizeof(uint32_t) * 3) + sketch_total_size;
        size_t bucket_total_size = (cell_size * n) + sizeof(std::vector<Cell_sp>) + sizeof(uint16_t);
        size_t divider_overhead = sizeof(uint16_t) * 2 + sizeof(uint8_t) + sizeof(float) + sizeof(uint32_t) + sizeof(std::vector<Bucket2>);
        
        return static_cast<double>(bf_size + (bucket_total_size * m) + divider_overhead) / (1024.0 * 1024.0);
    }

    void outputToCSV_Simple(const std::vector<ComparisonResult2> &results, uint32_t thp, uint16_t sampling_rate,
                            uint16_t m, uint16_t n, float theta, uint16_t bf_n) {
        std::string filename = "../result/log/" + traceData.filename+"Error_2_m"  + std::to_string(m * n) + "theta_" +
                               std::to_string(theta) + "s_" + std::to_string(sampling_rate) + ".csv";
        
        std::filesystem::create_directories("../result");
        bool is_new_file = (!std::filesystem::exists(filename) || std::filesystem::file_size(filename) == 0);
        std::ofstream csv_file(filename, std::ios_base::app);

        if (is_new_file) csv_file << "Memory,THP,ExtractorMetric,IdealMetric,Error\n";

        double sum_metric_e = 0.0, sum_metric_i = 0.0;
        int num = 0;
        for (const auto &r : results) {
            if (!r.extractor_available || r.extractor_metrics.empty()) continue;
            num++;
            sum_metric_e += r.extractor_metrics[r.extractor_metrics.size() - 2];
            sum_metric_i += r.extractor_metrics.back();
        }

        if (num > 0) {
            csv_file << calculate_logsketch_memory_usage_mb(bf_n, m, n) << ","
            << thp << "," << sum_metric_e / num << "," << sum_metric_i / num
                << "," << std::abs(sum_metric_e-sum_metric_i) / sum_metric_i << "\n";
        }
        csv_file.close();
    }
};

#endif //BENCHMARK2_H