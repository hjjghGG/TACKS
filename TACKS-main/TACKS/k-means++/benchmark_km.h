#pragma once
#ifndef BENCHMARK_KM_H
#define BENCHMARK_KM_H

#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "km_sketch.h" // 替换为刚刚生成的 km_sketch 头文件
#include "BF.h"
#include "Tracedata.h"

struct ComparisonResultKM {
    uint32_t flow_id;
    std::vector<double> metrics;
};

struct MetricsResultKM {
    double extractor_metric;
    double ideal_metric;
    double cv_gt_mean;
    double cv_est_mean;
};

class BenchMarkKM {
public:
    explicit BenchMarkKM(TraceData& trace_data) : traceData(trace_data) {}

    void Error_test(uint16_t sampling_rate, float theta, uint16_t m, uint16_t n) {
        std::cout << "=== Running KM_Sketch Divider Test ===" << std::endl;

        int bf_n = m * 3;
        DecayingBF bf(bf_n, 3);
        L8 L8(sampling_rate, m, n, theta);

        auto start = std::chrono::high_resolution_clock::now();
        long lines = 0;
        for (const auto &p : traceData.allPackets) {
            if (++lines % 100000 == 1 && lines > 100000) bf.decay();
            if (!bf.insert(p.first)) continue;
            L8.insert(p.first, p.second);
        }
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::unordered_map<uint32_t, OutputResult5> output = L8.output_top_k((TOP));

        double sum_e = 0, sum_i = 0, sum_cv_gt = 0, sum_cv_est = 0;
        int count = 0;

        for (const auto &it : traceData.flowLatencies) {
            if (output.find(it.first) == output.end()) continue;

            std::vector<uint32_t> intervals;
            for (size_t i = 1; i < it.second.size(); ++i) {
                uint32_t diff = it.second[i] - it.second[i-1];
                if (diff <= 65535) intervals.push_back(diff);
            }
            if (intervals.empty()) continue;

            const auto& res = output[it.first];
            std::vector<double> centers;
            std::vector<uint16_t> b_cnts;
            std::vector<double> b_matchings;
            std::vector<float> b_sum_sq;
            for(auto& b : res.top_k_buckets) {
                centers.push_back(b.b_center);
                b_cnts.push_back(b.b_total);
                b_matchings.push_back(b.b_matching);
                b_sum_sq.push_back(b.b_sum_sq);
            }

            auto mres = calculate_metrics(b_cnts, b_matchings, b_sum_sq, intervals, centers);
            sum_e += mres.extractor_metric;
            sum_i += mres.ideal_metric;
            sum_cv_gt += mres.cv_gt_mean;
            sum_cv_est += mres.cv_est_mean;
            count++;
        }

        double memory = calculate_memory_mb(bf_n, m, n);
        double thp = (elapsed > 0) ? (lines / (elapsed / 1000000.0)) : 0;

        save_to_csv(memory, thp, sum_e/count, sum_i/count,
                    sum_cv_gt/count, sum_cv_est/count, sampling_rate, m, n, theta);
        
        std::cout << "Done. elapsed: " << elapsed / 1000000.0 << "s, THP: " << thp/1e6 << " Mpps" << std::endl;
    }

private:
    TraceData& traceData;

    MetricsResultKM calculate_metrics(const std::vector<uint16_t>& cnt, const std::vector<double>& match,
                                      const std::vector<float>& sum_sq,
                                      const std::vector<uint32_t>& data, const std::vector<double>& centers) {
        MetricsResultKM mr = {0.0, 0.0, 0.0, 0.0};
        size_t k = centers.size();
        if (k == 0) return mr;

        double s_match = 0, s_cnt = 0;
        for(size_t i=0; i<cnt.size() && i<k; ++i) { s_match += match[i]; s_cnt += cnt[i]; }

        std::vector<double> gt_M2(k, 0.0);
        std::vector<int> gt_count(k, 0);

        int total_assigned = 0, stable_assigned = 0;
        for (uint32_t p : data) {
            double dp = (double)p;
            double min_dist = 1e9;
            int best = 0;
            for(size_t j=0; j<k; ++j) {
                double d = std::abs(dp - centers[j]);
                if(d < min_dist) { min_dist = d; best = (int)j; }
            }
            total_assigned++;
            if (min_dist < beta || (centers[best] > 0 && (min_dist/centers[best]) < alpha)) stable_assigned++;

            gt_count[best]++;
            double delta = dp - centers[best];
            gt_M2[best] += delta * delta;
        }

        mr.extractor_metric = (s_cnt > 0 ? s_match/s_cnt : 0.0);
        mr.ideal_metric = (total_assigned > 0 ? (double)stable_assigned/total_assigned : 0.0);

        double cv_gt_sum = 0.0, cv_est_sum = 0.0;
        int cv_buckets = 0;
        for (size_t i = 0; i < k; ++i) {
            if (gt_count[i] > 1 && centers[i] > 0) {
                double std_gt = std::sqrt(gt_M2[i] / gt_count[i]);
                double cv_gt = std_gt / centers[i];
                cv_gt_sum += cv_gt;

                if (i < sum_sq.size() && i < cnt.size() && cnt[i] > 0 && centers[i] > 0) {
                    double est_var = (double)sum_sq[i] / cnt[i];
                    double cv_est = std::sqrt(std::max(0.0, est_var)) / centers[i];
                    cv_est_sum += cv_est;
                }
                cv_buckets++;
            }
        }
        if (cv_buckets > 0) {
            mr.cv_gt_mean = cv_gt_sum / cv_buckets;
            mr.cv_est_mean = cv_est_sum / cv_buckets;
        }

        return mr;
    }

    double calculate_memory_mb(uint32_t bf_n, uint16_t m, uint8_t n) {
        size_t bf_mem = (bf_n + 1) / 2;
        size_t sketch_mem = sizeof(std::vector<KMBucket>) + (S_NUM * sizeof(KMBucket));
        size_t cell_mem = sizeof(uint32_t)*3 + sketch_mem;
        size_t total = bf_mem + m * (sizeof(std::vector<Cell_km>) + n * cell_mem + sizeof(uint16_t));
        return (double)total / (1024.0 * 1024.0);
    }

    void save_to_csv(double mem, double thp, double e, double i, double cv_gt, double cv_est, int s, int m, int n, float theta) {
        std::string dir_path = "../result/km/";
        std::string filename = dir_path + traceData.filename + "KM_Error_m" + std::to_string(m*n) + "theta_"+std::to_string(theta)+".csv";
        std::filesystem::create_directories(dir_path);

        bool is_new_file = (!std::filesystem::exists(filename) || std::filesystem::file_size(filename) == 0);
        std::ofstream csv_file(filename, std::ios_base::app);

        double cv_are = (cv_gt > 0) ? std::abs(cv_est - cv_gt) / cv_gt : 0.0;
        if (is_new_file) csv_file << "Memory,THP,ExtractorMetric,IdealMetric,Error,CV_Mean,CV_ARE\n";
        csv_file << mem << "," << (int)thp << "," << e << "," << i << ","
            << std::abs(e-i)/i << "," << cv_gt << "," << cv_are << "\n";
    }
};

#endif