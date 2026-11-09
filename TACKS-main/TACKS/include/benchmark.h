#ifndef PROJ_BENCHMARK_H
#define PROJ_BENCHMARK_H
#include <chrono>
#include <vector>
#include <climits>
#include <fstream>
#include <sstream>
#include <divider.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <filesystem>
#include "k-means.h"
#include "BF.h"
#include "Tracedata.h"


struct ComparisonResult {
    uint32_t flow_id;
    int cluster_count;
    std::vector<double> kmeans_metric;
    std::vector<double> extractor_metrics; //steady flow
    bool extractor_available;
    std::string extractor_centers;
    std::string kmeans_centers;
    double cv_gt_mean = 0.0;
    double cv_est_mean = 0.0;
};

class BenchMark {
public:
    explicit BenchMark(std::string &file) {
        tracedata = readTraceFile(file);
    }

    explicit BenchMark(TraceData trace_data) {
        tracedata = trace_data;
    }

    void savePacketsToFile() {
        std::ofstream outFile("case1.txt");
        if (!outFile.is_open()) {
            // Handle error
            return;
        }

        for (const auto &packet: tracedata.allPackets) {
            outFile << packet.first << " " << packet.second << "\n";
        }
        std::cout << "=== save end ===" << std::endl;
        outFile.close();
    }

    void Error_test(uint16_t sampling_rate, float theta, uint16_t m, uint16_t n) {
        std::cout << "=== Initializing Divider ===" << std::endl;
        std::cout << "Parameters: sampling_rate=" << sampling_rate
                << ", m=" << m << ", n=" << (int) n << ", theta=" << theta << std::endl;
        int bf_n = m * 4;
        DecayingBF bf = DecayingBF(bf_n, 3);
        divider div = divider(sampling_rate, m, n, theta);

        TP start = now();
        long lines = 0;
        for (std::pair<uint32_t, int> p: tracedata.allPackets) {
            ++lines;
            if (lines % 100000 == 1 && lines > 100000) {
                bf.decay();
            }
            if (!bf.insert(p.first)) continue;
            div.insert(p.first, p.second);
        }
        TP end = now();
        auto elapsed = durationms(end, start);

        std::cout << "process " << lines << " packets" << std::endl;
        std::cout << "Elapsed time: " << elapsed / 1000000.0 << " s" << std::endl;
        div.output();

        std::unordered_map<uint32_t, OutputResult> output = div.output2();

        std::vector<ComparisonResult> all_results;


        int total_flows = 0;
        int extractor_available_count = 0;

        int lines_ = 0;
        std::cout << "\n=== Starting Flow Analysis ===" << std::endl;
        std::cout << tracedata.flowLatencies.size() << std::endl;
        for (const auto &l: tracedata.flowLatencies) {
            total_flows++;

            uint32_t id = l.first;

            if (output.find(id) == output.end()) {
                continue;
            }

            ComparisonResult result;
            result.flow_id = l.first;

            // KMeans
            std::vector<uint32_t> latencies = l.second;
            std::vector<uint32_t> intervals;

            if (!latencies.empty()) {
                for (size_t i = 1; i < latencies.size(); ++i) {
                    lines_++;

                    if (latencies[i] - latencies[i - 1] > 65535)continue;
                    intervals.push_back(latencies[i] - latencies[i - 1]);
                }
            }
            if (intervals.empty())continue;

            ClusterResult kmeans_res = KMeans::run(intervals, 2, true);
            result.cluster_count = kmeans_res.centroids.size();

            result.kmeans_metric = kmeans_res.metric;

            std::stringstream kmeans_centers_ss;
            for (size_t i = 0; i < kmeans_res.centroids.size(); ++i) {
                if (i > 0) kmeans_centers_ss << "|";
                kmeans_centers_ss << static_cast<int>(round(kmeans_res.centroids[i]));
            }
            result.kmeans_centers = kmeans_centers_ss.str();

            // Extractor
            result.extractor_available = false;
            result.extractor_centers = "";


            OutputResult ors = output[id];
            std::vector<double> extractor_centers;
            std::vector<double> extractor_counts;
            for (int i = 0; i < C_NUM; ++i) {
                if (ors.centers[i].first > 0) {
                    extractor_centers.push_back(ors.centers[i].first);
                }
                if (ors.centers[i].second > 0) {
                    extractor_counts.push_back(ors.centers[i].second);
                }
            }

            if (!extractor_centers.empty()) {
                result.extractor_available = true;
                extractor_available_count++;

                // Build sum_sq_diffs vector
                std::vector<float> sum_sq_diffs;
                for (int i = 0; i < C_NUM; ++i) {
                    if (ors.centers[i].first > 0) {
                        sum_sq_diffs.push_back(ors.sum_sq_diffs[i]);
                    }
                }

                // Extractor Metric
                result.extractor_metrics = calculate_center_metric(ors.cnts, extractor_counts, intervals,
                                                                   extractor_centers, sum_sq_diffs,
                                                                   result.cv_gt_mean, result.cv_est_mean);

                std::stringstream extractor_centers_ss;
                for (size_t i = 0; i < extractor_centers.size(); ++i) {
                    if (i > 0) extractor_centers_ss << "|";
                    extractor_centers_ss << static_cast<int>(round(extractor_centers[i]));
                }
                result.extractor_centers = extractor_centers_ss.str();
            }

            all_results.push_back(result);

            if (total_flows % 1000 == 0) {
                std::cout << "Processed " << total_flows << " flows..." << std::endl;
            }
        }

        std::cout << "\n=== Analysis Complete ===" << std::endl;
        std::cout << "Total flows analyzed: " << total_flows << std::endl;
        std::cout << "Flows with Extractor available: " << extractor_available_count << std::endl;
        std::cout << "Coverage rate: " << (100.0 * extractor_available_count / total_flows) << "%" << std::endl;

        double throughput = (elapsed > 0) ? (1000000.0 * lines / elapsed) : 0;
        outputToCSV_Simple(all_results, static_cast<uint32_t>(throughput), sampling_rate, m, n, theta, bf_n);
       // outputToCSV(all_results, sampling_rate, m, n, theta);
    }


    static double calculate_memory_usage_mb(uint32_t bf_n, uint16_t m, uint8_t n) {
        // 1. DecayingBF size (bytes)
        size_t bf_size = (bf_n + 1) / 2;

        // 2. PSquare size (bytes)
        int marker_count = 2 * Q_NUM + 3;
        size_t psquare_size = sizeof(PSquare);

        // 3. List3 size (bytes)
        size_t list3_size = sizeof(List3);

        // 4. Cell size
        size_t cell_size = (sizeof(uint32_t) * 3) + // id, cnt, t_x
                           psquare_size +
                           list3_size;

        // 5. Bucket size
        size_t bucket_size = (cell_size * n) + sizeof(uint16_t); // cells + cnt_n

        // Total size in bytes
        size_t total_bytes = bf_size + (bucket_size * m);

        return static_cast<double>(total_bytes) / (1024.0 * 1024.0);
    }

private:
    TraceData tracedata;

    typedef std::chrono::high_resolution_clock::time_point TP;

    inline TP now() {
        return std::chrono::high_resolution_clock::now();
    }

    inline double durationms(TP finish, TP start) {
        return std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1, 1000000> > >(finish - start).
                count();
    }

    std::vector<double> calculate_center_metric(
        const std::vector<uint16_t> &cnt,
        const std::vector<double> &counts,
        const std::vector<uint32_t> &data,
        const std::vector<double> &centers,
        const std::vector<float> &sum_sq_diffs,
        double &cv_gt_out,
        double &cv_est_out) {
        size_t k = centers.size();
        cv_gt_out = 0.0;
        cv_est_out = 0.0;
        if (k == 0 || cnt.size() < k || counts.size() < k)
            return {};

        std::vector<double> result;

        double sum_counts = 0.0;
        double sum_cnts = 0.0;
        for (size_t i = 0; i < k; ++i) {
            double ratio = 0.0;
            sum_counts += counts[i];
            sum_cnts += cnt[i];
            if (cnt[i] > 0) {
                ratio = static_cast<double>(counts[i]) / cnt[i];
            }
            result.push_back(ratio);
        }
        result.push_back(sum_cnts == 0 ? 0 : sum_counts / sum_cnts);

        std::vector<int> total_assigned(k, 0);
        std::vector<int> stable_assigned(k, 0);
        // Welford accumulators for ground-truth CV per bucket
        std::vector<double> gt_M2(k, 0.0);

        for (uint32_t point: data) {
            double dp = static_cast<double>(point);
            int best = 0;
            double min_dist = std::abs(dp - centers[0]);
            for (size_t j = 1; j < k; ++j) {
                double dist = std::abs(dp - centers[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best = static_cast<int>(j);
                }
            }

            total_assigned[best]++;

            double center_val = centers[best];
            bool is_stable = false;
            if (min_dist < beta) {
                is_stable = true;
            } else if (center_val > 0.0) {
                double ratio = min_dist / center_val;
                if (ratio < alpha) is_stable = true;
            }

            if (is_stable) stable_assigned[best]++;

            // Welford update for ground-truth CV
            double delta = dp - centers[best];
            gt_M2[best] += delta * delta;
        }

        for (size_t i = 0; i < k; ++i) {
            double stab = 0.0;
            if (total_assigned[i] > 0) {
                stab = static_cast<double>(stable_assigned[i]) / total_assigned[i];
            }
            result.push_back(stab);
        }

        int sum_total = 0, sum_stable = 0;
        for (size_t i = 0; i < k; ++i) {
            sum_total += total_assigned[i];
            sum_stable += stable_assigned[i];
        }
        // global stability
        double global_stab = (sum_total > 0) ? static_cast<double>(sum_stable) / sum_total : 0.0;
        result.push_back(global_stab);

        // Compute CV per bucket (ground truth from Welford, estimated from sketch sum_sq)
        double cv_gt_sum = 0.0, cv_est_sum = 0.0;
        int cv_buckets = 0;
        for (size_t i = 0; i < k; ++i) {
            if (total_assigned[i] > 1 && centers[i] > 0) {
                double std_gt = std::sqrt(gt_M2[i] / total_assigned[i]);
                double cv_gt = std_gt / centers[i];
                cv_gt_sum += cv_gt;

                if (i < sum_sq_diffs.size() && i < cnt.size() && cnt[i] > 0 && centers[i] > 0) {
                    double est_var = static_cast<double>(sum_sq_diffs[i]) / cnt[i];
                    double cv_est = std::sqrt(std::max(0.0, est_var)) / centers[i];
                    cv_est_sum += cv_est;
                }
                cv_buckets++;
            }
        }
        if (cv_buckets > 0) {
            cv_gt_out = cv_gt_sum / cv_buckets;
            cv_est_out = cv_est_sum / cv_buckets;
        }

        return result;
    }

    void outputToCSV_Simple(const std::vector<ComparisonResult> &results, uint32_t thp, uint16_t sampling_rate,
                            uint16_t m, uint16_t n, float theta, uint16_t bf_n) {
        std::string filename = "../result/qsquare/"+tracedata.filename+"_"+"Error_m" + std::to_string(m * n) + "theta_" +
                               std::to_string(theta) + "s_" + std::to_string(sampling_rate) + ".csv";

        std::filesystem::create_directories("../result");

        std::ofstream csv_file(filename, std::ios_base::app);
        if (!csv_file.is_open()) {
            std::cout << "Failed to save CSV file: " << filename << std::endl;
            return;
        }

        double sum_mem = 0.0;
        double sum_metric_k = 0.0;
        double sum_metric_e = 0.0;
        double sum_metric_i = 0.0;
        double sum_cv_gt = 0.0;
        double sum_cv_est = 0.0;
        int num = 0;

        for (const auto &r: results) {
            if (!r.extractor_available) continue;

            if (r.kmeans_metric.empty()) {
                continue;
            }
            if (r.extractor_metrics.size() < 2) {
                continue;
            }

            num++;
            float mem = calculate_memory_usage_mb(bf_n, m, n);
            sum_mem += mem;
            sum_metric_k += r.kmeans_metric.back();

            size_t half_size = r.extractor_metrics.size() / 2;
            if (half_size > 0) {
                sum_metric_e += r.extractor_metrics[half_size - 1];
            } else {
                sum_metric_e += 0.0;
            }

            sum_metric_i += r.extractor_metrics.back();
            sum_cv_gt += r.cv_gt_mean;
            sum_cv_est += r.cv_est_mean;
        }

        if (num == 0) {
            std::cout << "Warning: No valid flows to output" << std::endl;
            csv_file.close();
            return;
        }

        bool is_new_file = (!std::filesystem::exists(filename) ||
                            std::filesystem::file_size(filename) == 0);
        if (is_new_file) {
            csv_file << "Memory,THP,KMeansMetric,ExtractorMetric,IdealMetric,Error,CV_Mean,CV_ARE\n";
        }

        double avg_cv_gt = sum_cv_gt / num;
        double avg_cv_est = sum_cv_est / num;
        double cv_are = (avg_cv_gt > 0) ? std::abs(avg_cv_est - avg_cv_gt) / avg_cv_gt : 0.0;
        csv_file << sum_mem / num << ","
                << thp << "," << sum_metric_k / num
                << "," << sum_metric_e / num
                << "," << sum_metric_i / num
                << "," << std::abs(sum_metric_e-sum_metric_i) / sum_metric_i
                << "," << avg_cv_gt
                << "," << cv_are;
        csv_file << "\n";
        csv_file.close();
        std::cout << "\nResults saved to: " << filename << std::endl;
    }


    void outputToCSV(const std::vector<ComparisonResult> &results, uint16_t sampling_rate, uint16_t m, uint16_t n,
                     float theta) {
        std::string filename = "../result/qsquare/Error_m" + std::to_string(m) + "n_" + std::to_string(n) + "theta_" +
                               std::to_string(theta)
                               + "s_" + std::to_string(sampling_rate) + ".csv";
        std::ofstream csv_file(filename);
        if (!csv_file.is_open()) {
            std::cout << "Failed to create CSV file: " << filename << std::endl;
            return;
        }

        csv_file << "FlowID,ClusterCount,KMeansCenters,KMeansMetrics,"
                << "ExtractorAvailable,ExtractorCenters,ExtractorMetrics\n";

        for (const auto &r: results) {
            csv_file << r.flow_id << ","
                    << r.cluster_count << ",";

            // KMeans centers
            csv_file << "\"" << r.kmeans_centers << "\",";

            // KMeans metrics
            csv_file << "\"";
            for (size_t i = 0; i < r.kmeans_metric.size(); ++i) {
                if (i > 0) csv_file << "|";
                csv_file << r.kmeans_metric[i];
            }
            csv_file << "\",";

            // Extractor
            csv_file << (r.extractor_available ? 1 : 0) << ",";

            // Extractor centers
            csv_file << "\"" << r.extractor_centers << "\",";

            // Extractor metrics（ratio1|ratio2|…|stab1|stab2|…|global_stab）
            csv_file << "\"";
            for (size_t i = 0; i < r.extractor_metrics.size(); ++i) {
                if (i > 0) csv_file << "|";
                csv_file << r.extractor_metrics[i];
            }
            csv_file << "\"\n";
        }

        csv_file.close();
        std::cout << "\nResults saved to: " << filename << std::endl;
    }
};

#endif //PROJ_BENCHMARK_H
