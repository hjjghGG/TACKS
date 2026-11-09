#include <iostream>
#include <fstream>
#include <string>
#include <benchmark.h>

#include "Log_layer/benchmark2.h"
#include  "SimpleLayer/benchmark3.h"
#include  "DSSketch/benchmark4.h"
#include "k-means++/benchmark_km.h"


// void estimate_max_q_c_num(double target_mb, uint32_t bf_n, uint16_t m, uint8_t n) {
//     double target_bytes = target_mb * 1024.0 * 1024.0;
//     size_t bf_size = (bf_n + 1) / 2; // proposed
//
//     double remaining_bytes = target_bytes - bf_size;
//     if (remaining_bytes <= 0) {
//         std::cout << "Target memory too small for the given Bloom Filter." << std::endl;
//         return;
//     }
//
//     // (remaining_bytes / m - 2) / n
//     double bytes_per_cell = (remaining_bytes / m - sizeof(uint16_t)) / n;
//
//     int best_q = 0;
//     int best_c = 0;
//
//     for (int q = 1; q < 100; ++q) {
//         int c = (q + 1) / 2;
//
//         int marker_count = 2 * q + 3;
//         size_t ps_size = 4 + (4 * marker_count) + (2 * marker_count) + (4 * marker_count);
//         size_t l3_size = (2 * c * 3) + (2 * c) + (1 * c) + 1;
//         size_t current_cell_size = 12 + ps_size + l3_size;
//
//         if (current_cell_size <= bytes_per_cell) {
//             best_q = q;
//             best_c = c;
//         } else {
//             break;
//         }
//     }
//
//     std::cout << "For " << target_mb << " MB limit:" << std::endl;
//     std::cout << "Max Q_NUM: " << best_q << std::endl;
//     std::cout << "Max C_NUM: " << best_c << std::endl;
// }

int main() {


    std::string filelist[] = {
        "MIRAI_EI.txt",
       "MIRAI.txt",
       "case.txt",
        "MAWI.txt"
    };
    for ( auto file : (filelist))
    {
        TraceData tracedata = readTraceFile(file);

         auto b = BenchMark(tracedata);
         auto c = BenchMark2(tracedata);
         auto d = BenchMark3(tracedata);
         auto e = BenchMark4(tracedata);
        auto g = BenchMarkKM(tracedata);
        int list[] = {
            100,
            200,
            500
        };
        for (auto i : list) {
             b.Error_test(1,0.1,i,5);
        //  c.Error_test(1,0.1,i,5);
        //      d.Error_test(1,0.1,i,5);
            e.Error_test(1,0.1,i,5);
           g.Error_test(1,0.1,i,5);
        }
    }
    return 0;

}

/**
void test_hq_ps(HybridQuantile hq, PSquare ps)
{

    std::string filename = "flow_301347700_intervals.txt";
    std::ifstream infile(filename);
    uint32_t curr = 0;
    std::string line;
    int lines = 0;
    while (getline(infile, line)
        //&& data.totalLines < 3000000
    )
    {
        if (line.empty()) continue;
        if (lines>1000000)break;
        try
        {
            lines ++;
            auto latency = std::stoul(line);

            if (latency > 0 && latency < 65535)
            {
                hq.update(latency);
                float m_p[Q_NUM] = {0.25, 0.5, 0.75, 0.9, 0.98};
                ps.update(latency,m_p);
            }
        } catch (const std::exception &e) {
            std::cout << "line: " <<  " transform error " << e.what() << ", content: " << line <<
                    std::endl;
            continue;
        }
    }
    std::cout << "hq"  << std::endl;
    for(int j=1; j<6; ++j) {
        std::cout << j << ": " << hq.get_quantile(j) << "|";
    }
    std::cout << std::endl;

    std::cout << "ps"  << std::endl;
    for(int j=0; j<7; ++j) {
        std::cout << j << ": " << ps.get_quantile(j) << "|";
    }
    std::cout << std::endl;
}
**/
