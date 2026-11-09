//
// Created by 13770 on 2026/5/7.
//

#ifndef PROJ___TRACEDATA_H
#define PROJ___TRACEDATA_H
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <iomanip>
#include <filesystem>

struct TraceData {
    std::string filename;
    int totalLines = 0;
    int nonZerosLatency = 0;
    std::unordered_map<uint32_t, std::vector<uint32_t> > flowLatencies;
    std::unordered_map<uint32_t, uint32_t> flowCounters;
    std::vector<std::pair<uint32_t, uint32_t> > allPackets;
};


    // Read trace data file
    static TraceData readTraceFile(const std::string &filename) {
        TraceData data;
        data.filename = filename;
        std::ifstream infile(filename);
        if (!infile.is_open()) {
            std::cout << "can't open file " << filename << std::endl;
            return data;
        }
        uint32_t curr = 0;
        std::string line;
        while (getline(infile, line)
            // data.totalLines < 3000000
        ) {
            if (line.empty()) continue;
            data.totalLines++;

            if (data.totalLines % 1000000 == 0) {
                std::cout << "processed " << data.totalLines << " lines" << std::endl;
            }

            std::vector<std::string> splits;
            std::string split;
            std::istringstream iss(line);
            char delim = ',';

            if (filename.find("case") != std::string::npos ||
                filename.find("MIRAI") != std::string::npos) delim = ' ';

            while (std::getline(iss, split, delim)) {
                splits.push_back(split);
            }

            if (splits.size() < 2) continue;

            try {
                auto ID = (uint32_t) stoull(splits[0]);
                auto latency = std::stoul(splits[1]);

                if (latency < 100000) {
                    if (filename.find("case") != std::string::npos) {
                        latency /= 1000;
                    }
                    curr += latency;
                    latency = curr;

                    data.nonZerosLatency++;
                    data.allPackets.emplace_back(ID, latency);

                    data.flowLatencies[ID].push_back(latency);

                    data.flowCounters[ID]++;
                }
            } catch (const std::exception &e) {
                std::cout << "line: " << data.totalLines << " transform error " << e.what() << ", content: " << line <<
                        std::endl;
                continue;
            }
        }

        infile.close();
        std::cout << "Read " << data.totalLines << " lines, Non-zero: " << data.nonZerosLatency << std::endl;

        return data;
    }
#endif //PROJ___TRACEDATA_H