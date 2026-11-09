//
// Created by 13770 on 2026/4/17.
//

#ifndef PROJ___BF_H
#define PROJ___BF_H

#include <vector>
#include <cstdint>
#include <hash.h>
#include <algorithm>

/**
 * @brief bloom filter [High 4 bits | Low 4 bits]
 */
class DecayingBF {
private:

    uint32_t n;
    uint32_t m_bytes;
    std::vector<uint8_t> table;
    uint16_t fill;
    uint8_t k; // hash func


    static uint32_t hash_func(uint32_t id, int seed) {
        return Hash::BOBHash32((uint8_t*)&id, sizeof(uint32_t), seed);
    }


public:
    DecayingBF(uint32_t n, uint8_t k)
        : n(n), k(k), m_bytes((n + 1) / 2), table(m_bytes, 0) {}

    uint8_t increment(uint32_t idx) {
        uint32_t byte_idx = idx >> 1; // idx / 2
        if (idx & 1) {
            // high 4
            uint8_t val = (table[byte_idx] >> 4) & 0x0F;
            if (val < 15) table[byte_idx] += 0x10;
            return  val;
        } else {
            // low 4
            uint8_t val = table[byte_idx] & 0x0F;
            if (val < 15) table[byte_idx] += 0x01;
            return val;
        }
    }


    bool insert(uint32_t id) {
        uint8_t min = 16;
        for (uint8_t i = 0; i < k; ++i) {
            min = std::min(min,increment(hash_func(id, i*12345)%n));
        }
        if (min >= 15) {
            ++fill;
            return true;
        }
        return false;
    }



    void decay() {

        if (fill <= 1000) {
            fill = 0;
            return;
        }
        fill = 0;
        for (uint32_t i = 0; i < m_bytes; ++i) {
            uint8_t b = table[i];
            if (b == 0) continue;

            uint8_t high = (b >> 4) & 0x0F;
            uint8_t low  = b & 0x0F;

            if (high > 0) high/=2;
            if (low > 0)  low/=2;

            table[i] = (high << 4) | low;
        }
    }


};

#endif //PROJ___BF_H