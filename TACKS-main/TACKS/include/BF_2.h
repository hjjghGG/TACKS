//
// Created by 13770 on 2026/4/17.
//

#ifndef PROJ___BF_2_H
#define PROJ___BF_2_H

#include <vector>
#include <cstdint>
#include <hash.h>
#include <algorithm>

/**
 * @brief Bloom filter with Flag bit per counter
 *
 * Each 4-bit counter is reorganized as: [Flag bit (1 bit) | Counter value (3 bits)]
 * - Flag bit: 1 bit, used to store a flag F (0 or 1)
 * - Counter: 3 bits, value range 0-7 (instead of original 0-15)
 *
 */
class DecayingBF {
private:

    uint32_t n;
    uint32_t m_bytes;
    std::vector<uint8_t> table;
    uint16_t fill;
    uint8_t k; // number of hash functions

    static uint32_t hash_func(uint32_t id, int seed) {
        return Hash::BOBHash32((uint8_t*)&id, sizeof(uint32_t), seed);
    }

    /**
     * @brief Get flag and counter from a byte at given index
     * @param idx bit index (0-based, each 4-bit cell)
     * @param flag output: flag bit (0 or 1)
     * @param counter output: counter value (0-7)
     */
    void get_cell(uint32_t idx, uint8_t& flag, uint8_t& counter) const {
        uint32_t byte_idx = idx >> 1; // idx / 2
        if (idx & 1) {
            // High 4 bits: [Flag (bit7) | Counter (bits6-4)]
            uint8_t cell = (table[byte_idx] >> 4) & 0x0F;
            flag = (cell >> 3) & 0x01;
            counter = cell & 0x07;
        } else {
            // Low 4 bits: [Flag (bit3) | Counter (bits2-0)]
            uint8_t cell = table[byte_idx] & 0x0F;
            flag = (cell >> 3) & 0x01;
            counter = cell & 0x07;
        }
    }

    /**
     * @brief Set flag and counter for a given cell
     * @param idx bit index (0-based, each 4-bit cell)
     * @param flag flag bit (0 or 1)
     * @param counter counter value (0-7)
     */
    void set_cell(uint32_t idx, uint8_t flag, uint8_t counter) {
        uint32_t byte_idx = idx >> 1; // idx / 2
        uint8_t cell = (flag << 3) | (counter & 0x07);

        if (idx & 1) {
            // High 4 bits
            uint8_t low = table[byte_idx] & 0x0F;
            table[byte_idx] = (cell << 4) | low;
        } else {
            // Low 4 bits
            uint8_t high = (table[byte_idx] >> 4) & 0x0F;
            table[byte_idx] = (high << 4) | cell;
        }
    }

    /**
     * @brief Increment counter for a cell based on flag F
     * @param idx bit index
     * @param F input flag (0 or 1)
     * @return counter value before increment (or current value if not incremented)
     */
    uint8_t increment_with_flag(uint32_t idx, uint8_t F) {
        uint8_t old_flag, old_counter;
        get_cell(idx, old_flag, old_counter);

        if (old_flag == F) {
            // Same flag: do not increment counter
            return old_counter;
        } else {
            // Different flag: increment counter and update flag
            uint8_t new_counter = old_counter;
            if (new_counter < 7) {  // Max counter value is 7 for 3 bits
                new_counter++;
            }
            set_cell(idx, F, new_counter);
            return old_counter;
        }
    }

public:
    DecayingBF(uint32_t n, uint8_t k)
        : n(n), k(k), m_bytes((n + 1) / 2), table(m_bytes, 0), fill(0) {}

    /**
     * @brief Insert an element with flag F
     * @param id element ID to insert
     * @param F flag bit (0 or 1)
     * @return true if all k counters reached threshold (7), false otherwise
     */
    bool insert(uint32_t id, uint8_t F) {
        uint8_t min_counter = 8;
        for (uint8_t i = 0; i < k; ++i) {
            uint32_t idx = hash_func(id, i * 12345) % n;
            uint8_t old_val = increment_with_flag(idx, F & 0x01);  // Ensure F is 0 or 1
            min_counter = std::min(min_counter, old_val);
        }
        if (min_counter >= 7) {  // Threshold is 7 (max value for 3 bits)
            ++fill;
            return true;
        }
        return false;
    }

    /**
     * @brief Check if an element might exist (without modifying counters)
     * @param id element ID to check
     * @param F flag bit to check against
     * @return true if all counters have value >= 1 AND flag matches, false otherwise
     */
    bool might_contain(uint32_t id, uint8_t F) const {
        for (uint8_t i = 0; i < k; ++i) {
            uint32_t idx = hash_func(id, i * 12345) % n;
            uint8_t flag, counter;
            get_cell(idx, flag, counter);
            if (counter == 0 || flag != (F & 0x01)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Decay all counters (divide by 2)
     * Flag bits remain unchanged during decay
     */
    void decay() {
        if (fill <= 1000) {
            fill = 0;
            return;
        }
        fill = 0;
        for (uint32_t i = 0; i < m_bytes; ++i) {
            uint8_t b = table[i];
            if (b == 0) continue;

            // Process high 4 bits
            uint8_t high_cell = (b >> 4) & 0x0F;
            if (high_cell != 0) {
                uint8_t high_flag = (high_cell >> 3) & 0x01;
                uint8_t high_counter = high_cell & 0x07;
                if (high_counter > 0) high_counter = high_counter / 2;
                high_cell = (high_flag << 3) | high_counter;
            }

            // Process low 4 bits
            uint8_t low_cell = b & 0x0F;
            if (low_cell != 0) {
                uint8_t low_flag = (low_cell >> 3) & 0x01;
                uint8_t low_counter = low_cell & 0x07;
                if (low_counter > 0) low_counter = low_counter / 2;
                low_cell = (low_flag << 3) | low_counter;
            }

            table[i] = (high_cell << 4) | low_cell;
        }
    }

    /**
     * @brief Reset all counters to 0 (preserves flag bits? use reset_all for full reset)
     */
    void reset_counters() {
        for (uint32_t i = 0; i < m_bytes; ++i) {
            uint8_t b = table[i];
            uint8_t high_flag = (b >> 7) & 0x01;      // Flag from high cell (bit7)
            uint8_t low_flag = (b >> 3) & 0x01;       // Flag from low cell (bit3)

            uint8_t high_cell = (high_flag << 3) | 0;  // Counter = 0
            uint8_t low_cell = (low_flag << 3) | 0;    // Counter = 0

            table[i] = (high_cell << 4) | low_cell;
        }
    }

    /**
     * @brief Full reset: clear both flags and counters
     */
    void full_reset() {
        std::fill(table.begin(), table.end(), 0);
        fill = 0;
    }
};

#endif //PROJ___BF_2_H