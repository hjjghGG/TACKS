#ifndef P_SQUARE_H
#define P_SQUARE_H

#include <cstdint>
#include <algorithm>
#include <cmath>
#include <Parameter.h>


class PSquare {
public:
    PSquare() : m_count(0) {
    }

    void update(uint16_t x) {
        float val = static_cast<float>(x);

        if (m_count < MARKER_COUNT) {
            m_q[m_count++] = val;
            if (m_count == MARKER_COUNT) {
                std::sort(m_q, m_q + MARKER_COUNT);
                for (int i = 0; i < MARKER_COUNT; ++i) {
                    m_n[i] = i + 1;
                }
                update_target_pos();
            }
            return;
        }

        int k;
        if (val < m_q[0]) {
            m_q[0] = val;
            k = 0;
        } else if (val >= m_q[MARKER_COUNT - 1]) {
            m_q[MARKER_COUNT - 1] = val;
            k = MARKER_COUNT - 2;
        } else {
            for (k = 0; k < MARKER_COUNT - 1; ++k) {
                if (val < m_q[k + 1]) break;
            }
        }

        for (int i = k + 1; i < MARKER_COUNT; ++i) {
            m_n[i]++;
        }
        m_count++;
        update_target_pos();

        for (int i = 1; i < MARKER_COUNT - 1; ++i) {
            float d = m_np[i] - m_n[i];

            if ((d >= 1.0f && (m_n[i + 1] - m_n[i] > 1)) ||
                (d <= -1.0f && (m_n[i - 1] - m_n[i] < -1))) {

                int d_int = (d >= 0) ? 1 : -1;
                float qt = parabolic(i, d_int);

                if (m_q[i - 1] < qt && qt < m_q[i + 1]) {
                    m_q[i] = qt;
                } else {
                    m_q[i] += d_int * (m_q[i + d_int] - m_q[i]) / (m_n[i + d_int] - m_n[i]);
                }
                m_n[i] += d_int;
            }
        }
    }



    /**
 * @brief Adjusts an existing quantile and reshapes the sketch to maintain stability.
 * * @param q_idx The index of the quantile to change (1 to Q_NUM).
 * @param new_p The new probability value (e.g., 0.85).
 * @param m_p   The updated array of all target probabilities.
 */
void adjust_quantile_and_reset(int q_idx, float new_p) {

    int target_idx = 2 * q_idx;

    if (target_idx <= 0 || target_idx >= MARKER_COUNT - 1) return;

    // Interpolate to estimate the new m_q using adjacent quantile markers
    float p_prev = (q_idx == 1) ? 0.0f : m_p[q_idx - 2];
    float p_next = (q_idx == Q_NUM) ? 1.0f : m_p[q_idx];
    float q_prev = m_q[target_idx - 2];
    float q_next = m_q[target_idx + 2];

    // Log-linear interpolation for tail-end (p > 0.5) to handle skewness
    if (new_p > 0.5f && q_prev > 0 && q_next > 0) {
        float log_q_prev = std::log(q_prev);
        float log_q_next = std::log(q_next);
        float interpolated_log_q = log_q_prev + (log_q_next - log_q_prev) * ((new_p - p_prev) / (p_next - p_prev));
        m_q[target_idx] = std::exp(interpolated_log_q);
    } else {
        // Standard linear interpolation for lower quantiles
        m_q[target_idx] = q_prev + (q_next - q_prev) * ((new_p - p_prev) / (p_next - p_prev));
    }

    //Update the global target probability array
    //m_p[q_idx - 1] = new_p;

    //Synchronize target positions (m_np) based on the new probabilities
    update_target_pos();

    // RESET of actual positions (m_n) to current desired positions (m_np)
    for (int i = 0; i < MARKER_COUNT; ++i) {
        m_n[i] = static_cast<uint16_t>(std::round(m_np[i]));
    }

    // Reset intermediate markers (the odd indices)
    for (int i = 1; i < MARKER_COUNT - 1; i += 2) {
        m_q[i] = (m_q[i-1] + m_q[i+1]) / 2.0f;
    }
}

    float get_min() const { return m_q[0]; }
    float get_max() const { return m_q[MARKER_COUNT - 1]; }
    float get_quantile(int index) const {
        return m_q[0 + 2 * index];
    }
    uint32_t get_m_count() const { return m_count; }

private:

    uint32_t m_count;
    float m_q[MARKER_COUNT]{};     // value
    uint16_t m_n[MARKER_COUNT]{};  // actual
    float m_np[MARKER_COUNT]{};    // desire

    void update_target_pos() {
        m_np[0] = 1;
        m_np[MARKER_COUNT - 1] = m_count;

        for (int i = 0; i < Q_NUM; ++i) {
            m_np[2 * (i + 1)] = 1 + (m_count - 1) * m_p[i];
        }

        for (int i = 1; i < MARKER_COUNT - 1; i += 2) {
            m_np[i] = (m_np[i - 1] + m_np[i + 1]) / 2.0f;
        }
    }

    inline float parabolic(int i, int d) {
        float n_float = static_cast<float>(m_n[i]);
        float n_prev = static_cast<float>(m_n[i - 1]);
        float n_next = static_cast<float>(m_n[i + 1]);

        return m_q[i] + (d / (n_next - n_prev)) * (
            (n_float - n_prev + d) * (m_q[i + 1] - m_q[i]) / (n_next - n_float) +
            (n_next - n_float - d) * (m_q[i] - m_q[i - 1]) / (n_float - n_prev)
        );
    }
};

#endif