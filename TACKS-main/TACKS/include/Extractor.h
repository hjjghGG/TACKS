#ifndef PROJ_LAYER3_H
#define PROJ_LAYER3_H

#include <cstdint>
#include <vector>
#include <cmath>
#include <Parameter.h>
#include <numeric>
#define EMA

struct List3
{
    uint16_t centers[C_NUM]{};
    uint16_t counts[C_NUM]{};
    uint16_t cnt[C_NUM]{};
    float sum_sq_diff[C_NUM]{};

# ifndef EMA
    int16_t accumulate[C_NUM]{};
    uint8_t accumulatecount[C_NUM]{};
# endif

    uint8_t stable;

    List3() :  stable(0){}


    void update(uint16_t latency)
    {
        int best_idx = -1;
        float min_diff = -1.0f;
        for (int i = 0; i < C_NUM; ++i)
        {
            if (centers[i] != 0)
            {
                float diff = std::abs(static_cast<float>(latency) - centers[i]);
                if (best_idx == -1 || diff < min_diff)
                {
                    min_diff = diff;
                    best_idx = i;
                }
            }
        }

        if (best_idx != -1)
        {
            cnt[best_idx]++;

# ifndef EMA
            if (stable > STABLE_THRESHOLD)
            {
                int16_t offset = static_cast<int16_t>(latency) - static_cast<int16_t>(centers[best_idx]);

                if (accumulatecount[best_idx] >= 250 ||
                    accumulate[best_idx] > 30000 ||
                    accumulate[best_idx] < -30000)
                { //overflow
                    apply_correction(best_idx);
                }
                if (std::abs(offset) < centers[best_idx]*0.5)
                {
                    accumulate[best_idx] += offset;
                    accumulatecount[best_idx]++;
                }
            }

#else
            uint16_t old_center = centers[best_idx];

            float times = (double)latency/(double)old_center;
            if (times <1)times = 1.0/times;
            if (times < thresh) {
                auto proposed_center =(double)old_center +(static_cast<float>(latency) - old_center )
                *( alpha_ema/(alpha_ema + counts[best_idx]));

                // check bound
                if (best_idx > 0 && centers[best_idx - 1] != 0) {
                    if (proposed_center <= centers[best_idx - 1]) proposed_center = centers[best_idx - 1] + 1;
                }
                if (best_idx < C_NUM-1 && centers[best_idx + 1] != 0) {
                    if (proposed_center >= centers[best_idx + 1]) proposed_center = centers[best_idx + 1] - 1;
                }

                centers[best_idx] = proposed_center;
            }


#endif
            {
                float d = static_cast<float>(latency) - centers[best_idx];
                sum_sq_diff[best_idx] += d * d;
            }

            if ((min_diff / centers[best_idx]) < alpha || min_diff < beta)
            {
                counts[best_idx]++;
                return;
            }

        }

        //  nothing
    }

#ifndef EMA
    void apply_correction(int idx)
    {

        if (accumulatecount[idx] == 0) return;

        int16_t correction = accumulate[idx] / accumulatecount[idx];
        if (correction == 0) {
            accumulate[idx] = 0;
            accumulatecount[idx] = 0;
            return;
        }

        uint16_t old_center = centers[idx];
        uint16_t proposed_center = old_center + correction;

        // check bound
        if (idx > 0 && centers[idx - 1] != 0) {
            if (proposed_center <= centers[idx - 1]) proposed_center = centers[idx - 1] + 1;
        }
        if (idx < C_NUM-1 && centers[idx + 1] != 0) {
            if (proposed_center >= centers[idx + 1]) proposed_center = centers[idx + 1] - 1;
        }

        centers[idx] = proposed_center;

        // reset
        accumulate[idx] = 0;
        accumulatecount[idx] = 0;
    }
#endif


    void set_centers(const std::vector<uint16_t>& new_centers)
    {
        if (stable > STABLE_THRESHOLD) return;
        bool stb = std::accumulate(std::begin(cnt),std::end(cnt),0)> 50;
        for (size_t i = 0; i < C_NUM && i < new_centers.size(); ++i)
        {
            float diff = std::abs((double)new_centers[i] - (double)centers[i])
                / (double)centers[i];
            stb &= diff < 0.2;
            centers[i] = static_cast<uint16_t>(new_centers[i]);
            if (diff >= 0.2)
            {
                counts[i] = std::round(counts[i] * std::max((float)0.5,1-diff));
                cnt[i] -= std::round(counts[i] * std::min((float)0.5,diff));
            }
        }
        // set stable
        if (stb)
        {
            stable++;
        }
    }

    void output() const
    {

        std::cout << std::left << std::setw(12)
                       << std::setw(10) << ( stable> STABLE_THRESHOLD ? "Stable" : "Active") <<"|";

        for (int j = 0; j < C_NUM; ++j)
        {
            std::cout << "[" << centers[j] << ":" << counts[j] << "] ";
        }
        for (int j = 0; j < C_NUM; ++j)
        {
            std::cout << "[" << cnt[j] << "] ";
        }
        std::cout << std::endl;
        std::cout << std::string(60, '-') << std::endl;

    }
};

#endif //PROJ_LAYER3_H
