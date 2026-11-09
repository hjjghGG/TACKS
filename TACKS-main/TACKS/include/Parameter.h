//
// Created by 13770 on 2026/4/28.
//

#ifndef PROJ___PARAMETER_H
#define PROJ___PARAMETER_H

#include <random>


#define Q_NUM  4
#define C_NUM  4 //psquare

#define S_NUM  8 //for simple and k-means
#define M_RATIO 0.2 // simple
#define MAX_BIN 8 //log

#define D_NUM  2 //dstream

#define alpha 0.1
#define beta 1.0
#define SEED 1919
//#define SEED 54321

//#define Stb


#define alpha_ema 1.0f
#define thresh 0.5
#define stb_para 0.1

# define DEFAULT_BASE  1.414213562

//#define HOT_CACHE

constexpr int MARKER_COUNT = 2*Q_NUM+3;
//constexpr float m_p[Q_NUM] = {0.15, 0.3, 0.45, 0.6, 0.75, 0.85, 0.9, 0.95};
//constexpr float m_p[Q_NUM] = {0.2, 0.4, 0.6, 0.75, 0.88, 0.92,0.98};
//constexpr float m_p[Q_NUM] = {0.25, 0.5, 0.75, 0.9, 0.95};
constexpr float m_p[Q_NUM] = {0.3, 0.6,  0.9, 0.95};
constexpr int STABLE_THRESHOLD = 2;


#define K_MAX_MODES  6
#define TOP  1000
#endif //PROJ___PARAMETER_H
