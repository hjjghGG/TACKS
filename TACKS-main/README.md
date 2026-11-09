# Build Guide

## Introduction

OASIS leverages a synergistic two-tiered pipeline: a dynamic filtering stage with an optimized bucket-cell architecture first screens for significant flows, followed by a log-binned center array structure that tracks their hierarchical stability. To ensure robustness under dynamic traffic, we introduce adaptive stability calibration to counteract cluster evolution, and an optional inline Coefficient of Variation estimation module to profile granular temporal volatility.

## Prerequisites

Before building, ensure your system meets the following requirements:

- **CMake** ≥ 4.0
- **C++ Compiler**: GCC ≥ 15.2.0
- **Build Tools**: make, ninja, or Clion

## Quick Start


```bash
# 1. Create and enter build directory
mkdir build && cd build

# 2. Configure the project
cmake ..

# 3. Compile with all available cores
cmake --build . -j $(nproc)

```

Initialize each benchmark class in main.cpp to test the corresponding algorithm:

* **BenchMark**: Corresponds to PSF
* **BenchMark2**: Corresponds to OASIS
* **BenchMark3**: Corresponds to Baseline
* **BenchMark4**: Corresponds to DS-Sketch
* **BenchMarkKM**: Corresponds to KMeans+
