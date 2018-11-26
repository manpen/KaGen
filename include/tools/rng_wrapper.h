/*******************************************************************************
 * include/tools/rng_wrapper.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _RNG_WRAPPER_H_
#define _RNG_WRAPPER_H_

#include <random>

#include "methodR.hpp"

namespace kagen {

class RNGWrapperImpl {
 public:
  RNGWrapperImpl(const PGeneratorConfig &config)
      : config_(config)
  {};

  HPFloat GenerateHypergeometric(SInt seed, HPFloat n, HPFloat m, HPFloat N) {
    if (config_.use_binom) {
      return GenerateBinomial(seed, n, (double) (m / N));
    } else {
      if (m < 1) return 0;
      hyp_.seed(seed);
      return hyp_(n, N-n, m);
    }
  }

  SInt GenerateBinomial(SInt seed, SInt n, double p) {
    rng_.seed(seed);
    std::binomial_distribution<SInt> bin(n, p);
    return bin(rng_);
  }

  template <typename F>
  void GenerateSample(SInt seed, HPFloat N, SInt n, F &&callback) {
    sampling::HashSampling<> hs(seed, config_.base_size);
    sampling::SeqDivideSampling<> sds(hs, config_.base_size, seed, config_.use_binom);
    sds.sample(N, n, callback);
  }

 private:
  const PGeneratorConfig &config_;

  std::mt19937_64 rng_;
  sampling::hypergeometric_distribution<> hyp_;
};

using RNGWrapper = RNGWrapperImpl;

}
#endif
