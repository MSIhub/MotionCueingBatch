#ifndef RESAMPLER_H
#define RESAMPLER_H

#include <cstdint>

//
template <typename T> class Resampler {
public:
  Resampler() : Resampler(T(0.9)) {}
  Resampler(T beta_) : beta{beta_} {}

  void addSample(T sample, T time)
  {
    // printf(" idx: %d \n", idx);
    ++idx %= bufSize;
    sampleBuf[idx] = sample;
    tBuf[idx] = time;
  }

  T filt(T t)
  {
    auto prev = idx > 0 ? idx - 1 : bufSize - 1;
    auto sample = (tBuf[idx] == tBuf[prev])
                      ? sampleBuf[prev]
                      : (sampleBuf[prev]) + (sampleBuf[idx] - sampleBuf[prev]) /
                                                (tBuf[idx] - tBuf[prev]) *
                                                (t - tBuf[prev]);
    smoothData = smoothData - (beta * (smoothData - sample));
    return smoothData;
  }

public:
  static constexpr int bufSize{3};
  const T beta;
  T sampleBuf[bufSize]{T(0.0)};
  T tBuf[bufSize]{T(0.0)};
  uint32_t idx{0};
  T smoothData{T()};
  int prev;
};
//
#endif