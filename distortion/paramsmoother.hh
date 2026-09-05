#pragma once

#include <cmath>
#include <cassert>
#include <algorithm>

enum class SmootherType {
  linear,
  logarithmic
};

template <SmootherType smoother_type>
class ParamSmoother
{
public:
  ParamSmoother (float start_value)
  {
    set_target (start_value, true);
  }

  /**
   * Initializes the ramp length in samples, with a minimum of one sample.
   */
  void reset (double sample_rate, float ramp_time_sec)
  {
    assert (sample_rate > 0);
    ramp_samples_ = static_cast<unsigned int> (std::max (ramp_time_sec * sample_rate, 1.0));
  }

  /**
   * Sets a new destination target for the parameter.
   */
  void set_target (float target, bool now = false)
  {
    if constexpr (smoother_type == SmootherType::logarithmic)
      assert (target > 0.0f);

    if (target == target_ && !now)
      return;

    if (target == current_)
      now = true;

    target_ = target;

    if (now)
      {
        current_ = target_;
        if constexpr (smoother_type == SmootherType::linear)
          {
            step_or_factor_ = 0.0f;
          }
        else
          {
            step_or_factor_ = 1.0f;
          }
        ramp_counter_ = 0;
        return;
      }

    // user must initialize ramp length using reset() before using set_target (..., false)
    assert (ramp_samples_ > 0);
    ramp_counter_ = ramp_samples_;

    if constexpr (smoother_type == SmootherType::linear)
      {
        step_or_factor_ = (target_ - current_) / ramp_samples_;
      }
    else
      {
        assert (current_ > 0.0f);
        step_or_factor_ = std::pow (target_ / current_, 1.0f / ramp_samples_);
      }
  }

  /**
   * Get next smoothed value
   */
  inline float
  get_next()
  {
    if (ramp_counter_ > 0)
      {
        if constexpr (smoother_type == SmootherType::linear)
          {
            current_ += step_or_factor_;
          }
        else
          {
            current_ *= step_or_factor_;
          }

        ramp_counter_--;

        if (ramp_counter_ == 0)
          current_ = target_;
      }
    return current_;
  }
  /**
   * Check wether we're currently smoothing or if the output is a constant value
   */
  bool
  is_constant()
  {
    return ramp_counter_ == 0;
  }

  /**
   * Fills an entire buffer array block.
   */
  void
  process_block (float* output_buffer, unsigned int num_samples)
  {
    for (unsigned int i = 0; i < num_samples; ++i)
      output_buffer[i] = get_next();
  }

private:
  unsigned int ramp_samples_ = 0;

  float target_ = 0.0f;
  float current_ = 0.0f;
  float step_or_factor_ = 0.0f;
  unsigned int ramp_counter_ = 0;
};
