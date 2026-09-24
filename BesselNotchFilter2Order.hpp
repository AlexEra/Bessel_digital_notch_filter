#pragma once

#include <numbers> // for pi
#include <concepts>

namespace BesselNotchFilter {

/**
 * *_cont variables for calculations of the filter's coefficients
    these variables are connected with analog lowpass Bessel filter with cutoff frequency 1 rad/s
 */

constexpr float b_0_cont = 1.0;
constexpr float a_cont[] = {1.0, 1.7320508075688772, 1.0};

template <typename T>
concept vals_to_filter = requires (T value) {
  value = value;
  value += value;
  value -= value;
  value + value;
  value - value;
  value / value;
  value * value;
  value < value;
  value > value;
  value <= value;
  value >= value;
};

template <typename F>
concept float_or_double = std::same_as<F, float> || std::same_as<F, double>;

/**
 * @brief Class of the digital Bessel filter
 * 
 * Notch type, only 2th order (continious, 5th order with discrete)
 */
template<vals_to_filter T, float_or_double F>
class BesselNotch2Order {
public:
  /**
   * @brief Used to filtrate for one step in accordance to sample frequency
   * @param new_val New filter input value
   * @return Filtered value
   */
  T step(T new_val) {
    T y = new_val * b[0] + x_prev[0] * b[1] + x_prev[1] * b[2] + 
      x_prev[2] * b[3] + x_prev[3] * b[4] - y_prev[0] * a[1] - 
      y_prev[1] * a[2] - y_prev[2] * a[3] - y_prev[3] * a[4];

    /* delay for the x and y values */
    for (unsigned char i{(sizeof(x_prev) / sizeof(T)) - 1}; i > 0; i--) {
      x_prev[i] = x_prev[i - 1];
      y_prev[i] = y_prev[i - 1];
    }
    x_prev[0] = new_val;
    y_prev[0] = y;

    return y;
  }

  /**
   * @param new_q_f New quality factor value
   */
  void set_quality_factor(F new_q_f) {
    q_f = 1 / new_q_f;
    coefficients_calculating();
  }

  /**
   * @return Inverse quality factor value (1 / Q)
   */
  F get_quality_factor(void) {
    return q_f;
  }

  /**
   * @param new_f_m New medium frequency value
   */
  void set_medium_frequency(F new_f_m) {
    f_m = 2 * std::numbers::pi * new_f_m; // transforming to the rad/s
    coefficients_calculating();
  }

  /**
   * @return Medium frequency value in rad/s
   */
  F get_medium_frequency(void) {
    return f_m;
  }

  /**
   * @param new_f_s New sample frequency value
   */
  void set_sample_frequency(F new_f_s) {
    f_s = new_f_s;
    coefficients_calculating();
  }

  /**
   * @return Sample frequency value
   */
  F get_sample_frequency(void) {
    return f_s;
  }

  /**
   * @brief Setup new filter parameters
   * @param new_f_s New sample frequency value [Hz]
   * @param new_f_m New medium frequency [Hz]
   * @param new_q_f New quality factor (1 / Q)
   */
  void setup(F new_f_s, F new_f_m, F new_q_f) {
    f_s = new_f_s;
    f_m = 2 * std::numbers::pi * new_f_m; // transforming to the rad/s
    q_f = 1 / new_q_f;
    coefficients_calculating();
  }

  T operator() (T new_value) {
    return step(new_value);
  }

protected:
  /**
   * @brief Calculating the filter coefficients
   * 
   * Firstly, the analog 2nd order filter coefficients according to medium frequency are obtained.
   * 
   * Finally, the digital filter coefficients are calculated.
   * 
   * Coefficients equations were obtained with bilinear transformation of the analog filter transfer function.
   */
  void coefficients_calculating(void) {
    F a_tmp[5], b_tmp[5]; // temporary values to keep coefficients
    F a_0; // for divide operation
    F f_2, f_3, f_4; // for keeping the powers of the frequency

    /* temporary powers of the medium frequency */
    f_2 = pow(f_m, 2);
    f_3 = pow(f_m, 3);
    f_4 = pow(f_m, 4);

    /* calculating the analog notch filter coefficients */
    a_0 = a_cont[2] * pow(q_f, 2) / f_4;
    a_tmp[0] = 1.0;
    a_tmp[1] = (a_cont[1] * q_f / f_3) / a_0;
    a_tmp[2] = ((2 * a_cont[2] * pow(q_f, 2) + a_cont[0]) / f_2) / a_0;
    a_tmp[3] = (a_cont[1] * q_f / f_m) / a_0;
    a_tmp[4] = (a_cont[2] * pow(q_f, 2)) / a_0;

    b_tmp[0] = (b_0_cont * pow(q_f, 2) / f_4) / a_0;
    b_tmp[1] = 0.0;
    b_tmp[2] = (2 * b_0_cont * pow(q_f, 2) / f_2) / a_0;
    b_tmp[3] = 0.0;
    b_tmp[4] = (b_0_cont * pow(q_f, 2)) / a_0;

    /* temporary powers of the sample frequency */
    f_2 = pow(f_s, 2);
    f_3 = pow(f_s, 3);
    f_4 = pow(f_s, 4);

    a_0 = 16 * a_tmp[0] * f_4 + 8 * a_tmp[1] * f_3 + 4 * a_tmp[2] * f_2 + 2 * a_tmp[3] * f_s + a_tmp[4];
    a[0] = 1.0;
    a[1] = (-64 * a_tmp[0] * f_4 - 16 * a_tmp[1] * f_3 + 4 * a_tmp[3] * f_s + 4 * a_tmp[4]) / a_0;
    a[2] = (96 * a_tmp[0] * f_4 - 8 * a_tmp[2] * f_2 + 6 * a_tmp[4]) / a_0;
    a[3] = (-64 * a_tmp[0] * f_4 + 16 * a_tmp[1] * f_3 - 4 * a_tmp[3] * f_s + 4 * a_tmp[4]) / a_0;
    a[4] = (16 * a_tmp[0] * f_4 - 8 * a_tmp[1] * f_3 + 4 * a_tmp[2] * f_2 - 2 * a_tmp[3] * f_s + a_tmp[4]) / a_0;

    b[0] = (16 * b_tmp[0] * f_4 + 4 * b_tmp[2] * f_2 + b_tmp[4]) / a_0;
    b[1] = (-64 * b_tmp[0] * f_4 + 4 * b_tmp[4]) / a_0;
    b[2] = (96 * b_tmp[0] * f_4 - 8 * b_tmp[2] * f_2 + 6 * b_tmp[4]) / a_0;
    b[3] = b[1];
    b[4] = b[0];
  }

protected:
  F f_s{1.0}; // sample frequency
  F f_m{2.0 * std::numbers::pi}; // medium frequency
  F q_f{1.0}; // quality factor
  F a[5], b[5]; // filter coefficients
  T y_prev[4] = {0,}; // for keeping output values with delay
  T x_prev[4] = {0,}; // for keeping input values with delay
};

} /* namespace BesselNotchFilter */
