#pragma once

/**
 * *_cont variables for calculations of the filter's coefficients
    these variables are connected with analog lowpass Bessel filter with cutoff frequency 1 rad/s
 */
constexpr float b_0_cont = 1.0;
constexpr float a_cont[] = {1.0, 1.7320508075688772, 1.0};

/**
 * @brief Class of the digital Bessel filter
    Notch type, only 2th order (continious, 5th order with discrete)
 */
class BesselNotchFilter2Order {
    public:
        BesselNotchFilter2Order();
        BesselNotchFilter2Order(float freq_sample, float freq_medium, float quality_factor);
        ~BesselNotchFilter2Order();
        float step(float new_val);
        void set_quality_factor(float new_q_f);
        float get_quality_factor(void);
        void set_medium_frequency(float new_f_m);
        float get_medium_frequency(void);
        void set_sample_frequency(float new_f_s);
        float get_sample_frequency(void);
    protected:
        void coefficients_calculating(void);
        float f_s; // sample frequency
        float f_m; // medium frequency
        float q_f; // quality factor
        float a[5], b[5]; // filter coefficients
        float y_prev[4], x_prev[4]; // arrays for keeping previous x (input) and y (output) values
};
