#include <math.h>
#include <stdint.h>
#include "BesselNotchFilter2Order.h"


/**
 * @brief Creates basic filter with cutoff frequency 1 rad/s and sample frequency 1 Hz
 */
BesselNotchFilter2Order::BesselNotchFilter2Order() : f_s{1.0}, f_m{2 * M_PI}, q_f{1.0} {}

/**
 * @brief Constructs filter object
 * @param freq_sample sample frequency in Hz
 * @param freq_medium filter frequency in Hz (is automatically transformed to rad/s), frequency to be rejected
 * @param quality_factor quality factor from 0.0 to 1.0, width of the stopband
 */
BesselNotchFilter2Order::BesselNotchFilter2Order(float freq_sample, float freq_medium, float quality_factor) : f_s{freq_sample} {
    f_m = 2 * M_PI * freq_medium;
    q_f = quality_factor;
    coefficients_calculating();
}

BesselNotchFilter2Order::~BesselNotchFilter2Order() {}

/**
 * @brief Used to filtrate for one step in accordance to sample frequency
 * @param new_val new filter input value
 * @return filtered value y
 */
float BesselNotchFilter2Order::step(float new_val) {
    float y = new_val * b[0] + x_prev[0] * b[1] + x_prev[1] * b[2] + 
            x_prev[2] * b[3] + x_prev[3] * b[4] - y_prev[0] * a[1] - 
            y_prev[1] * a[2] - y_prev[2] * a[3] - y_prev[3] * a[4];

    /* delay for the x and y values */
    for (uint8_t i = (sizeof(x_prev) / sizeof(float)) - 1; i > 0; i--) {
        x_prev[i] = x_prev[i - 1];
        y_prev[i] = y_prev[i - 1];
    }
    x_prev[0] = new_val;
    y_prev[0] = y;

    return y;
}

/**
 * @param new_q_f new quality factor value
 */
void BesselNotchFilter2Order::set_quality_factor(float new_q_f) {
    q_f = 1 / new_q_f;
    coefficients_calculating();
}

/**
 * @return the inverse quality factor value (1/Q)
 */
float BesselNotchFilter2Order::get_quality_factor(void) {
    return q_f;
}

/**
 * @param new_f_m new medium frequency value
 */
void BesselNotchFilter2Order::set_medium_frequency(float new_f_m) {
    f_m = 2 * M_PI * new_f_m; // transforming to the rad/s
    coefficients_calculating();
}

/**
 * @return medium frequency value in rad/s
 */
float BesselNotchFilter2Order::get_medium_frequency(void) {
    return f_m;
}

/**
 * @param new_f_s the new sample frequency value
 */
void BesselNotchFilter2Order::set_sample_frequency(float new_f_s) {
    f_s = new_f_s;
    coefficients_calculating();
}

/**
 * @return the sample frequency value
 */
float BesselNotchFilter2Order::get_sample_frequency(void) {
    return f_s;
}

/**
 * @brief Calculating the filter coefficients
 * Firstly, the analog 2nd order filter coefficients according to medium frequency are obtained,
 * Finally, the digital filter coefficients are calculated
 * Coefficients equations were obtained with bilinear transformation of the analog filter transfer function
 */
void BesselNotchFilter2Order::coefficients_calculating(void) {
    float a_tmp[5], b_tmp[5]; // temporary values to keep coefficients
    float a_0; // for divide operation
    float f_2, f_3, f_4; // for keeping the powers of the frequency

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
