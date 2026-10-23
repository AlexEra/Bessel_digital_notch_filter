#include <iostream>
#include <fstream>
#include <cmath>
#include <numbers>
#include <array>
#include <random>
#include "BesselNotchFilter2Order.hpp"

using BesselNotchFilter::BesselNotch2Order;

int main() {
  const double dt{0.001}; // seconds
  const double f_sampling = 1.0 / dt; // Hz
  const double f_signal = 60; // Hz
  const double amp = 1.0;
  const double max_noise_percentage = 0.36;
  const double max_amp_noise = amp * max_noise_percentage;
  const double noise_freq = 70.0;
  const size_t max_noise_counter = 1 / (dt * noise_freq);

  std::random_device seed;
  std::array<double, 100> data;
  std::array<double, 100> f_data;

  BesselNotch2Order<double, double> bes;

  /* setup noise generator */
  std::mt19937 generator(seed());
  std::uniform_real_distribution<double> distribution(-max_amp_noise, max_amp_noise);

  /* preparing */
  {
    // fill data array
    double t{0};
    double noise{distribution(generator)};
    size_t noise_counter{0};
    for (auto &value : data) {
      if (++noise_counter == max_noise_counter) {
        noise_counter = 0;
        noise = distribution(generator);
      } else {
        noise = 0;
      }
      value = amp * sin(2 * std::numbers::pi * f_signal * t) + noise;
      t += dt;
    }
  }

  bes.setup(1000.0, 70.0, 0.69);

  // open file to save filtered data
  std::ofstream file("data.txt", std::ofstream::binary);
  // save header
  file << '#' << '\t' << 'n' << '\t' << "s_0" << '\t' << "s_1" << "\r\n";

  /* filtering */
  for (auto i{0}; i < data.size(); i++) {
    f_data[i] = bes(data[i]);
    // save data to file
    file << i << '\t' << data[i] << '\t' << f_data[i] << "\r\n";
  }
  file.close();

  /* plot */
  // plot two columns
  system("gnuplot -e \"set grid; plot 'data.txt' using 1:2 with\
    lines title 'raw', 'data.txt' using 1:3 with lines\
    title 'filtered'; pause mouse close\""
  );

  // for one column
  // system("gnuplot -e \"set size ratio 0.5; set grid; plot 'data.txt' using 1:2 with lines title 'raw'; pause mouse close\"");
  return 0;
}
