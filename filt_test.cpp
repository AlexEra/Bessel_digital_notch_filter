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

  std::random_device seed;
  std::array<double, 100> data;
  std::array<double, 100> f_data;

  BesselNotch2Order<double, double> bes;

  // TODO: generate noise on frequency about 70 Hz
  /* setup noise generator */
  std::mt19937 generator(seed());
  std::uniform_real_distribution<double> distribution(-max_amp_noise, max_amp_noise);

  /* preparing */
  {
    // fill data array
    double t{0};
    double noise{distribution(generator)};
    for (auto &value : data) {
      noise = distribution(generator);
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
    bes(data[i]);
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
