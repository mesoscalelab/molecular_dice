#include <iostream>
#include <cstddef>
#include <chrono>
#include <md_rng.h>

// types of distributions
enum class dist {uniform, normal, exp};

template <dist P>
void
calc_md_rng_rate(const std::size_t samples, unsigned long seed)
{
  // setup molecular dice RNG
  md::rng r(seed);

  // calculate the rate of random numbers generated per second
  // also calculate the mean while generating the numbers so that
  // compiler doesn't remove the sampling loop during optimization
  double mean = 0;
  using time_pt = std::chrono::system_clock::time_point;
  time_pt start = std::chrono::system_clock::now();
  for (std::size_t i = 0; i < samples; i++) {
    switch(P)
    {
      case dist::uniform : mean += r.uniform(); break;
      case dist::normal  : mean += r.normal(); break;
      case dist::exp     : mean += r.exp(); break;
      default            : break;
    }
  }
  time_pt end = std::chrono::system_clock::now();
  mean /= static_cast<double>(samples);

  // print results
  std::string rng_name = std::string("molecular_dice");
  std::string dist_name;
  switch(P)
  {
    case dist::uniform : dist_name = "uniform"; break;
    case dist::normal  : dist_name = "normal"; break;
    case dist::exp     : dist_name = "exponential"; break;
    default            : break;
  }
  std::chrono::duration<double> time_taken = end - start;
  const double rate = samples / time_taken.count();
  std::cout << std::scientific;
  std::cout << rng_name << ",";
  std::cout << dist_name << ",";
  std::cout << rate << ",";
  std::cout << mean << ",";
  std::cout << static_cast<double>(samples) << std::endl;
  std::cout << std::defaultfloat;

  return;
}

int
main()
{
  unsigned long int seed    = 1234;
  const std::size_t samples = 1e9;

  calc_md_rng_rate<dist::uniform>(samples, seed);
  calc_md_rng_rate<dist::normal>(samples, seed);
  calc_md_rng_rate<dist::exp>(samples, seed);

  return 0;
}
