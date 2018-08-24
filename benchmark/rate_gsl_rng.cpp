#include <iostream>
#include <cstddef>
#include <chrono>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// types of distributions
enum class dist {uniform, normal_boxm, normal_zigg, exp};

template <dist P>
void
calc_gsl_rng_rate(const std::size_t samples, unsigned long seed)
{
  // setup GSL RNG
  gsl_rng_env_setup();
  const gsl_rng_type* T = gsl_rng_default;
  gsl_rng* r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);

  // calculate the rate of random numbers generated per second
  // also calculate the mean while generating the numbers so that
  // compiler doesn't remove the sampling loop during optimization
  double mean = 0;
  using time_pt = std::chrono::system_clock::time_point;
  time_pt start = std::chrono::system_clock::now();
  for (std::size_t i = 0; i < samples; i++) {
    switch(P)
    {
      case dist::uniform     : mean += gsl_rng_uniform(r); break;
      case dist::normal_boxm : mean += gsl_ran_gaussian(r, 1.); break;
      case dist::normal_zigg : mean += gsl_ran_gaussian_ziggurat(r, 1.); break;
      case dist::exp         : mean += gsl_ran_exponential(r, 1.); break;
      default                : break;
    }
  }
  time_pt end = std::chrono::system_clock::now();
  mean /= static_cast<double>(samples);

  // print results
  std::string rng_name = std::string("gsl_") + std::string(gsl_rng_name(r));
  std::string dist_name;
  switch(P)
  {
    case dist::uniform     : dist_name = "uniform"; break;
    case dist::normal_boxm : dist_name = "normal_boxm"; break;
    case dist::normal_zigg : dist_name = "normal_zigg"; break;
    case dist::exp         : dist_name = "exponential"; break;
    default                : break;
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

  // free GSL RNG
  gsl_rng_free(r);
}

int
main()
{
  unsigned long int seed    = 1234;
  const std::size_t samples = 1e9;

  calc_gsl_rng_rate<dist::uniform>(samples, seed);
  calc_gsl_rng_rate<dist::normal_boxm>(samples, seed);
  calc_gsl_rng_rate<dist::normal_zigg>(samples, seed);
  calc_gsl_rng_rate<dist::exp>(samples, seed);

  return 0;
}
