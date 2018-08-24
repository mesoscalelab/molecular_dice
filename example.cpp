#include <iostream>
#include <md_rng.h>

int main()
{
  // create an instance of Molecular Dice RNG
  md::rng r;

  // generate a uniformly distributed random number
  std::cout << r.uniform() << std::endl;

  // generate a normally distributed random number
  std::cout << r.normal() << std::endl;

  // generate an exponentially distributed random number
  std::cout << r.exp() << std::endl;

  return 0;
}
