Molecular Dice: Random Number Generator
=======================================
Molecular Dice (MD) random number generator (RNG) samples the
position, velocity and energy properties of particles at
equilibrium in order to generate uniform, normal and exponential
random variates, respectively. This package contains a header-only
C++ implementation of the molecular dice RNG.

Usage
=====
`example.cpp` demonstrates the simple manner in which an instance
of the RNG can be created and used to generate the random variates
following different distributions.

Compile and Run
===============
```
$ g++ -std=c++11 -O3 -march=native -I include/ example.cpp -o example
$ ./example
```

Benchmarks
==========
The rate of random number generation for the above distributions using
molecular dice is found to be faster than currently available implementations
of RNGs in the standard C++ `random` module as well as the GNU Scientific
Library (GSL) RNG module. The driver code for performance comparison is
contained in the `benchmark` folder. The benchmark code for molecular
dice RNG can be compiled and run using:
```
$ g++ -std=c++11 -O3 -march=native -I ../include/ rate_md_rng.cpp -o rate_md
$ ./rate_md
```

The benchmark code for GSL RNG can be compiled and run using:
```
$ g++ -std=c++11 -O3 -march=native -I ../include/ rate_gsl_rng.cpp -o rate_gsl -lgsl -lgslcblas
$ ./rate_gsl
```

These benchmark codes print the rate of double precision random number
generation per second as well as the mean of the generated numbers for each
distribution, along with the type of RNG and number of samples used for
the calculation - all in a comma-separated value format.

For example, on a machine with Intel(R) Core(TM) i7-6700HQ 2.60GHz CPU,
we obtained the following RNG rates from molecular dice, C++ Library RNG
and GSL RNG.

|     RNG Name   |    Distribution     | Rate (doubles/sec) |      Mean     |    Samples   |
|:--------------:|:-------------------:|:------------------:|:-------------:|:------------:|
| molecular_dice | uniform             | 1.929989e+08       |  4.999870e-01 | 1.000000e+09 |
|  gsl_mt19937   | uniform             | 1.188178e+08       |  5.000021e-01 | 1.000000e+09 |
|  cpp_mt19937   | uniform             | 8.127211e+07       |  4.999967e-01 | 1.000000e+09 |

|     RNG Name   |    Distribution     | Rate (doubles/sec) |      Mean     |    Samples   |
|:--------------:|:-------------------:|:------------------:|:-------------:|:------------:|
| molecular_dice | normal              | 2.524733e+08       | -1.640548e-05 | 1.000000e+09 |
|  gsl_mt19937   | normal (ziggurat)   | 7.950162e+07       | -1.725520e-05 | 1.000000e+09 |
|  cpp_mt19937   | normal              | 2.341855e+07       | -3.711367e-05 | 1.000000e+09 |
|  gsl_mt19937   | normal (box-muller) | 1.628116e+07       | -1.831839e-05 | 1.000000e+09 |

|     RNG Name   |    Distribution     | Rate (doubles/sec) |      Mean     |    Samples   |
|:--------------:|:-------------------:|:------------------:|:-------------:|:------------:|
| molecular_dice | exponential         | 2.424187e+08       |  1.000023e+00 | 1.000000e+09 |
|  gsl_mt19937   | exponential         | 3.696452e+07       |  9.999986e-01 | 1.000000e+09 |
|  cpp_mt19937   | exponential         | 2.420209e+07       |  1.000024e+00 | 1.000000e+09 |
