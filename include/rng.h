// Copyright (c) 2018 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

#include <cstddef>
#include <array>
#include <string>
#include "rotation_matrix.h"
#include "rng_state.h"

namespace md {

class rng
{
public:
  // dimension of particle system
  static const std::size_t dim = rng_state::dim;

  // constructor
  // arguments::
  // seed   : seed passed to external RNG for initializing
  //          particles to an equilibrium state
  // num    : number of particles in the RNG state
  // dt     : time gap between successive collisions
  rng(unsigned long     seed = 1234,
      const std::size_t num  = 131072,
      const double      dt   = 0.1);

  // random number generation calls
  // ------------------------------

  // generates a random real uniformly distributed in (0,1]
  double uniform();

  // generates a random real normally distributed with
  // mean 0 and variance 1
  double normal();

  // generates a random real exponentially distributed as exp(-x)
  // where x lies in [0,inf)
  double exp();

private:
  // generates a random real uniformly distributed in (0,1]
  // note: this function is only for internal use for setting
  // random parameters private to the rng class
  double uniform_private();

  // calculate values of fixed parameters
  std::size_t
  calc_max_unip_buffers_filled(const std::size_t num) const;

  std::size_t
  calc_max_pairs_collided(const std::size_t num) const;

  // assignment of randomized parameters
  void refresh_unip_pool();
  void refresh_rand_rot_matrix_params();
  void refresh_rand_pair_select_params();
  void refresh_rand_params();
  void refresh_collision_pair();

  // refill RNG buffers
  void refill_unip_buffer();
  void refill_unif_buffer();
  void refill_norm_buffer();
  void refill_expo_buffer();

  // particle system which acts as the RNG state
  rng_state m_state;

  // buffers for storing multiple random numbers
  // generated during one collision process
  std::array<double, 2 * dim> m_unip_buffer;
  std::array<double, 1 * dim> m_norm_buffer;
  std::array<double, 2 * dim> m_unif_buffer;
  std::array<double, 1 * dim> m_expo_buffer;

  // counts of random numbers used from each buffer
  std::size_t m_num_unips_used = 0;
  std::size_t m_num_norms_used = 0;
  std::size_t m_num_unifs_used = 0;
  std::size_t m_num_expos_used = 0;

  // count of buffers used up during internal
  // uniform RNG process, determines when the
  // internal RNG pool should be refreshed
  const std::size_t m_max_unip_buffers_filled;
  std::size_t m_num_unip_buffers_filled = 0;

  // 3D rotation matrix for pair collision 
  rotation_matrix m_rot_matrix;

  // count of pairs collided during RNG process,
  // determines when a new set of randomized
  // parameters must be brought in
  const std::size_t m_max_pairs_collided;
  std::size_t m_num_pairs_collided = 0;

  // collision pair selection scheme and parameters
  std::size_t m_start = 0;
  std::size_t m_shift = 0;
  std::size_t m_jump  = 0;

  // indices of particle pair used for most recent collision
  std::size_t m_idx_a;
  std::size_t m_idx_b;

  // time gap between consecutive collisions
  const double m_dt;
};

} // namespace md
