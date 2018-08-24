// Copyright (c) 2018 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

#include <iostream>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <random>
#include "equilibriate.h"
#include "rng.h"

namespace md {

// constructor
rng::rng(unsigned long     seed,
         const std::size_t num,
         const double      dt)
: m_max_unip_buffers_filled(calc_max_unip_buffers_filled(num)),
  m_max_pairs_collided(calc_max_pairs_collided(num)),
  m_dt(dt)
{
  // check validity of arguments
  if (m_max_pairs_collided < 2) {
    throw std::invalid_argument("use more particles for RNG state");
  }

  // initialize particle system
  m_state.initialize(num);

  // instantly equilibriate the particle system by
  // initializing positions and velocities of particles with
  // equilibrium distribution values using an external RNG
  std::mt19937 xr(seed);
  equilibriate_positions(m_state, xr);
  equilibriate_velocities(m_state, xr);
  m_state.update_all_pos(m_dt);
 
  // fill internal uniform RNG buffer and use
  // it to initialize randomized parameters
  refresh_rand_rot_matrix_params();
  refresh_rand_pair_select_params();
}

// random number generation calls
// ------------------------------
// in each case, the RNG call refills it's respective buffer with new values
// if all values present in the buffer have been used up; then serves the
// first unused value in the buffer

double
rng::uniform()
{
  if (m_num_unifs_used == 0 or m_num_unifs_used >= m_unif_buffer.size()) {
    refill_unif_buffer();
    m_num_unifs_used = 0;
  }
  return m_unif_buffer[m_num_unifs_used++];
}

double
rng::normal()
{
  if (m_num_norms_used == 0 or m_num_norms_used >= m_norm_buffer.size()) {
    refill_norm_buffer();
    m_num_norms_used = 0;
  }
  return m_norm_buffer[m_num_norms_used++];
}

double
rng::exp()
{
  if (m_num_expos_used == 0 or m_num_expos_used >= m_expo_buffer.size()) {
    refill_expo_buffer();
    m_num_expos_used = 0;
  }
  return m_expo_buffer[m_num_expos_used++];
}

double
rng::uniform_private()
{
  if (m_num_unips_used == 0 or m_num_unips_used >= m_unip_buffer.size()) {
    refill_unip_buffer();
    m_num_unips_used = 0;
  }
  return m_unip_buffer[m_num_unips_used++];
}

// calculate values of constant parameters
// ---------------------------------------
std::size_t
rng::calc_max_unip_buffers_filled(const std::size_t num) const
{
  return (dim * num) / m_unip_buffer.size();
}

std::size_t
rng::calc_max_pairs_collided(const std::size_t num) const
{
  return num / 8;
}

// assignment of randomized parameters
// -----------------------------------

// update positions of all particles in the particle system if
// all position coordinates have already been used as random
// numbers, so that the new updated positions can be used as a
// source of uniform real RNGs for the private uniform RNG
void
rng::refresh_unip_pool()
{
  if (m_num_unip_buffers_filled >= m_max_unip_buffers_filled) {
    m_state.update_all_pos(m_dt);
    m_num_unip_buffers_filled = 0;
  }
  return;
}

// fill elements of the 3D rotation matrix with
// values corresponding to a random triplet of
// Eulerian angles
void
rng::refresh_rand_rot_matrix_params()
{
  const double alpha    = 2. * M_PI * uniform_private();
  const double theta    = 1. * M_PI * uniform_private();
  const double phi      = 2. * M_PI * uniform_private();
  const double nx       = std::sin(theta) * std::cos(phi);
  const double ny       = std::sin(theta) * std::sin(phi);
  const double nz       = std::cos(theta);
  const double defl_cos = std::cos(alpha); 
  const double defl_sin = std::sin(alpha); 
  m_rot_matrix.xx       = 0.5 * (nx * nx * (1. - defl_cos) + 1. * defl_cos);
  m_rot_matrix.xy       = 0.5 * (nx * ny * (1. - defl_cos) - nz * defl_sin);
  m_rot_matrix.xz       = 0.5 * (nx * nz * (1. - defl_cos) + ny * defl_sin);
  m_rot_matrix.yx       = 0.5 * (ny * nx * (1. - defl_cos) + nz * defl_sin);
  m_rot_matrix.yy       = 0.5 * (ny * ny * (1. - defl_cos) + 1. * defl_cos);
  m_rot_matrix.yz       = 0.5 * (ny * nz * (1. - defl_cos) - nx * defl_sin);
  m_rot_matrix.zx       = 0.5 * (nz * nx * (1. - defl_cos) - ny * defl_sin);
  m_rot_matrix.zy       = 0.5 * (nz * ny * (1. - defl_cos) + nx * defl_sin);
  m_rot_matrix.zz       = 0.5 * (nz * nz * (1. - defl_cos) + 1. * defl_cos);
  return;
}

// assign randomized values to collision pair selection parameters
// according to the pair selection scheme
void
rng::refresh_rand_pair_select_params()
{
  const std::size_t num = m_state.num_particles();
  m_start = static_cast<int>(uniform_private() * num);
  m_shift = static_cast<int>(uniform_private() * (num / (m_max_pairs_collided - 1.) - 1.)) + 1;
  m_jump  = static_cast<int>(uniform_private() * (num - 1)) + 1;
  return;
}

// refill the rotation matrix and pair selection parameters
// with a new set of randomized values if the maximum threshold
// for number of pairs collided in the process of random number
// generation has been exceeded
void
rng::refresh_rand_params()
{
  if (m_num_pairs_collided >= m_max_pairs_collided) {
    refresh_rand_rot_matrix_params();
    refresh_rand_pair_select_params();
    m_num_pairs_collided = 0;
  }
  return;
}

// set indices for a new pair of particles which will be used
// for the next collision event
void
rng::refresh_collision_pair()
{
  const std::size_t num = m_state.num_particles();
  m_idx_a  = m_start + m_num_pairs_collided * m_shift;
  m_idx_a -= num * (m_idx_a >= num);
  m_idx_b  = m_idx_a + m_jump;
  m_idx_b -= num * (m_idx_b >= num);
  return;
}

// refill RNG buffers
// ------------------

// sample position coordinates of two successive particles
// as uniformly distributed random variates for internal use
void
rng::refill_unip_buffer()
{
  refresh_unip_pool();
  const std::size_t idx_a = 2 * m_num_unip_buffers_filled + 0;
  const std::size_t idx_b = 2 * m_num_unip_buffers_filled + 1;
  m_num_unip_buffers_filled++;

  m_unip_buffer[0] = m_state.pos(idx_a).x;
  m_unip_buffer[1] = m_state.pos(idx_a).y;
  m_unip_buffer[2] = m_state.pos(idx_a).z;
  m_unip_buffer[3] = m_state.pos(idx_b).x;
  m_unip_buffer[4] = m_state.pos(idx_b).y;
  m_unip_buffer[5] = m_state.pos(idx_b).z;
  return;
}

// sample position coordinates of collided particle pair
// as uniformly distributed random variates
void
rng::refill_unif_buffer()
{
  refresh_rand_params();
  refresh_collision_pair();
  m_state.update(m_rot_matrix, m_idx_a, m_idx_b, true, m_dt);
  m_num_pairs_collided++;

  m_unif_buffer[0] = m_state.pos(m_idx_a).x;
  m_unif_buffer[1] = m_state.pos(m_idx_a).y;
  m_unif_buffer[2] = m_state.pos(m_idx_a).z;
  m_unif_buffer[3] = m_state.pos(m_idx_b).x;
  m_unif_buffer[4] = m_state.pos(m_idx_b).y;
  m_unif_buffer[5] = m_state.pos(m_idx_b).z;
  return;
}

// sample components of relative outgoing velocity between
// collided particle pair as normally distributed random variates
void
rng::refill_norm_buffer()
{
  refresh_rand_params();
  refresh_collision_pair();
  m_state.update(m_rot_matrix, m_idx_a, m_idx_b, false, 0);
  m_num_pairs_collided++;
  const velocity out_vel_rel = 0.5 * (m_state.vel(m_idx_a) - m_state.vel(m_idx_b));
  m_norm_buffer[0] = out_vel_rel.vx;
  m_norm_buffer[1] = out_vel_rel.vy;
  m_norm_buffer[2] = out_vel_rel.vz;
  return;
}

// sample, along each axis, the average kinetic energy between
// collided particle pair as exponentially distributed random variates
void
rng::refill_expo_buffer()
{
  refresh_rand_params();
  refresh_collision_pair();
  m_state.update(m_rot_matrix, m_idx_a, m_idx_b, false, 0);
  m_num_pairs_collided++;

  const velocity vel_a = m_state.vel(m_idx_a);
  const velocity vel_b = m_state.vel(m_idx_b);
  m_expo_buffer[0] = 0.25 * (vel_a.vx * vel_a.vx + vel_b.vx * vel_b.vx);
  m_expo_buffer[1] = 0.25 * (vel_a.vy * vel_a.vy + vel_b.vy * vel_b.vy);
  m_expo_buffer[2] = 0.25 * (vel_a.vz * vel_a.vz + vel_b.vz * vel_b.vz);
  return;
}

} // namespace md
