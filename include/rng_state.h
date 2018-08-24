// Copyright (c) 2018 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

#include <cstddef>
#include <vector>
#include "position.h"
#include "velocity.h"
#include "rotation_matrix.h"

namespace md {

class rng_state
{
public:
  // dimension of the particle system
  static const std::size_t dim = 3;

  // constructor
  rng_state()
  { initialize(0); }

  rng_state(const std::size_t num)
  { initialize(num); }

  std::size_t
  num_particles() const
  { return m_vel.size(); }

  // accessors for velocity and position of each particle
  position&
  pos(const std::size_t idx)
  { return m_pos[idx]; }

  const position&
  pos(const std::size_t idx) const
  { return m_pos[idx]; }

  velocity&
  vel(const std::size_t idx)
  { return m_vel[idx]; }

  const velocity&
  vel(const std::size_t idx) const
  { return m_vel[idx]; }

  void
  initialize(const std::size_t num)
  {
    m_pos.resize(num);
    m_vel.resize(num);
  }

  // wrap coordinates to lie within [0,1]
  double
  periodic_wrap(const double x) const
  { return x - 1.0 * (x > 1.0) + 1.0 * (x < 0.0); }

  // update position of a particle
  void
  update_pos(const std::size_t idx, const double dt)
  {
    m_pos[idx].x = periodic_wrap(m_pos[idx].x + m_vel[idx].vx * dt);
    m_pos[idx].y = periodic_wrap(m_pos[idx].y + m_vel[idx].vy * dt);
    m_pos[idx].z = periodic_wrap(m_pos[idx].z + m_vel[idx].vz * dt);
    return;
  }

  // update positions of all particles
  void
  update_all_pos(const double dt)
  {
    for (std::size_t i = 0; i < num_particles(); i++) {
      update_pos(i, dt);
    }
    return;
  }

  // update velocities of a pair of colliding particles
  void
  update_vel(const rotation_matrix& R,
             const std::size_t      idx_a,
             const std::size_t      idx_b)
  {
    const velocity ua   = m_vel[idx_a];
    const velocity ub   = m_vel[idx_b];
    const velocity urel = ua - ub;
    const velocity vrel = R * urel;
    const velocity ucm  = 0.5 * (ua + ub);
    m_vel[idx_a]        = ucm + vrel;
    m_vel[idx_b]        = ucm - vrel;
    return;
  }

  // update state by colliding two particles
  void
  update(const rotation_matrix& R,
         const std::size_t      idx_a,
         const std::size_t      idx_b,
         const bool             update_positions,
         const double           dt)
  {
    update_vel(R, idx_a, idx_b);
    if (update_positions) {
      update_pos(idx_a, dt);
      update_pos(idx_b, dt);
    }
    return;
  }

private:
  // position and velocity of each particle in the system
  std::vector<position> m_pos;
  std::vector<velocity> m_vel;
};

} // namespace md
