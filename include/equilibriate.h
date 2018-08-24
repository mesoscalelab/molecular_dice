// Copyright (c) 2018 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

#include <cstddef>
#include <random>
#include "velocity.h"
#include "rng_state.h"

namespace md {

// initialize positions
template <typename XR>
void
equilibriate_positions(rng_state& s, XR& xr)
{
  // generate uniform distribution of positions
  std::uniform_real_distribution<double> uniform(0., 1.);
  for (std::size_t i = 0; i < s.num_particles(); i++) {
    s.pos(i).x = uniform(xr);
    s.pos(i).y = uniform(xr);
    s.pos(i).z = uniform(xr);
  }
  return;
}

// initialize velocities
template <typename XR>
void
equilibriate_velocities(rng_state& s, XR& xr)
{
  const double temperature = 2.;
  const double stddev = std::sqrt(temperature);
  
  // generate normal distribution of velocities
  std::normal_distribution<double> normal(0., stddev);
  for (std::size_t i = 0; i < s.num_particles(); i++) {
    s.vel(i).vx = normal(xr);
    s.vel(i).vy = normal(xr);
    s.vel(i).vz = normal(xr);
  }
  // force center of mass velocity to zero
  velocity avg_vel_cm;
  avg_vel_cm.vx = 0;
  avg_vel_cm.vy = 0;
  avg_vel_cm.vz = 0;
  for (std::size_t i = 0; i < s.num_particles(); i++) {
    avg_vel_cm += s.vel(i);
  }
  avg_vel_cm /= (1.0 * s.num_particles());
  for (std::size_t i = 0; i < s.num_particles(); i++) {
    s.vel(i) -= avg_vel_cm;
  }

  // enforce chosen temperature value
  double avg_energy = 0;
  for (std::size_t i = 0; i < s.num_particles(); i++) {
    avg_energy += s.vel(i) * s.vel(i);
  }
  avg_energy /= (1. * rng_state::dim * s.num_particles());
  for (std::size_t i = 0; i < s.num_particles(); i++) {
    s.vel(i) = (stddev * s.vel(i)) / std::sqrt(avg_energy);
  }
  return;
}

} // namespace md
