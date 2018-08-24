// Copyright (c) 2018 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

namespace md {

// matrix for rotations in 3D space
struct rotation_matrix
{
  double xx;
  double xy;
  double xz;

  double yx;
  double yy;
  double yz;

  double zx;
  double zy;
  double zz;
};

// result of rotating a velocity vector using the
// rotation matrix
velocity
operator*(const rotation_matrix& R, const velocity& u)
{
  velocity v;
  v.vx = R.xx * u.vx + R.xy * u.vy + R.xz * u.vz;
  v.vy = R.yx * u.vx + R.yy * u.vy + R.yz * u.vz;
  v.vz = R.zx * u.vx + R.zy * u.vy + R.zz * u.vz;
  return v;
}

} // namespace md
