// Copyright (c) 2018 Santosh Ansumali @ JNCASR
// See LICENSE

#pragma once

namespace md {

// velocity components of a particle
class velocity
{
public:
  velocity&
  operator+=(const velocity& rhs)
  {
    this->vx += rhs.vx;
    this->vy += rhs.vy;
    this->vz += rhs.vz;
    return *this;
  }

  velocity&
  operator-=(const velocity& rhs)
  {
    this->vx -= rhs.vx;
    this->vy -= rhs.vy;
    this->vz -= rhs.vz;
    return *this;
  }

  velocity&
  operator*=(const double rhs)
  {
    this->vx *= rhs;
    this->vy *= rhs;
    this->vz *= rhs;
    return *this;
  }

  velocity&
  operator/=(const double rhs)
  {
    this->vx /= rhs;
    this->vy /= rhs;
    this->vz /= rhs;
    return *this;
  }

  double vx;
  double vy;
  double vz;
};

velocity
operator+(const velocity& a, const velocity& b)
{
  velocity c;
  c.vx = a.vx + b.vx;
  c.vy = a.vy + b.vy;
  c.vz = a.vz + b.vz;
  return c;
}

velocity
operator-(const velocity& a, const velocity& b)
{
  velocity c;
  c.vx = a.vx - b.vx;
  c.vy = a.vy - b.vy;
  c.vz = a.vz - b.vz;
  return c;
}

double
operator*(const velocity& a, const velocity& b)
{
  return a.vx * b.vx + a.vy * b.vy + a.vz * b.vz;
}

velocity
operator*(const velocity& a, const double b)
{
  velocity c;
  c.vx = a.vx * b;
  c.vy = a.vy * b;
  c.vz = a.vz * b;
  return c;
}

velocity
operator*(const double b, const velocity& a)
{
  return a * b;
}

velocity
operator/(const velocity& a, const double b)
{
  velocity c;
  c.vx = a.vx / b;
  c.vy = a.vy / b;
  c.vz = a.vz / b;
  return c;
}

} // namespace md
