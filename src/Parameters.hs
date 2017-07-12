module Parameters where

data Par = Par { e, h, d, gx, gy, th, rhX, rhY, a0, r0, g0 :: Double }

par :: Par
par = Par {
  e   = 1.0,
  h   = 0.1,
  d   = 0.1,
  gx  = 0.01,
  gy  = 0.01,
  th  = 10.0,
  rhX = 2.0,
  rhY = 2.0,
  a0  = 1.0,
  r0  = 1.0,
  g0  = 1.0
}

