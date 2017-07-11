{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveGeneric #-}

module MougiIwasa (runMougiIwasa) where

import Graphics.Rendering.Chart.Easy ((.=), layout_title, plot, line)

import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, fromList)
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)

import Util (diff, getNowTimeString, makePlottableTuples, writePlot)

import qualified Dhall as D

data Par = Par {
  e, h, d, gx, gy, th, rhX, rhY, a0, r0, g0 :: Double
} deriving (Show, D.Generic)

instance D.Interpret Par

-- constant parameters
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

-- initial conditions
initVals :: Vector Double
initVals = fromList [ 0.5, 0.5, 0.5, 0.5 ]

-- sub functions
a :: Par -> Double -> Double -> Double
a Par{..} u v = a0 / (1 + exp (th*(u-v)))

r, g :: Par -> Double -> Double
r Par{..} u = r0 * (1 - u**rhX)
g Par{..} v = g0 * (1 - v**rhY)

-- differential equations
dx, dy, du, dv, wx, wy :: Par -> Double -> Double -> Double -> Double -> Double
dx Par{..} x y u v = (r par u - e*x - a par u v *y/(1 + a par u v *h*x))*x
dy Par{..} x y u v = (g par v *       a par u v *x/(1 + a par u v *h*x) - d)*y
du Par{..} x y u v = gx * diff (flip (wx par x y) v) u -- dWx_du
dv Par{..} x y u v = gy * diff       (wy par x y u)  v -- dWy_dv

-- fitness defined as per capita growth rate of x and y
wx Par{..} x y u v = 1/x * dx par x y u v
wy Par{..} x y u v = 1/y * dy par x y u v

eqSystem :: Double -> Vector Double -> Vector Double
eqSystem t vars = fromList [ dx par x y u v, dy par x y u v
                           , du par x y u v, dv par x y u v ]
  where
    [x, y, u, v] = toList vars

-- the time steps for which the result is given
times :: Vector Double
times = linspace 2000 (0, 999 :: Double)

-- the solutions matrix
solution :: Matrix Double
solution = odeSolveV
  RK8pd    -- ODE Method
  1E-8     -- initial step size
  1E-8     -- absolute tolerance for the state vector
  0        -- relative tolerance for the state vector
  eqSystem -- differential eqations: xdot(t,x), ...
  initVals -- inital conditions [ x0, y0, α0, β0 ]
  times    -- desired solution times

timePlot = do
  layout_title .= "Mougi / Iwasa – time series"
  plot $ line "prey"     [makePlottableTuples times solution !! 0]
  plot $ line "predator" [makePlottableTuples times solution !! 1]
  plot $ line "trait u"  [makePlottableTuples times solution !! 2]
  plot $ line "trait v"  [makePlottableTuples times solution !! 3]

phasePlot = do
  layout_title .= "Mougi / Iwasa – phase space"
  plot $ line "prey - predator" [ map (\ [x, y, _, _] -> (x, y)) $ toLists solution ]

runMougiIwasa :: IO ()
runMougiIwasa = do
  timeStr <- getNowTimeString
  writePlot ("plots/MI_timePlot_" ++ timeStr ++ ".pdf") timePlot
  writePlot ("plots/MI_phasePlot_" ++ timeStr ++ ".pdf") phasePlot
