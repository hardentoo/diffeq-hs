{-# LANGUAGE RecordWildCards #-}

module MougiIwasa (runMougiIwasa) where

import Foundation
import Data.List ((!!))
import Numeric  ((**), exp)

import Graphics.Rendering.Chart.Easy ((.=), layout_title, plot, line)

import qualified Numeric.LinearAlgebra as LA
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)

import Util (diff, getNowTimeString, makePlottableTuples, writePlot)
import Parameters (Par(..), par)

-- initial conditions: [ x0, y0, u0, v0 ]
initVals :: LA.Vector Double
initVals = LA.fromList [ 0.5, 0.5, 0.1, 0.1 ]

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

eqSystem :: Double -> LA.Vector Double -> LA.Vector Double
eqSystem t vars = LA.fromList [ dx par x y u v, dy par x y u v
                              , du par x y u v, dv par x y u v ]
  where
    [x, y, u, v] = LA.toList vars

-- the time steps for which the result is given
times :: LA.Vector Double
times = LA.linspace 10000 ( 0, 9999 :: Double )

-- the solutions matrix
solution :: LA.Matrix Double
solution = odeSolveV
  RKf45    -- ODE Method
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
  plot $ line "prey - predator" [ fmap (\ [x, y, _, _] -> (x, y)) $ LA.toLists solution ]


-- to make a bifurcation we solve the system of ODEs for a range of parameters.
-- usually one parameter is varied in a specified range and with a step size
-- that allows the computation to terminate in a reasonable time.

-- the approach is to create a vector of values for the parameter we want
-- to create the bifurcation diagram for.
-- mapping 
bifurcate = undefined


runMougiIwasa :: IO ()
runMougiIwasa = do
  timeStr <- getNowTimeString
  writePlot ("plots/MI_timePlot_" <> timeStr <> ".pdf") timePlot
  writePlot ("plots/MI_phasePlot_" <> timeStr <> ".pdf") phasePlot
