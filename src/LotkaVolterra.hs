module LotkaVolterra where

import Data.List ((!!))
import Foundation

import Graphics.Rendering.Chart.Easy ((.=), layout_title, setColors, opaque, blue, red, plot, line)

import qualified Numeric.LinearAlgebra as LA
import Numeric.GSL.ODE

import Util

type D = Double

x0, y0, α, β, γ, δ :: D
-- initial conditions --
x0 = 0.5 -- initial value of x
y0 = 0.5 -- initial value of y
-- parameters --
α  = 0.2 -- reproduction rate
β  = 0.5 -- attack rate
γ  = 0.1 -- mortality rate
δ  = 0.1 -- conversion efficiency

-- the differential equations themselves
dx, dy :: D -> D -> D
dx x y = α*x   - β*x*y -- prey
dy x y = δ*x*y - γ*y   -- predator

-- the differential equations as the system to solve
eqs :: D -> [D] -> [D]
eqs t [x, y] = [ dx x y, dy x y ]

-- the time (steps)
time :: LA.Vector D
time = LA.linspace 1000 (0, 999 :: D)

-- the solutions matrix
sol :: LA.Matrix D
sol = odeSolve eqs [x0, y0] time

timePlot = do
  layout_title .= "Lotka-Volterra – time series"
  setColors [opaque blue, opaque red]
  plot $ line "prey"     [makePlottableTuples time sol !! 0]
  plot $ line "predator" [makePlottableTuples time sol !! 1]

phasePlot = do
  layout_title .= "Lotka-Volterra – phase space"
  setColors [opaque blue]
  plot $ line "prey - predator" [ fmap (\[x, y] -> (x, y)) $ LA.toLists sol]

main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plots/LV_timePlot_" <> timeStr <> ".pdf") timePlot
  writePlot ("plots/LV_phasePlot_" <> timeStr <> ".pdf") phasePlot
