import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import qualified Graphics.Rendering.Chart.Backend.Cairo as BC
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, setColors, opaque, blue, red, plot, line)

import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, toColumns)
import Numeric.GSL.ODE

type D = Double

-- TODO: bifurcation diagram
-- TODO: bifurcation diagram in parallel
-- TODO: write bifurcation data to CSV
-- TODO: more concise ways to do different plots etc... map write function?
-- TODO: find phase relationship (Fourier? min/max correlation)
-- TODO: nicer way to manage parameters, instead of passing globals?
-- TODO: factor out models

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
time :: Vector D
time = linspace 1000 (0, 999 :: D)

-- the solutions matrix
sol :: Matrix D
sol = odeSolve eqs [x0, y0] time

-- takes the time vector and the solution matrix
-- creates a list of lists containing tuples
-- for each variable and the corresponding time
-- necessary to satisfy plotting funtion
makePlottableTuples :: Vector D -> Matrix D -> [[(D, D)]]
makePlottableTuples v m = map (zip (toList v) . toList) (toColumns m)

-- create a string like "2017-06-09_131211"
-- as a timestamp to use for writing files
getNowTimeString :: IO String
getNowTimeString = do
  now <- getCurrentTime
  return (formatTime defaultTimeLocale "%F_%H%M%S" now)

timePlot = do
  layout_title .= "Lotka-Volterra – time series"
  setColors [opaque blue, opaque red]
  plot $ line "prey"     [makePlottableTuples time sol !! 0]
  plot $ line "predator" [makePlottableTuples time sol !! 1]

phasePlot = do
  layout_title .= "Lotka-Volterra – phase space"
  setColors [opaque blue]
  plot $ line "prey - predator" [ map (\[x, y] -> (x, y)) $ toLists sol]

writePlot filePath plot = BC.toFile def {BC._fo_format=BC.PDF} filePath plot

main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plots/LV_timePlot_" ++ timeStr ++ ".pdf") timePlot
  writePlot ("plots/LV_phasePlot_" ++ timeStr ++ ".pdf") phasePlot
