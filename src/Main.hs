import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import Graphics.Rendering.Chart.Backend.Diagrams(toFile)
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, setColors, opaque, blue, red, plot, line)

import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toColumns)
import Numeric.GSL.ODE

type D = Double

-- TODO: bifurcation diagram
-- TODO: write to CSV

-- TODO: nicer way to manage parameters, instead of passing globals?
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
dx x y = α*x   - β*x*y -- prey: dx = ax - bxy
dy x y = δ*x*y - γ*y   -- pred: dy = gxy - dy

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

writePlot :: FilePath -> IO ()
writePlot filePath = toFile def filePath $ do
  layout_title .= "Lotka-Volterra"
  setColors [opaque blue, opaque red]
  plot $ line "prey"     [makePlottableTuples time sol !! 0]
  plot $ line "predator" [makePlottableTuples time sol !! 1]

main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plot_" ++ timeStr ++ ".svg")
