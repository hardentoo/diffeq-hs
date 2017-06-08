import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Diagrams(toFile)

-- import qualified Data.Vector as V
import qualified Data.Vector.Storable as V
import Numeric.GSL.ODE
import qualified Numeric.LinearAlgebra as LA

xdot :: Double -> [Double] -> [Double]
xdot t [x, y] =
  [ 0.2*x - 0.5*x*y -- prey: dx = ax - bxy
  , 0.5*x*y - 0.1*y -- pred: dy = gxy - dy
  ]

time :: LA.Vector Double
time = LA.linspace 1000 (0, 1000 :: Double)

-- only the solutions
sol :: LA.Matrix Double
sol = odeSolve xdot [0.5, 0.5] time

-- takes
makePlottableTuples :: LA.Vector Double -> LA.Matrix Double -> [[(Double, Double)]]
makePlottableTuples v m = map (zip (LA.toList v) . LA.toList) (LA.toColumns m)

writePlot :: String -> IO ()
writePlot filePath = toFile def filePath $ do
  -- create plot here
  layout_title .= "Lotka-Volterra"
  setColors [opaque blue, opaque red]
  plot $ line "prey"     [makePlottableTuples time sol !! 0]
  plot $ line "predator" [makePlottableTuples time sol !! 1]

getNowTimeString :: IO String
getNowTimeString = do
  now <- getCurrentTime
  return (formatTime defaultTimeLocale "%F_%H%M%S" now)

-- main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plot_" ++ timeStr ++ ".svg")

