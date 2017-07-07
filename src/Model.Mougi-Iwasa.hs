import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import Graphics.Rendering.Chart.Backend.Cairo (toFile, FileFormat(..), _fo_format)
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, plot, line)
import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, fromList, toColumns)
import Numeric.GSL.ODE (odeSolveV, ODEMethod(..))
import Numeric.AD

-- data Par = Par {} deriving (Show, Generic)

-- instance DC.Interpret Par

-- constant parameters

-- initial conditions
initVals :: Vector Double
initVals = fromList [ 0.5, 0.5, 0.5, 0.5 ]

-- the differential equations themselves
-- TODO: some equations can be factored out and
-- TODO: differentiation related boilerplate reduced

dx :: Double -> Double -> Double -> Double -> Double
dx x y a b = undefined

dy :: Double -> Double -> Double -> Double -> Double
dy x y a b = undefined

da :: Double -> Double -> Double -> Double -> Double
da x y a b = undefined

db :: Double -> Double -> Double -> Double -> Double
db x y a b = undefined

-- HOW CAN I MAKE ALL FUNCTIONS NOT USE GLOBALS?
eqSystem :: Double -> Vector Double -> Vector Double
eqSystem t vars = fromList [ dx x y a b
                           , dy x y a b
                           , da x y a b
                           , db x y a b ]
  where
    [x, y, a, b] = toList vars

-- the time (steps)
times :: Vector Double
times = linspace 5000 (0, 499 :: Double)

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

-- takes the time vector and the solution matrix
-- creates a list of lists containing tuples
-- for each variable and the corresponding time
-- necessary to satisfy plotting funtion
makePlottableTuples :: Vector Double -> Matrix Double -> [[(Double, Double)]]
makePlottableTuples v m = map (zip (toList v) . toList) (toColumns m)

-- create a string like "2017-06-09_131211"
-- as a timestamp to use for writing files
getNowTimeString :: IO String
getNowTimeString = do
  now <- getCurrentTime
  pure (formatTime defaultTimeLocale "%F_%H%M%S" now)

timePlot = do
  layout_title .= "Cortez / Weitz – time series"
  plot $ line "prey"     [makePlottableTuples times solution !! 0]
  plot $ line "predator" [makePlottableTuples times solution !! 1]
  plot $ line "trait α"  [makePlottableTuples times solution !! 2]
  plot $ line "trait β"  [makePlottableTuples times solution !! 3]

phasePlot = do
  layout_title .= "Cortez / Weitz – phase space"
  plot $ line "prey - predator" [ map (\ [x, y, _, _] -> (x, y)) $ toLists solution ]

writePlot filePath plot = toFile def {_fo_format=PDF} filePath plot

main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plots/CW_timePlot_" ++ timeStr ++ ".pdf") timePlot
  writePlot ("plots/CW_phasePlot_" ++ timeStr ++ ".pdf") phasePlot
