import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import Graphics.Rendering.Chart.Backend.Cairo (toFile, FileFormat(..), _fo_format)
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, plot, line)
import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, fromList, toColumns)
import Numeric.GSL.ODE (odeSolveV, ODEMethod(..))
import Numeric.AD hiding (du)

-- constant parameters
e   = 1.0
h   = 0.1
d   = 0.1
gx  = 0.01
gy  = 0.01
th  = 10.0
rhX = 2.0
rhY = 2.0
a0  = 1.0
r0  = 1.0
g0  = 1.0

-- initial conditions
initVals :: Vector Double
initVals = fromList [ 0.1, 0.1, 0.1, 0.1 ]

-- the differential equations themselves
-- TODO: some equations can be factored out and
-- TODO: differentiation related boilerplate reduced

dx :: Double -> Double -> Double -> Double -> Double
dx x y u v = (r u v - e*x - a u v *y/(1 + a u v *h*x))*x

dy :: Double -> Double -> Double -> Double -> Double
dy x y u v = (g* a u v *y/(1 + a u v *h*x) - d)*y

du :: Double -> Double -> Double -> Double -> Double
du x y u v = gx * dWx_du
  where
    dWx_du x y u v = grad (\ [x y u v] -> auto (r u v) - auto e * x - auto (a u v) )

-- dgamma_b = last $ grad (\ [a, b] -> auto g0 + auto g1*a + auto g2*b + auto g3*a*b + auto g4*b**2) [a, b]

dv :: Double -> Double -> Double -> Double -> Double
dv x y u v = gy * dWy_dv
  where
    dWy_dv = undefined

-- sub functions
a u v = a0 / (1 + exp (th*(u-v)))
r u v = r0 * (1 - u^rhX)
g u v = g0 * (1 - v^rhY)

wx :: Double -> Double -> Double -> Double -> Double
wx x y u v = r u v - e*x - a u v *y/(1 + a u v *h*x)

wy :: Double -> Double -> Double -> Double -> Double
wy x y u v = g* a u v *y/(1 + a u v *h*x) - d

-- HOW CAN I MAKE ALL FUNCTIONS NOT USE GLOBALS?
eqSystem :: Double -> Vector Double -> Vector Double
eqSystem t vars = fromList [ dx x y u v
                           , dy x y u v
                           , du x y u v
                           , dv x y u v ]
  where
    [x, y, u, v] = toList vars

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
