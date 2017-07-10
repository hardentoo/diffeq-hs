module Model.MougiIwasa where

import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import qualified Graphics.Rendering.Chart.Backend.Cairo as BC
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, plot, line)
import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, fromList, toColumns)
import Numeric.GSL.ODE (ODEMethod(..), odeSolveV)
import Numeric.GSL.Differentiation (derivCentral)

-- alias differentiation function to central derivative with initial step size 0.01
diff :: (Double -> Double) -> Double -> Double
diff fun point = fst $ derivCentral 0.01 fun point

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
initVals = fromList [ 0.5, 0.5, 0.5, 0.5 ]

-- sub functions
a :: Double -> Double -> Double
a u v = a0 / (1 + exp (th*(u-v)))

r, g :: Double -> Double
r u = r0 * (1 - u**rhX)
g v = g0 * (1 - v**rhY)

-- differential equations
dx, dy, du, dv, wx, wy :: Double -> Double -> Double -> Double -> Double
dx x y u v = (r u - e*x - a u v *y/(1 + a u v *h*x))*x
dy x y u v = (g v *       a u v *x/(1 + a u v *h*x) - d)*y
du x y u v = gx * diff (flip (wx x y) v) u -- dWx_du
dv x y u v = gy * diff       (wy x y u)  v -- dWy_dv

-- fitness defined as per capita growth rate of x and y
wx x y u v = 1/x * dx x y u v
wy x y u v = 1/y * dy x y u v

eqSystem :: Double -> Vector Double -> Vector Double
eqSystem t vars = fromList [ dx x y u v, dy x y u v
                           , du x y u v, dv x y u v ]
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
  layout_title .= "Mougi / Iwasa – time series"
  plot $ line "prey"     [makePlottableTuples times solution !! 0]
  plot $ line "predator" [makePlottableTuples times solution !! 1]
  plot $ line "trait u"  [makePlottableTuples times solution !! 2]
  plot $ line "trait v"  [makePlottableTuples times solution !! 3]

phasePlot = do
  layout_title .= "Mougi / Iwasa – phase space"
  plot $ line "prey - predator" [ map (\ [x, y, _, _] -> (x, y)) $ toLists solution ]

writePlot filePath plot = BC.toFile def {BC._fo_format=BC.PDF} filePath plot

main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plots/MI_timePlot_" ++ timeStr ++ ".pdf") timePlot
  writePlot ("plots/MI_phasePlot_" ++ timeStr ++ ".pdf") phasePlot
