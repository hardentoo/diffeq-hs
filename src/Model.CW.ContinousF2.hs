import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import qualified Graphics.Rendering.Chart.Backend.Cairo as BC
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, plot, line)
import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, fromList, toColumns)
import Numeric.GSL.ODE
import Numeric.AD

-- constant parameters
aMin, aMax, bMin, bMax, s0, s1, s2, k0, k1, g0, g1, g2, g3, g4, c0, c1, d0, d1, h, e :: Double
aMin = 0.1
aMax = 1.0
bMin = 0.1
bMax = 1.0
s0 = 7.0
s1 = 2.0
s2 = 1.0
k0 = 1.0
k1 = 0.01
g0 = 0.8
g1 = 0.8
g2 = 2.2
g3 = -0.3
g4 = 0.3
c0 = 1.0
c1 = 0.01
d0 = 0.486
d1 = 1.052
h  = 1.0
e  = 1.0

-- initial conditions
initVals :: Vector Double
initVals = fromList [ 0.5, 0.5, 0.5, 0.5 ]

-- the differential equations themselves
-- TODO: some equations can be factored out and
-- TODO: differentiation related boilerplate reduced

dx :: Double -> Double -> Double -> Double -> Double
dx x y a b = x*sigma*(1 - x/kapa) - (gamma*x*y / (1 + h*x))
  where
    kapa  = k0 + k1*a
    sigma = s0 + s1*a + s2*a**2
    gamma = g0 + g1*a + g2*b + g3*a*b + g4*b**2

dy :: Double -> Double -> Double -> Double -> Double
dy x y a b = c * (gamma*x*y / (1 + h*x)) - y*delta
  where
    c = c0 + c1*b
    gamma = g0 + g1*a + g2*b + g3*a*b + g4*b**2
    delta = d0 + d1*b

da :: Double -> Double -> Double -> Double -> Double
da x y a b = (1/e) * aa * (dsigma_a * (1 - x/k) + dgamma_a * (y / (1 + h*x)))
  where
    aa = (a - aMin)*(aMax - a)
    k  = k0 + k1*a
    -- derivatives: dζ/dα and ∂γ/∂α
    dsigma_a = diff (\ a -> auto s0 + auto s1*a + auto s2*a**2) a
    dgamma_a = head $ grad (\ [a, b] -> auto g0 + auto g1*a + auto g2*b + auto g3*a*b + auto  g4*b**2) [a, b]

db :: Double -> Double -> Double -> Double -> Double
db x y a b = (1/e) * bb * (c * dgamma_b * (x*y/(1 + h*x)) - ddelta_b)
  where
    bb = (b - bMin)*(bMax - b)
    c  = c0 + c1*b
    -- derivatives: ∂y/∂β and dδ/dβ
    dgamma_b = last $ grad (\ [a, b] -> auto g0 + auto g1*a + auto g2*b + auto g3*a*b + auto g4*b**2) [a, b]
    ddelta_b = diff (\ b -> auto d0 + auto d1*b) b

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
  return (formatTime defaultTimeLocale "%F_%H%M%S" now)

timePlot = do
  layout_title .= "Cortez / Weitz – time series"
  plot $ line "prey"     [makePlottableTuples times solution !! 0]
  plot $ line "predator" [makePlottableTuples times solution !! 1]
  plot $ line "trait α"  [makePlottableTuples times solution !! 2]
  plot $ line "trait β"  [makePlottableTuples times solution !! 3]

phasePlot = do
  layout_title .= "Cortez / Weitz – phase space"
  plot $ line "prey - predator" [ map (\ [x, y, _, _] -> (x, y)) $ toLists solution ]

writePlot filePath plot = BC.toFile def {BC._fo_format=BC.PDF} filePath plot

main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plots/CW_timePlot_" ++ timeStr ++ ".pdf") timePlot
  writePlot ("plots/CW_phasePlot_" ++ timeStr ++ ".pdf") phasePlot
