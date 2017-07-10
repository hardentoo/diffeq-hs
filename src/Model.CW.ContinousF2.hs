import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import qualified Graphics.Rendering.Chart.Backend.Cairo as BC
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, plot, line)
import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, fromList, toColumns)
import Numeric.GSL.ODE
import Numeric.GSL.Differentiation
import Numeric.AD

-- alias differentiation function to central derivative with initial step size 0.01
diff_ fun point = fst $ derivCentral 0.01 fun point

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

-- sub functions
g  a b = g0 + g1*a + g2*b + g3*a*b + g4*b**2
s  a   = s0 + s1*a + s2*a**2
k  a   = k0 + k1*a
c  b   = c0 + c1*b
d  b   = d0 + d1*b
aa a   = (a - aMin)*(aMax - a)
bb b   = (b - bMin)*(bMax - b)

-- the differential equations themselves
dx, dy, da, db :: Double -> Double -> Double -> Double -> Double

dx x y a b = x* s a *(1 - x/ k a) - (g a b *x*y / (1 + h*x))

dy x y a b = c b * (g a b *x*y / (1 + h*x)) - y* d b

da x y a b = (1/e) * aa a * (ds_a * (1 - x/ k a) + dg_a * (y / (1 + h*x)))
  where
    -- derivatives: dζ/dα and ∂γ/∂α
    ds_a = diff_ (\ a -> s0 + s1*a + s2*a**2) a
    dg_a = head $ grad (\ [a, b] -> auto g0 + auto g1*a + auto g2*b + auto g3*a*b + auto  g4*b**2) [a, b]

db x y a b = (1/e) * bb b * (c b *dg_b*(x*y/(1 + h*x)) - dd_b)
  where
    -- derivatives: ∂y/∂β and dδ/dβ
    dg_b = last $ grad (\ [a, b] -> auto g0 + auto g1*a + auto g2*b + auto g3*a*b + auto g4*b**2) [a, b]
    dd_b = diff_ (\ b -> d0 + d1*b) b

-- the equation system that is passed to the solver
-- because the solver expects a function of this type,
-- a vector of doubles from the differential equation
-- functions, applied to their arguments is constructed.
-- the arguments are unpacked from a vector for the same reason.
eqSystem :: Double -> Vector Double -> Vector Double
eqSystem t vars = fromList [ dx x y a b, dy x y a b
                           , da x y a b, db x y a b ]
  where
    [x, y, a, b] = toList vars

-- the time steps for which the result is given
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
