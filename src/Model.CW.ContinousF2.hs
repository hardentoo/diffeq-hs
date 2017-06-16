{-# LANGUAGE RecordWildCards #-}

import Data.Time.Clock (getCurrentTime)
import Data.Time.Format (defaultTimeLocale, formatTime)

import qualified Graphics.Rendering.Chart.Backend.Cairo as BC
import Graphics.Rendering.Chart.Easy ((.=), def, layout_title, plot, line)
import Numeric.LinearAlgebra (Vector, Matrix, linspace, toList, toLists, fromList, toColumns)
import Numeric.GSL.ODE
import Numeric.AD

type D = Double

αMin, αMax, βMin, βMax, ζ0, ζ1, ζ2, k0, k1, γ0, γ1, γ2, γ3, γ4, c0, c1, δ0, δ1, h, e :: Double
αMin = 0.1
αMax = 1.0
βMin = 0.1
βMax = 1.0
ζ0 = 7.0
ζ1 = 2.0
ζ2 = 1.0
k0 = 1.0
k1 = 0.01
γ0 = 0.8
γ1 = 0.8
γ2 = 2.2
γ3 = -0.3
γ4 = 0.3
c0 = 1.0
c1 = 0.01
δ0 = 0.486
δ1 = 1.052
h  = 1.0
e  = 1.0

-- initial conditions
initVals :: Vector D
initVals = fromList [ 0.5, 0.5, 0.5, 0.5 ]

-- the differential equations themselves
-- TODO: some equations can be factored out and
-- TODO: differentiation related boilerplate reduced

dx :: D -> D -> D -> D -> D
dx x y α β = x*ζ*(1 - x/k) - (γ*x*y / (1 + h*x))
  where
    ζ = ζ0 + ζ1*α + ζ2*α**2
    k = k0 + k1*α
    γ = γ0 + γ1*α + γ2*β + γ3*α*β + γ4*β**2

dy :: D -> D -> D -> D -> D
dy x y α β = c * (γ*x*y / (1 + h*x)) - y*δ
  where
    c = c0 + c1*β
    γ = γ0 + γ1*α + γ2*β + γ3*α*β + γ4*β**2
    δ = δ0 + δ1*β

dα :: D -> D -> D -> D -> D
dα x y α β = (1/e) * a * (dζ_α * (1 - x/k) + dγ_α * (y / (1 + h*x)))
  where
    a = (α - αMin)*(αMax - α)
    k = k0 + k1*α
    -- derivatives: dζ/dα and ∂γ/∂α
    dζ_α = diff (\α -> auto ζ0 + auto ζ1*α + auto ζ2*α**2) α
    dγ_α = head $ grad (\[α, β] -> auto γ0 + auto γ1*α + auto γ2*β + auto γ3*α*β + auto  γ4*β**2) [α, β]

dβ :: D -> D -> D -> D -> D
dβ x y α β = (1/e) * b * (c * dγ_β * (x*y/(1 + h*x)) - dδ_β)
  where
    b = (β - βMin)*(βMax - β)
    c = c0 + c1*β
    k = k0 + k1*α
    -- derivatives: ∂y/∂β and dδ/dβ
    dγ_β = last $ grad (\[α, β] -> auto γ0 + auto γ1*α + auto γ2*β + auto γ3*α*β + auto γ4*β**2) [α, β]
    dδ_β = diff (\β -> auto δ0 + auto δ1*β) β

-- HOW CAN I MAKE ALL FUNCTIONS NOT USE GLOBALS?
eqSystem :: D -> Vector D -> Vector D
eqSystem t vars = fromList [ dx x y α β
                           , dy x y α β
                           , dα x y α β
                           , dβ x y α β ]
  where
    [x, y, α, β] = toList vars

-- the time (steps)
times :: Vector D
times = linspace 5000 (0, 499 :: D)

-- the solutions matrix
solution :: Matrix D
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
makePlottableTuples :: Vector D -> Matrix D -> [[(D, D)]]
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
  plot $ line "prey - predator" [ map (\[x, y, _, _] -> (x, y)) $ toLists solution ]

writePlot filePath plot = BC.toFile def {BC._fo_format=BC.PDF} filePath plot

main :: IO ()
main = do
  timeStr <- getNowTimeString
  writePlot ("plots/CW_timePlot_" ++ timeStr ++ ".pdf") timePlot
  writePlot ("plots/CW_phasePlot_" ++ timeStr ++ ".pdf") phasePlot
