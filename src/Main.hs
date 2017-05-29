import Prelude hiding (pi)

import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Diagrams(toFile)

pi :: Double
pi = 3.1415926535

signal :: [Double] -> [(Double,Double)]
signal xs = [ (x, (sin (x * pi / 45) + 1) / 2 * sin (x * pi / 5)) | x <- xs ]

main = toFile def "mychart.svg" $ do
    layout_title .= "Amplitude Modulation"
    setColors [opaque blue, opaque red]
    plot (line "am" [signal [0, 0.5 .. 400]])
    plot (points "am points" (signal [0,7..400]))
