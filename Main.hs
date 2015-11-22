import System.Random
import System.Environment
import Data.Vector as V
import Control.Monad.State
import Data.List
import Data.Map as M

import Types

spherePoints :: Double -> [Vector3]
spherePoints radius = do
  let r = ((ceiling radius)::Int) + 1
  x <- [-r .. r]
  y <- [-r .. r]
  return (x,y)
  let z_square_min = (radius - 2)^2 - (fromIntegral (x^2 + y^2))
  let z_square_max = (radius + 2)^2 - (fromIntegral (x^2 + y^2))
  let lower_bound = if z_square_min < 0 then 0 else ceiling $ sqrt z_square_min
  let upper_bound = if z_square_max < 0 then -1 else floor $ sqrt z_square_max
  abs_z <- [lower_bound .. upper_bound]
  z <- nub $ [abs_z, -abs_z]
  return $ Vector3 x y z

genSpace :: Double -> Space
genSpace radius = Data.List.foldr (\cords -> M.insert (middle + cords) Lamina) M.empty (spherePoints radius) where
  middle_point = (ceiling radius) `div` 2
  middle = Vector3 middle_point middle_point middle_point


loadChain :: FilePath -> FilePath -> IO (V.Vector Atom)
loadChain laminBS binderBS = undefined

recalculateEnergy :: SimulationState -> SimulationState
recalculateEnergy s = undefined


genSimState :: (RandomGen g) => g -> Int -> V.Vector Atom -> Space -> SimulationState
genSimState random numBinders beads space = undefined


createRandomDelta :: State SimulationState Move
createRandomDelta = undefined

energyFromDelta :: Move -> State SimulationState Double
energyFromDelta = undefined


applyDelta :: Move -> State SimulationState ()
applyDelta = undefined

simulateStep :: State SimulationState ()
simulateStep = do
        d <- createRandomDelta
        e <- energyFromDelta d
        return () -- TODO

simulate :: Int -> State SimulationState ()
simulate 0 = return ()
simulate n = simulateStep >> simulate (n - 1)

main :: IO ()
main = do
        [laminFile, binderFile, r, numBinders, steps] <- getArgs
        chain <- loadChain laminFile binderFile
        let space = genSpace $ read r
        random <- newStdGen
        let state = genSimState random (read numBinders) chain space
        let ret = execState (simulate (read steps)) state
        print ret
