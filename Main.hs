import System.Random
import System.Environment
import Data.Vector as V
import Control.Monad.State

import Types

genSpace :: Double -> Space
genSpace radius = undefined


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
