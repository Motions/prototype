{-# LANGUAGE RecordWildCards #-}

module Main where
import Criterion.Main
import Data.Either
import Control.Monad.Loops
import Bio.Motions.Prototype as Prototype
import System.Random
import Control.Monad.State.Strict
import Control.DeepSeq
import Control.Exception
import Control.Monad.Except
import Control.Monad.Random

globalParams :: Input
globalParams = Input
    { inputChainLength = 256
    , inputLamins = [1, 8, 10, 23, 99]
    , inputBinders = [7, 9, 11, 17, 63, 89]
    , inputRadius = 30
    , inputNumBinders = 20
    }

generate :: MonadRandom m => Int -> Int -> [Atom] -> m SimulationState
generate radius numBinders chain = do
    Right res <- iterateUntil isRight $ runExceptT $
        genSimState radius numBinders chain space
    return res
  where
    space = genSpace radius

runSim :: MonadRandom m => Int -> SimulationState -> m SimulationState
runSim steps st = flip execStateT st $ replicateM_ steps simulateStep

instance NFData SimulationState
instance NFData Atom

runBench :: ([Atom] -> SimulationState, [Atom])
         -> (SimulationState -> SimulationState, SimulationState)
         -> IO ()
runBench (gen, chain) (run, st) = defaultMain [
        bgroup "pure" [
            bench "genstate" $ nf gen chain,
            bench "sim-100k" $ nf run st]
        ]

main :: IO ()
main = do
    let Input{..} = globalParams
    let stg = generate inputRadius inputNumBinders
    chain <- evaluate $ force $ loadChain inputChainLength inputLamins inputBinders

    chainRan <- newStdGen
    let chainTest = (flip evalRand chainRan . stg, chain)

    -- use different random generators for different tests
    stateRan <- newStdGen
    st <- evaluate . force =<< evalRandIO (stg chain)
    let runTest = (flip evalRand stateRan . runSim (10^5), st)

    runBench chainTest runTest
