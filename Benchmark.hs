{-# LANGUAGE RecordWildCards #-}

module Main where
import Criterion.Main
import Bio.Motions.Prototype as Prototype
import System.Random
import Control.Monad.State.Strict
import Control.DeepSeq
import Control.Exception

globalParams :: Input
globalParams = Input
    { inputChainLength = 256
    , inputLamins = [1, 8, 10, 23, 99]
    , inputBinders = [7, 9, 11, 17, 63, 89]
    , inputRadius = 30
    , inputNumBinders = 20
    , inputNumSteps = 100000
    , inputRandGen = mkStdGen 42
    }

generate :: StdGen -> Double -> Int -> [Atom] -> SimulationState
generate gen radius numBinders chain =
    let (gen', gen'') = split gen
        space = genSpace radius
    in case genSimState gen' radius numBinders chain space of
           Left _ -> generate gen'' radius numBinders chain
           Right st -> st

runSim :: Int -> SimulationState -> SimulationState
runSim steps = execState (simulate steps)

instance NFData SimulationState where
    rnf (SimulationState s bi be en gr _) =
        rnf s `seq` rnf bi `seq` rnf be `seq` rnf en `seq` rnf gr
instance NFData Atom


runBench :: ([Atom] -> SimulationState, [Atom])
         -> (SimulationState -> SimulationState, SimulationState)
         -> IO ()
runBench (gen, chain) (run, state) = defaultMain [
        bgroup "io" [
            bench "100k-chain256" $ nf Prototype.run globalParams],
        bgroup "pure" [
            bench "genstate" $ nf gen chain,
            bench "sim-100k" $ nf run state]
        ]

main :: IO ()
main = do
    let Input{..} = globalParams
    let stg ran = generate ran inputRadius inputNumBinders
    chain <- evaluate $ force $ loadChain inputChainLength inputLamins inputBinders

    chainRan <- newStdGen
    let chainTest = (stg chainRan, chain)
    -- use different random generators for different tests
    stateRan <- newStdGen
    state <- evaluate . force $ stg stateRan chain
    let runTest = (runSim (10^5), state)
    runBench chainTest runTest
