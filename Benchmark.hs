{-# LANGUAGE RecordWildCards #-}

module Main where
import Criterion.Main
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

generate :: StdGen -> Int -> Int -> [Atom] -> (SimulationState, StdGen)
generate gen radius numBinders chain =
    let space = genSpace radius
        (est, gen') = flip runRand gen $ runExceptT $ genSimState radius numBinders chain space
    in case est of
           Left _ -> generate gen' radius numBinders chain
           Right st -> (st, gen')

runSim :: StdGen -> Int -> SimulationState -> SimulationState
runSim gen steps st = flip evalRand gen $ flip execStateT st $ replicateM_ steps simulateStep

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
    let stg ran = generate ran inputRadius inputNumBinders
    chain <- evaluate $ force $ loadChain inputChainLength inputLamins inputBinders

    chainRan <- newStdGen
    let chainTest = (fst . stg chainRan, chain)
    -- use different random generators for different tests
    stateRan <- newStdGen
    let (st, ran) = stg stateRan chain
    st' <- evaluate . force $ st
    let runTest = (runSim ran (10^5), st')
    runBench chainTest runTest
