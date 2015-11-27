module Main where
import Criterion.Main
import Types
import Prototype
import System.Random
import Control.Monad.State.Strict
import Control.DeepSeq
import Control.Exception

globalParams :: [String]
--TODO don't hardcode paths
globalParams = ["256", "../lamin_bsites.txt", "../regular_bsites.txt", "30", "200", "100000"]

-- hack
mainNoOut :: [String] -> IO SimulationState
mainNoOut args = do
        let [chain_length, laminFile, binderFile, r, numBinders, steps] = args
        chain <- loadChain (read chain_length) laminFile binderFile
        let radius = read r
            space = genSpace radius
        randGen <- newStdGen
        let state = genSimState randGen radius (read numBinders) chain space
            ret = execState (simulate (read steps)) state
        return ret

generate :: StdGen -> Double -> Int -> [Atom] -> SimulationState
generate ran radius numBinders chain =
    genSimState ran radius numBinders chain space
        where
            space = genSpace radius

runSim :: Int -> SimulationState -> SimulationState
runSim steps = execState (simulate steps)

instance NFData SimulationState where
    rnf (SimulationState s bi be en _) =
        rnf s `seq` rnf bi `seq` rnf be `seq` rnf en
instance NFData Atom


runBench :: ([Atom] -> SimulationState, [Atom])
         -> (SimulationState -> SimulationState, SimulationState)
         -> IO ()
runBench (gen, chain) (run, state) = defaultMain [
        bgroup "io" [
            bench "100k-chain256" $ nfIO (mainNoOut globalParams)],
        bgroup "pure" [
            bench "genstate" $ nf gen chain,
            bench "sim-100k" $ nf run state]
        ]

main :: IO ()
main = do
    let [chain_length, laminFile, binderFile, r, numBinders, steps] = globalParams
    let stg ran = generate ran (read r) (read numBinders)
    chain <- loadChain (read chain_length) laminFile binderFile >>= evaluate . force

    chainRan <- newStdGen
    let chainTest = (stg chainRan, chain)
    -- use different random generators for different tests
    stateRan <- newStdGen
    state <- evaluate . force $ stg stateRan chain
    let runTest = (runSim (10^5), state)
    runBench chainTest runTest
