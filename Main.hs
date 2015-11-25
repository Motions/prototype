import System.Random
import System.Environment
import Data.Maybe
import Data.List
import qualified Data.Vector as V
import Control.Monad.State
import Control.Monad.Trans.Maybe
import qualified Data.Map as M
import Linear

import Types

atomMoves :: [Vector3]
atomMoves = [V3 x y z | list@[x,y,z] <- replicateM 3 [-1,0,1], sum (map abs list) `elem` [1,2,3]]

atomMovesLen :: Int
atomMovesLen = length atomMoves

getRandWith :: (StdGen -> (a, StdGen)) -> State SimulationState a
getRandWith func = do
    st <- get
    let randGen = randgen st
        (randVal, newRandGen) = func randGen
    put st { randgen = newRandGen }
    return randVal

getRand :: Random a => State SimulationState a
getRand = getRandWith random

getRandRange :: Random a => (a, a) -> State SimulationState a
getRandRange range = getRandWith $ randomR range

spherePoints :: Double -> [Vector3]
spherePoints radius = do
  let r = (ceiling radius :: Int) + 2
  x <- [-r .. r]
  let y_max = ceiling $ sqrt $ (radius + 2)^2 - (fromIntegral x^2)
  y <- [-y_max .. y_max]
  return (x, y)
  let z_square_min = (radius - 2)^2 - fromIntegral (x^2 + y^2)
  let z_square_max = (radius + 2)^2 - fromIntegral (x^2 + y^2)
  let lower_bound = if z_square_min < 0 then 0 else ceiling $ sqrt z_square_min
  let upper_bound = if z_square_max < 0 then -1 else floor $ sqrt z_square_max
  abs_z <- [lower_bound .. upper_bound]
  z <- nub [abs_z, -abs_z]
  return $ V3 x y z

genSpace :: Double -> Space
genSpace radius =
    foldr (\cords -> M.insert (middle + cords) Lamina) M.empty (spherePoints radius)
        where
            middle_point = ceiling radius `div` 2
            middle = V3 middle_point middle_point middle_point


loadChain :: FilePath -> FilePath -> IO (V.Vector Atom)
loadChain laminBS binderBS = undefined

recalculateEnergy :: SimulationState -> SimulationState
recalculateEnergy state = state { energy = V.sum $ (V.map <$> (localEnergy . space) <*> beads) state }

genSimState :: (RandomGen g) => g -> Int -> V.Vector Atom -> Space -> SimulationState
genSimState randGen numBinders beads space = undefined


createRandomDelta :: MaybeT (State SimulationState) Move
createRandomDelta = do
        moveBinder <- lift getRand
        atoms <- gets $ if moveBinder then binders else beads
        atomIx <- lift $ getRandRange (0, length atoms - 1)
        whichMove <- lift $ getRandRange (0, atomMovesLen - 1)
        let delta = atomMoves !! whichMove
            move = (if moveBinder then MoveBinder else MoveBead) atomIx delta
        st <- get
        guard $ not $ collides move st
        if moveBinder
            then return move
            else do guard $ not $ breaksChain move st
                    guard $ not $ intersectsChain move st
                    return move


collides :: Move -> SimulationState -> Bool
collides (MoveBinder ix delta) st =
        let binderPos = (V.! ix) . binders $ st
            newPos = binderPos + delta
        in M.member newPos . space $ st
collides (MoveBead ix delta) st =
        let beadPos = (V.! ix) . beads $ st
            newPos = beadPos + delta
        in M.member newPos . space $ st

localNeighbors :: Move -> SimulationState -> [(Vector3, Vector3)]
localNeighbors (MoveBinder _ _) _ = []
localNeighbors (MoveBead ix delta) st =
        let chain = beads st
            chainLen = length chain
            localBeads = [chain V.! (ix - 1) | ix > 0]
                      ++ [(chain V.! ix) + delta]
                      ++ [chain V.! (ix + 1) | ix < chainLen - 1]
        in zip localBeads (tail localBeads)

breaksChain :: Move -> SimulationState -> Bool
breaksChain (MoveBinder _ _) _ = False
breaksChain move@(MoveBead ix delta) st = all goodNeighbor $ localNeighbors move st
    where goodNeighbor (b1, b2) = let d = dist b1 b2 in d > 0 && d <= sqrt 2


intersectsChain :: Move -> SimulationState -> Bool
intersectsChain (MoveBinder _ _) _ = False
intersectsChain move@(MoveBead ix delta) st = all goodNeighbor $ localNeighbors move st
    where goodNeighbor (b1@(V3 x1 y1 z1), b2@(V3 x2 y2 z2)) =
              let d = dist b1 b2 -- assume prior (not . breaksChain) check, i. e. 0 < d <= sqrt 2
              in d == 1 || (let crossPositions =
                                    case (x1 == x2, y1 == y2, z1 == z2) of
                                        (True, _, _) -> (V3 x1 y1 z2, V3 x1 y2 z1)
                                        (_, True, _) -> (V3 x1 y1 z2, V3 x2 y1 z1)
                                        (_, _, True) -> (V3 x1 y2 z1, V3 x2 y1 z1)
                                        (_, _, _   ) -> error "d > sqrt 2"
                            in areChainNeighbors crossPositions)
          areChainNeighbors (fstPos, sndPos) =
              let (a1, a2) = (M.lookup fstPos (space st), M.lookup sndPos (space st))
                  chainAtoms = [Just NormBead, Just LBBead, Just BBBead]
              in all (`elem` chainAtoms) [a1, a2]
                 && (let chain = beads st
                     in case V.elemIndex fstPos chain of
                            Nothing -> error "bead in space but not in chain"
                            Just ix -> sndPos `elem` [chain V.! (ix - 1) | ix > 0]
                                                  ++ [chain V.! (ix + 1) | ix < length chain - 1])

dist :: Vector3 -> Vector3 -> Double
dist u v = sqrt $ fromIntegral $ qd u v

energyFromDelta :: Move -> SimulationState -> Double
energyFromDelta move state = next - current
  where
    (from, diff) = case move of
        MoveBinder pos diff -> (binders state V.! pos, diff)
        MoveBead pos diff -> (beads state V.! pos, diff)
    to = from + diff
    spc = space state
    current = localEnergy spc from
    next = localEnergy (moveParticle from to spc) to


-- |Computes the sum of binding energies between the atom placed at 'pos' and its neighbours.
localEnergy :: Space -> Vector3 -> Double
localEnergy space pos = sum $ map energyTo $ neighbourPositions pos
  where
    energyTo pos' = fromMaybe 0 $ energyBetweenAtoms <$> atom <*> M.lookup pos' space
    energyBetweenAtoms x y
      | (x, y) `elem` bindings || (y, x) `elem` bindings = 1
      | otherwise = 0
    atom = M.lookup pos space
    bindings = [(Lamina, LBBead), (Binder, BBBead), (Binder, LBBead)]
    neighbourPositions pos = map (pos ^+^) $ [id, negated] <*> basis

moveParticle :: Vector3 -> Vector3 -> Space -> Space
moveParticle from to space = M.insert to (space M.! from) $ M.delete from space

applyDelta :: Move -> State SimulationState ()
applyDelta (MoveBinder number delta) = do
  state <- get
  let current = binders state V.! number
  energy_from_delta <- gets $ energyFromDelta (MoveBinder number delta)
  put SimulationState{
    space = moveParticle current (current + delta) (space state),
    binders = binders state V.// [(number, current + delta)],
    beads = beads state,
    energy = energy state + energy_from_delta,
    randgen = randgen state}
applyDelta (MoveBead number delta) = do
  state <- get
  let current = beads state V.! number
  energy_from_delta <- gets $ energyFromDelta (MoveBead number delta)
  put SimulationState{
    space = moveParticle current (current + delta) (space state),
    binders = binders state,
    beads = beads state V.// [(number, current + delta)],
    energy = energy state + energy_from_delta,
    randgen = randgen state}

simulateStep :: State SimulationState ()
simulateStep = do
        Just d <- runMaybeT createRandomDelta -- TODO
        e <- gets $ energyFromDelta d
        return () -- TODO

simulate :: Int -> State SimulationState ()
simulate = flip replicateM_ simulateStep

main :: IO ()
main = do
        [laminFile, binderFile, r, numBinders, steps] <- getArgs
        chain <- loadChain laminFile binderFile
        let space = genSpace $ read r
        randGen <- newStdGen
        let state = genSimState randGen (read numBinders) chain space
        let ret = execState (simulate (read steps)) state
        print ret
