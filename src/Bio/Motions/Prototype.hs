{-# LANGUAGE TupleSections, RecordWildCards, OverloadedStrings, FlexibleContexts, ViewPatterns #-}
module Bio.Motions.Prototype(
    module Types,
    module PDB,
    simulate,
    run,
    genSimState,
    genSpace,
    loadChain,
    simulateStep,
    initialize) where
import System.Random
import Data.Maybe
import Data.List
import qualified Data.Vector.Unboxed as V
import Control.Monad.State.Strict
import Control.Monad.Trans.Maybe
import Control.Monad.Except
import qualified Data.Map.Strict as M
import Linear
import Data.Ord
import Data.MonoTraversable

import Bio.Motions.Types as Types
import Bio.Motions.PDB as PDB

atomMoves :: V.Vector Vector3
atomMoves = V.fromList [V3 x y z | list@[x,y,z] <- replicateM 3 [-1,0,1],
                                                   sum (map abs list) `elem` [1,2]]

getRandWith :: MonadState SimulationState m => (StdGen -> (a, StdGen)) -> m a
getRandWith func = do
    st <- get
    let randGen = randgen st
        (randVal, newRandGen) = func randGen
    put st { randgen = newRandGen }
    return randVal

getRand :: (MonadState SimulationState m, Random a) => m a
getRand = getRandWith random

getRandRange :: (MonadState SimulationState m, Random a) => (a, a) -> m a
getRandRange range = getRandWith $ randomR range

getRandFromVec :: (MonadState SimulationState m, V.Unbox a) => V.Vector a -> m a
getRandFromVec vec = do
        idx <- getRandRange (0, olength vec - 1)
        return $ vec V.! idx


spherePoints :: Double -> [Vector3]
spherePoints radius = do
  let r = (ceiling radius::Int) + 2
  x <- [-r .. r]
  let y_max = ceiling $ sqrt $ sq (radius + 2) - fromIntegral (sq x)
  y <- [-y_max .. y_max]
  let z_square_min = sq (radius - 2) - fromIntegral (sq x + sq y)
  let z_square_max = sq (radius + 2) - fromIntegral (sq x + sq y)
  let lower_bound = if z_square_min < 0 then 0 else ceiling $ sqrt z_square_min
  let upper_bound = if z_square_max < 0 then -1 else floor $ sqrt z_square_max
  abs_z <- [lower_bound .. upper_bound]
  z <- nub [abs_z, -abs_z]
  return $ V3 x y z
  where sq x = x * x

genSpace :: Double -> Space
genSpace radius = M.fromList $ map (,Lamina) $ spherePoints radius

fillGaps :: Int -> [(Int, Atom)] -> [Atom]
fillGaps _ [] = []
fillGaps x ((y, a):tl) | x == y = a : fillGaps (x+1) tl
fillGaps x l = NormBead : fillGaps (x+1) l

loadChain :: Int -> [Int] -> [Int] -> [Atom]
loadChain len laminBS binderBS = init $ fillGaps 1 indexed_bss
  where
    laminBSs = map (,LBBead) laminBS
    binderBSs = map (,BBBead) binderBS
    indexed_bss = sortBy (comparing fst) $ laminBSs ++ binderBSs ++ [(len + 1, NormBead)]


recalculateEnergy :: SimulationState -> SimulationState
recalculateEnergy st@SimulationState{..} = st {
    energy = V.sum $ V.map (localEnergy space) beads,
    gyrationRadius = coef * V.sum (V.map prefix $ V.indexed beads) }
  where
    prefix (idx, bead) = V.sum $ V.map (dist bead) $ V.unsafeTake idx beads
    coef = 2 / fromIntegral (V.length beads * (V.length beads - 1))

genSimState :: MonadError String m => StdGen -> Double -> Int -> [Atom] -> Space -> m SimulationState
genSimState randGen radius numBinders (b:beads) space = recalculateEnergy <$> st'
  where
    st = SimulationState {
        space = M.insert zero b space,
        binders = V.empty,
        beads = V.singleton zero,
        energy = 0,
        gyrationRadius = 0,
        randgen = randGen }
    st' = flip execStateT st $ genBeads beads >> genBinders radius numBinders
genSimState _ _ _ [] _ = throwError "Empty beads list"

maxGenRetries :: Int
maxGenRetries = 100

genBeads :: (MonadState SimulationState m, MonadError String m) => [Atom] -> m ()
genBeads [] = return ()
genBeads (b:bs) = do
        pos <- tryGen maxGenRetries
        st@SimulationState{..} <- get
        let newSpace = M.insert pos b space
            newBeads = V.snoc beads pos
        put $ st { space = newSpace, beads = newBeads } 
        genBeads bs
    where tryGen 0 = throwError "Unable to find initialization (beads)"
          tryGen n = do
              idx <- getRandRange (0, olength atomMoves - 1)
              st <- get
              let lastPos = V.last . beads $ st
                  delta = atomMoves V.! idx
                  newPos = lastPos + delta
              if collides newPos st || intersectsChain lastPos newPos st
                  then tryGen (n - 1)
                  else return newPos

genBinders :: (MonadState SimulationState m, MonadError String m) => Double -> Int -> m ()
genBinders radius n0 = replicateM_ n0 $ tryGen maxGenRetries
    where tryGen 0 = throwError "Unable to find initialization (binders)"
          tryGen n = do
              [x, y, z] <- replicateM 3 $ getRandRange (-r, r)
              let v = V3 x y z
              st <- get
              if dist v zero > fromIntegral (r - 2) || collides v st
                  then tryGen (n - 1)
                  else let space' = M.insert v Binder $ space st
                           binders' = V.snoc (binders st) v
                       in put st { space = space', binders = binders' }
          r = ceiling radius :: Int


createRandomDelta :: MonadState SimulationState m => MaybeT m Move
createRandomDelta = do
        moveBinder <- lift getRand
        atoms <- gets $ if moveBinder then binders else beads
        atomIdx <- lift $ getRandRange (0, olength atoms - 1)
        delta <- lift $ getRandFromVec atomMoves
        let move = (if moveBinder then MoveBinder else MoveBead) atomIdx delta
        st <- get
        guard $ not $ moveCollides move st
        if moveBinder
            then return move
            else do guard $ not $ moveBreaksChain move st
                    guard $ not $ moveIntersectsChain move st
                    return move

-- |Checks whether a point is already occupied by some particle.
collides :: Vector3 -> SimulationState -> Bool
collides pos = M.member pos . space

intersectsChain :: Vector3 -> Vector3 -> SimulationState -> Bool
intersectsChain b1@(V3 x1 y1 z1) b2@(V3 x2 y2 z2) SimulationState{..} = d /= 1 && cpOnChain && cpNeighbours
  where
    d = qd b1 b2 -- assume 0 < d <= 2
    crossPositions@[fstCrossPos, sndCrossPos]
        | x1 == x2 = [V3 x1 y1 z2, V3 x1 y2 z1]
        | y1 == y2 = [V3 x1 y1 z2, V3 x2 y1 z1]
        | z1 == z2 = [V3 x1 y2 z1, V3 x2 y1 z1]
        | otherwise = error "d > sqrt 2"
    cpOnChain = all (`elem` [Just NormBead, Just LBBead, Just BBBead]) ((`M.lookup` space) <$> crossPositions)
    cpNeighbours = case V.elemIndex fstCrossPos beads of
        Nothing -> error "bead in space but not in chain"
        Just ix -> sndCrossPos `elem` [beads V.! (ix - 1) | ix > 0]
                                   ++ [beads V.! (ix + 1) | ix < olength beads - 1]

-- |Checks whether a move would cause a collision
moveCollides :: Move -> SimulationState -> Bool
moveCollides move st = collides (snd $ moveEndPoints move st) st

moveBreaksChain :: Move -> SimulationState -> Bool
moveBreaksChain (MoveBinder _ _) _ = False
moveBreaksChain move@(MoveBead _ _) st = any badNeighbors $ localNeighbors move st
    where badNeighbors (b1, b2) = let d = dist b1 b2 in d <= 0 || d > sqrt 2

moveIntersectsChain :: Move -> SimulationState -> Bool
moveIntersectsChain (MoveBinder _ _) _ = False
moveIntersectsChain move@(MoveBead _ _) st = any (\(b1, b2) -> intersectsChain b1 b2 st) $ localNeighbors move st

localNeighbors :: Move -> SimulationState -> [(Vector3, Vector3)]
localNeighbors (MoveBinder _ _) _ = []
localNeighbors (MoveBead idx delta) st =
        let chain = beads st
            chainLen = olength chain
            localBeads = [chain V.! (idx - 1) | idx > 0]
                      ++ [(chain V.! idx) + delta]
                      ++ [chain V.! (idx + 1) | idx < chainLen - 1]
        in zip localBeads (tail localBeads)

-- |Returns the Euclidean distance between two vectors.
dist :: Vector3 -> Vector3 -> Double
dist u v = sqrt $ fromIntegral $ qd u v

-- |Returns the previous and next positions of the moved atom.
moveEndPoints :: Move -> SimulationState -> (Vector3, Vector3)
moveEndPoints move st = (from, from ^+^ delta)
  where
    (from, delta) = case move of
        MoveBinder idx del -> (binders st V.! idx, del)
        MoveBead idx del -> (beads st V.! idx, del)

-- |Computes the energy gain caused by a move (may be negative).
energyFromDelta :: Move -> SimulationState -> Double
energyFromDelta move st = nxt - cur
  where
    (from, to) = moveEndPoints move st
    spc = space st
    cur = localEnergy spc from
    nxt = localEnergy (moveParticle from to spc) to


-- |Computes the sum of binding energies between the atom placed at 'pos' and its neighbours.
localEnergy :: Space -> Vector3 -> Double
localEnergy space pos = sum $ map energyTo neighbourPositions
  where
    energyTo pos' = fromMaybe 0 $ energyBetweenAtoms <$> atom <*> M.lookup pos' space
    energyBetweenAtoms x y
      | (x, y) `elem` bindings || (y, x) `elem` bindings = 1
      | otherwise = 0
    atom = M.lookup pos space
    bindings = [(Lamina, LBBead), (Binder, BBBead), (Binder, LBBead)]
    neighbourPositions = map (pos ^+^) $ [id, negated] <*> basis

moveParticle :: Vector3 -> Vector3 -> Space -> Space
moveParticle from to space = M.insert to (space M.! from) $ M.delete from space

applyDelta :: MonadState SimulationState m => Move -> m ()
applyDelta move = do
    st <- get
    let (from, to) = moveEndPoints move st
    let st' = st {
        space = moveParticle from to (space st),
        energy = energy st + energyFromDelta move st
    }
    case move of
        MoveBinder idx _ -> put $ st' { binders = binders st' V.// [(idx, to)] }
        MoveBead idx _ -> put $ st' {
            beads = beads st' V.// [(idx, to)],
            gyrationRadius = gyrationRadius st' + gyrationRadiusDiff idx from to (beads st')}
    where
        gyrationRadiusDiff idx from to beads =
            let (front, V.tail -> back) = V.splitAt idx beads
                diff bead = fromIntegral (qd to bead - qd from bead) / (dist to bead + dist from bead)
                func = V.sum . V.map diff
                coef = 2 / fromIntegral (V.length beads * (V.length beads - 1))
            in  coef * (func front + func back)

simulateStep :: MonadState SimulationState m => m ()
simulateStep = runMaybeT selectMove >>= mapM_ applyDelta
  where
    selectMove = do
        m <- createRandomDelta
        gets (energyFromDelta m) >>= checkE >>= guard
        return m
    checkE delta
        | delta >= 0 = return True
        | otherwise = do
            r <- getRandRange (0, 1)
            return $ r < exp (delta * _DELTA)
    _DELTA = 2

simulate :: MonadState SimulationState m => Int -> m ()
simulate = flip replicateM_ simulateStep

initialize :: Input -> Either String SimulationState
initialize Input{..} = runExcept st
  where
    st = genSimState inputRandGen inputRadius inputNumBinders chain space
    chain = loadChain inputChainLength inputLamins inputBinders
    space = genSpace inputRadius

run :: Input -> Either String SimulationState
run input@Input{..} = execState (replicateM_ inputNumSteps simulateStep) <$> initialize input
