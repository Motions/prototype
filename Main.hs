import System.Random
import System.Environment
import Data.Maybe
import Data.List
import qualified Data.Vector as V
import Control.Monad.State
import Control.Monad.Trans.Maybe
import Control.Monad.Loops
import qualified Data.Map as M
import Linear
import Data.Ord

import Types

atomMoves :: V.Vector Vector3
atomMoves = V.fromList [V3 x y z | list@[x,y,z] <- replicateM 3 [-1,0,1],
                                                   sum (map abs list) `elem` [1,2,3]]

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

getRandFromVec :: V.Vector a -> State SimulationState a
getRandFromVec vec = do
        ix <- getRandRange (0, length vec - 1)
        return $ vec V.! ix


spherePoints :: Double -> [Vector3]
spherePoints radius = do
  let r = (ceiling radius::Int) + 2
  x <- [-r .. r]
  let y_max = ceiling $ sqrt $ (radius + 2)^2 - (fromIntegral x^2)
  y <- [-y_max .. y_max]
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
            middle_point = 0
            middle = V3 middle_point middle_point middle_point

fillGaps :: Int -> [(Int, Atom)] -> [Atom]
fillGaps _ [] = []
fillGaps x ((y, a):tl) | x == y = a : fillGaps (x+1) tl
fillGaps x l = NormBead : fillGaps (x+1) l

loadChain :: Int -> FilePath -> FilePath -> IO [Atom]
loadChain length laminBS binderBS = do
  laminBSs <- flip zip (repeat LBBead) . map read . lines <$> readFile laminBS
  binderBSs <- flip zip  (repeat BBBead) . map read . lines <$> readFile binderBS
  let indexed_bss = sortBy (comparing fst) $ laminBSs ++ binderBSs ++ [(length + 1, NormBead)]
  return $ init $ fillGaps 1 indexed_bss


recalculateEnergy :: SimulationState -> SimulationState
recalculateEnergy state = state { energy = V.sum $ (V.map <$> (localEnergy . space) <*> beads) state }

genSimState :: StdGen -> Double -> Int -> [Atom] -> Space -> SimulationState
genSimState randGen radius numBinders (b:beads) space =
        let fstBead = V3 0 0 0
            st = SimulationState (M.insert fstBead b space) V.empty (V.singleton fstBead) 0 randGen
            st' = genBeads beads st
            st'' = genBinders radius numBinders st'
        in recalculateEnergy st''

genBeads :: [Atom] -> SimulationState -> SimulationState
genBeads [] st = st
genBeads (b:bs) st@(SimulationState space binders beads energy randGen) =
        let (pos, newRandGen) = tryGen 100 randGen
            newSpace = M.insert pos b space
            newBeads = V.snoc beads pos
        in genBeads bs (SimulationState newSpace binders newBeads energy newRandGen)
    where tryGen 0 _ = error "Unable to find initialization"
          tryGen n gen =
              let (ix, gen') = randomR (0, length atomMoves - 1) gen
                  delta = atomMoves V.! ix
                  lastPos = V.last beads
                  newPos = lastPos + delta
              in if collides newPos st || intersectsChain lastPos newPos st
                     then tryGen (n - 1) gen'
                     else (newPos, gen')

genBinders :: Double -> Int -> SimulationState -> SimulationState
genBinders radius n st = flip execState st $ replicateM_ n $ tryGen 100
    where tryGen :: Int -> State SimulationState ()
          tryGen 0 = fail "Unable to find initialization"
          tryGen n = do
              x <- getRandRange (-r, r)
              y <- getRandRange (-r, r)
              z <- getRandRange (-r, r)
              let v = V3 x y z
                  d = dist v (V3 0 0 0)
              st <- get
              if d > fromIntegral (r - 2) || collides v st
                  then tryGen (n - 1)
                  else let space' = M.insert v Binder $ space st
                           binders' = V.snoc (binders st) v
                       in put st { space = space', binders = binders' }
          r = ceiling radius :: Int


createRandomDelta :: MaybeT (State SimulationState) Move
createRandomDelta = do
        moveBinder <- lift getRand
        atoms <- gets $ if moveBinder then binders else beads
        atomIx <- lift $ getRandRange (0, length atoms - 1)
        delta <- lift $ getRandFromVec atomMoves
        let move = (if moveBinder then MoveBinder else MoveBead) atomIx delta
        st <- get
        guard $ not $ moveCollides move st
        if moveBinder
            then return move
            else do guard $ not $ moveBreaksChain move st
                    guard $ not $ moveIntersectsChain move st
                    return move

collides :: Vector3 -> SimulationState -> Bool
collides pos = M.member pos . space

intersectsChain :: Vector3 -> Vector3 -> SimulationState -> Bool
intersectsChain b1@(V3 x1 y1 z1) b2@(V3 x2 y2 z2) st =
        let d = dist b1 b2 -- assume 0 < d <= sqrt 2
        in d == 1 || (let crossPositions =
                              case (x1 == x2, y1 == y2, z1 == z2) of
                                  (True, _, _) -> (V3 x1 y1 z2, V3 x1 y2 z1)
                                  (_, True, _) -> (V3 x1 y1 z2, V3 x2 y1 z1)
                                  (_, _, True) -> (V3 x1 y2 z1, V3 x2 y1 z1)
                                  (_, _, _   ) -> error "d > sqrt 2"
                      in areChainNeighbors crossPositions)
    where areChainNeighbors (fstPos, sndPos) =
              let (a1, a2) = (M.lookup fstPos (space st), M.lookup sndPos (space st))
                  chainAtoms = [Just NormBead, Just LBBead, Just BBBead]
              in all (`elem` chainAtoms) [a1, a2]
                 && (let chain = beads st
                     in case V.elemIndex fstPos chain of
                            Nothing -> error "bead in space but not in chain"
                            Just ix -> sndPos `elem` [chain V.! (ix - 1) | ix > 0]
                                                  ++ [chain V.! (ix + 1) | ix < length chain - 1])

moveCollides :: Move -> SimulationState -> Bool
moveCollides (MoveBinder ix delta) st =
        let binderPos = (V.! ix) . binders $ st
            newPos = binderPos + delta
        in collides newPos st
moveCollides (MoveBead ix delta) st =
        let beadPos = (V.! ix) . beads $ st
            newPos = beadPos + delta
        in collides newPos st

moveBreaksChain :: Move -> SimulationState -> Bool
moveBreaksChain (MoveBinder _ _) _ = False
moveBreaksChain move@(MoveBead ix delta) st = any badNeighbors $ localNeighbors move st
    where badNeighbors (b1, b2) = let d = dist b1 b2 in d <= 0 || d > sqrt 2

moveIntersectsChain :: Move -> SimulationState -> Bool
moveIntersectsChain (MoveBinder _ _) _ = False
moveIntersectsChain move@(MoveBead ix delta) st =
        any (\(b1, b2) -> intersectsChain b1 b2 st) $ localNeighbors move st

localNeighbors :: Move -> SimulationState -> [(Vector3, Vector3)]
localNeighbors (MoveBinder _ _) _ = []
localNeighbors (MoveBead ix delta) st =
        let chain = beads st
            chainLen = length chain
            localBeads = [chain V.! (ix - 1) | ix > 0]
                      ++ [(chain V.! ix) + delta]
                      ++ [chain V.! (ix + 1) | ix < chainLen - 1]
        in zip localBeads (tail localBeads)

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
simulateStep = gets energy >>= selectMove >>= applyDelta
    where
        selectMove oldEnergy = loop
            where loop = do
                    m <- findDelta
                    r <- gets (energyFromDelta m) >>= checkE oldEnergy
                    if r then return m
                         else loop
        findDelta = liftM fromJust $ iterateWhile isNothing $ runMaybeT createRandomDelta
        -- TODO zastanowić się: stałe numeryczne (_DELTA) się typują na Double, bo MonomorphismRestriction
        -- albo wyłączamy, albo wymuszamy Double.
        {-checkE :: (Floating f) => f -> f -> State SimulationState Bool-}
        checkE oldEnergy newEnergy
            | newEnergy >= oldEnergy = return True
            | otherwise = do
                r <- getRandRange (0, 1)
                return $ r < exp (newEnergy - oldEnergy) * _DELTA
        _DELTA = 2 :: Double


simulate :: Int -> State SimulationState ()
simulate = flip replicateM_ simulateStep

main :: IO ()
main = do
        [chain_length, laminFile, binderFile, r, numBinders, steps] <- getArgs
        chain <- loadChain (read chain_length) laminFile binderFile
        let radius = read r
            space = genSpace radius
        randGen <- newStdGen
        let state = genSimState randGen radius (read numBinders) chain space
            ret = execState (simulate (read steps)) state
        print ret
