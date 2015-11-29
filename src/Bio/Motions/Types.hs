{-# LANGUAGE RecordWildCards, DeriveGeneric #-}
module Bio.Motions.Types where

import System.Random
import qualified Data.Map.Strict as M
import qualified Data.Vector.Unboxed as V
import Linear as Lin
import GHC.Generics (Generic)

data Atom = NormBead | LBBead | BBBead | Lamina | Binder
          deriving (Eq, Show, Generic)

type Vector3 = Lin.V3 Int

type Space = M.Map Vector3 Atom

data SimulationState = SimulationState {
                     space :: !Space,
                     binders :: !(V.Vector Vector3),
                     beads :: !(V.Vector Vector3),
                     energy :: {-# UNPACK #-} !Double,
                     randgen :: !StdGen }

instance Show SimulationState where
        show SimulationState{..} = unlines [
                            {-show space,-}
                            show binders,
                            show beads,
                            show energy]

data Move = MoveBinder {-# UNPACK #-} !Int !Vector3 | MoveBead {-# UNPACK #-}!Int !Vector3
          deriving (Eq, Show)

data Input = Input
    { inputChainLength :: Int
    , inputLamins :: [Int]
    , inputBinders :: [Int]
    , inputRadius :: Double
    , inputNumBinders :: Int
    , inputNumSteps :: Int
    , inputRandGen :: StdGen
    }
