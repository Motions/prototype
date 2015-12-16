{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
module Bio.Motions.Types where

import qualified Data.Map.Strict as M
import qualified Data.Vector.Unboxed as V
import Linear as Lin
import GHC.Generics (Generic)
import Control.Monad.State
import Control.Monad.Random
import qualified Data.Serialize as Serialize
import Data.Vector.Serialize()

data Atom = NormBead | LBBead | BBBead | Lamina | Binder
          deriving (Eq, Show, Generic)

instance Serialize.Serialize Atom

type Vector3 = Lin.V3 Int

type Space = M.Map Vector3 Atom

type MonadSimulation m = (MonadState SimulationState m, MonadRandom m)

data SimulationState = SimulationState
    { space :: !Space
    , binders :: !(V.Vector Vector3)
    , beads :: !(V.Vector Vector3)
    , energy :: {-# UNPACK #-} !Double
    , gyrationRadius :: {-# UNPACK #-} !Double
    } deriving Generic

instance Serialize.Serialize SimulationState

instance Show SimulationState where
        show SimulationState{..} = unlines [
                            {-show space,-}
                            show binders,
                            show beads,
                            show energy]

data Move = MoveBinder {-# UNPACK #-} !Int !Vector3 | MoveBead {-# UNPACK #-} !Int !Vector3
          deriving (Eq, Show)

data Input = Input
    { inputChainLength :: Int
    , inputLamins :: [Int]
    , inputBinders :: [Int]
    , inputRadius :: Int
    , inputNumBinders :: Int
    }
