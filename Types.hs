{-# LANGUAGE RecordWildCards #-}
module Types where

import System.Random
import qualified Data.Map as M
import qualified Data.Vector as V
import Linear as Lin

data Atom = NormBead | LBBead | BBBead | Lamina | Binder
          deriving (Eq, Show)

type Vector3 = Lin.V3 Int

type Space = M.Map Vector3 Atom

data SimulationState = SimulationState {
                     space :: Space,
                     binders :: V.Vector Vector3,
                     beads :: V.Vector Vector3,
                     energy :: Double,
                     randgen :: StdGen }

instance Show SimulationState where
        show SimulationState{..} = unlines [
                            {-show space,-}
                            show binders,
                            show beads,
                            show energy]

data Move = MoveBinder Int Vector3 | MoveBead Int Vector3
          deriving (Eq, Show)
