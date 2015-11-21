{-# LANGUAGE RecordWildCards #-}
module Types where

import System.Random
import Data.Map as M
import Data.Vector as V

data Atom = NormBead | LBBead | BBBead | Lamina | Binder
          deriving (Eq, Show)

data Vector3 = Vector3 Int Int Int
             deriving (Eq, Show, Ord)
instance Num Vector3 where
        --- TODO
        --

type Space = M.Map Vector3 Atom

data SimulationState = SimulationState {
                     space :: Space,
                     binders :: V.Vector Vector3,
                     beads :: V.Vector Vector3,
                     energy :: Double,
                     random :: StdGen }

instance Show SimulationState where
        show SimulationState{..} = unlines [
                            {-show space,-}
                            show binders,
                            show beads,
                            show energy]

data Move = MoveBinder Int Vector3 | MoveBead Int Vector3
          deriving (Eq, Show)
