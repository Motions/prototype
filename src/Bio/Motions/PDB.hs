{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
module Bio.Motions.PDB(writePDB, Header(..)) where

import qualified Bio.PDB.EventParser.PDBEvents as PE
import qualified Bio.PDB.EventParser.PDBEventPrinter as PP
import qualified Data.ByteString.Char8 as BS
import Data.Foldable
import Data.MonoTraversable
import qualified Data.Map.Strict as M
import qualified Data.Vector.Unboxed as V
import Linear
import System.IO

import Bio.Motions.Types

data Header = Header
    { headerSequenceNumber :: Int
    , headerStep :: Int
    , headerTitle :: String
    }

writePDB :: Handle -> Header -> SimulationState -> IO ()
writePDB handle Header{..} SimulationState{..} =
    writeHeader >> doWrite 0 beads >> doWrite chainLen binders >> writeConect (chainLen - 1)
  where
    writeHeader = PP.print handle header >> PP.print handle title
    writeConect n = traverse_ (\i -> PP.print handle $ PE.CONECT [i, i+1]) [1..n]
    doWrite offset = traverse_ (PP.print handle . atomMap) . getAtoms offset
    
    getAtoms :: Int -> V.Vector Vector3 -> [(Int, Vector3, Atom)]
    getAtoms offset = zipWith (\i pos -> (i, pos, space M.! pos)) [offset..] . otoList

    atomMap (i, V3 x y z, atom) = PE.ATOM {
      no = i,
      atomtype = getName atom,
      restype = getRes atom,
      chain = ' ', --na pewno puste
      resid = i, --oni tu mają drugi raz at_nr, trochę dziwnie
      resins = ' ', --chyba
      altloc = ' ', --na pewno puste
      coords = PE.Vector3 (fromIntegral x) (fromIntegral y) (fromIntegral z),
      occupancy = 0,  --to i następne to te 2 zera u nich na końcu
      bfactor = 0,  --to jest 'tempFactor' z PDB spec, ustawiają
      segid = "",   --te 3 rzeczy u nich w ogóle nie istnieją
      elt = "",
      charge = "",
      hetatm = False --musi być false
      }

    getName Binder = "O"
    getName Lamina = "P"
    getName _ = "C"   --beads

    getRes BBBead = "BOU"
    getRes LBBead = "LAM"
    getRes Binder = "BIN"
    getRes Lamina = "LAM"
    getRes NormBead = "UNB"

    chainLen = olength beads

    header = PE.HEADER
        (BS.pack $ show headerSequenceNumber)
        (BS.pack $ show $ olength beads + olength binders)
        (BS.pack $ "step " ++ show headerStep)
    title = PE.TITLE 0 $ BS.pack headerTitle
