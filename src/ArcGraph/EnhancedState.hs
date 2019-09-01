{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}

------------------------------------------------
-- |
-- Module    :  ArcGraph.EnhancedState
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Implementing manipulations on enhanced states
--
------------------------------------------------

module ArcGraph.EnhancedState where

import GHC.Generics (Generic)

import Control.DeepSeq
import Control.Parallel.Strategies

import Control.Applicative
import Control.Monad
import Control.Monad.ST

import Data.Bifunctor

import Data.STRef
import Data.Maybe as M
import Data.Foldable
import qualified Data.List as L
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import qualified Data.BitArray as BA

import qualified Numeric.LinearAlgebra as LA

import ArcGraph
import ArcGraph.Component
import ArcGraph.State

import Numeric.Algebra.FreeModule as FM
import Numeric.Algebra.Frobenius as Frob
import Numeric.Algebra.IntMatrix
import Numeric.Algebra.IntChain

{-- for debug
import Debug.Trace
--}

---------------------
-- Enhanced states --
---------------------
class (DState ds, Ord e) => Enhancement ds e where
  listEnh :: (Alternative f) => Int -> ArcGraph -> ds -> f e
  diffEnh :: ArcGraph -> ds -> e -> FreeMod Int (ds,e)

newtype MapEState pc = MEState (Map.Map pc SL2B)
  deriving (Eq, Ord, Show, Generic, NFData)

instance (DState ds, PComponent pc) => Enhancement ds (MapEState pc) where
  listEnh i ag st = 
    let comps = getComponents (smoothing ag st)
        deg' = i + L.length comps
        subs = filter (\sub -> 2*L.length sub == deg') $ L.subsequences comps
        mapS sub = Map.fromList $! map (\c -> (c, if elem c sub then SLI else SLX)) comps
    in MEState <$> foldr' (\x xs -> pure (mapS x) <|> xs) empty subs

  diffEnh ag st (MEState mp) =
    FM.sumFM $ (diffState ag st :: V.Vector (Int,ds)) >>= \dstp -> do
      let (sign,dst) = dstp
      return $! ((,) dst . MEState . Map.fromList) FM.@$>% sign FM.@*% Frob.tqftZ (hasIntersection ag) (Map.toList mp) (getComponents (smoothing ag dst))

----------------------
-- WORK IN PROGRESS --
----------------------{--
newtype BArray = BArray BA.BitArray
  deriving (Eq, Ord, Show, Generic)

instance NFData BArray where
  rnf (BArray x) = x `seq` ()

data BitEState pc = BEState pc BArray
  deriving (Eq, Show, Generic, NFData)
--}

data ArcGraphE ds e = AGraphE {
  arcGraph :: ArcGraph,
  state :: ds,
  enhancement :: e
  } deriving (Eq,Show,Generic,NFData)

enhancements :: (DState ds, PComponent pc, Alternative f) => Int -> ArcGraph -> ds -> f (Map.Map pc SL2B)
enhancements deg ag st =
  let comps = getComponents (smoothing ag st)
      deg' = deg + L.length comps
      subs = filter (\sub -> 2*L.length sub == deg') $ L.subsequences comps
      mapS sub = Map.fromList $! map (\c -> (c, if elem c sub then SLI else SLX)) comps
  in foldr' (\x xs -> pure (mapS x) <|> xs) empty subs

----------------------------------------------------------
-- The computation of (unnormalized) Khovanov homology
-----------------------------------------------------------
-- | The type to carry the data of Khovanov homologies
data KHData ds e = KHData {
  subject :: ArcGraph,
  rank :: Int,
  tors :: [Int],
  cycleV :: [FreeMod Int (ds, e)],
  bndryV :: Maybe [FreeMod Int (ds, e)] }
  deriving (Show,Eq,Generic, NFData)

vecToSum :: (Integral a, Num a, LA.Element a, Ord b) => [b] -> LA.Vector a -> FreeMod Int b
vecToSum bs v = sumFM $! zipWith (@*@%) (fromIntegral <$!> LA.toList v) bs

cohomologyToKH :: (DState ds, Enhancement ds e) => ArcGraph -> [(ds, e)] -> Bool -> IntHomology -> KHData ds e
cohomologyToKH ag basis hasBndry hdata =
  let basisAGE = uncurry (AGraphE ag) <$> basis
  in KHData {
    subject = ag,
    rank = L.length (freeCycs hdata),
    tors = fmap (fromIntegral . snd) (torCycs hdata),
    cycleV = fmap (vecToSum basis) (freeCycs hdata ++ fmap fst (torCycs hdata)),
    bndryV = if hasBndry
             then Just (fmap (vecToSum basis) (bndries hdata))
             else Nothing }

-- | Compute Khovanov homology for given range of cohomological degrees and a given quantum-degree
computeKhovanov :: (DState ds, Enhancement ds e) => ArcGraph -> Int -> Int -> Int -> [ds] -> Bool -> Map.Map (Int,Int) (KHData ds e)
computeKhovanov ag minhdeg maxhdeg qdeg states hasBndry =
  let numCrs = countCross ag
      minhdeg' = max 0 (minhdeg-1)
      maxhdeg' = min numCrs (maxhdeg+1)
      hdegsNHD = [minhdeg' .. maxhdeg']
      slimAG = slimCross ag
      basis i = L.concatMap (\st -> (,) st <$> listEnh (qdeg-i) slimAG st) (filter ((==i). degree slimAG) states)
      basisMap = Map.fromSet basis (Set.fromList hdegsNHD)
      diffs = flip (parMap rdeepseq) [minhdeg'..maxhdeg'-1] $ \i ->
        let sbasis = basisMap Map.! i
            tbasis = basisMap Map.! (i+1)
        in force $! if null sbasis || null tbasis
                    then (length tbasis LA.>< length sbasis) []
                    else matiDataToLA $ genMatrix (uncurry (diffEnh ag)) sbasis tbasis
  in if maxhdeg' <= minhdeg'
     then -- The case where there is no crossing point;
       Map.mapKeysMonotonic (\i -> (i,qdeg)) $ Map.map (\v -> KHData slimAG (L.length v) [] (fmap (1@*@%) v) Nothing) basisMap
     else -- The case where there is at least one crossing point;
       let !hdata = force $ intHomology (L.head diffs) (L.tail diffs)
           hdataMap = Map.fromList $ filter (\hdti -> fst hdti >= minhdeg && fst hdti <= maxhdeg && not (null (freeCycs (snd hdti)) && null (torCycs (snd hdti))) ) $ zip [minhdeg..maxhdeg] hdata
           khMap = flip Map.mapWithKey hdataMap $ \i hdt -> cohomologyToKH slimAG (basisMap Map.! i) hasBndry hdt
       in Map.mapKeysMonotonic (\i -> (i,qdeg {- -2*i -})) khMap
