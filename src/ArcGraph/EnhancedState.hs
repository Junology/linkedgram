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
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IMap
import qualified Data.Set as Set
import qualified Data.BitArray as BA

import qualified Numeric.LinearAlgebra as LA

import ArcGraph
import ArcGraph.Component
import ArcGraph.State

import Numeric.Matrix.Integral
import Numeric.Algebra.FreeModule as FM
import Numeric.Algebra.Frobenius as Frob
import Numeric.Algebra.Presentation
import Numeric.Algebra.Homology

{-- for debug
import Debug.Trace

traceCond :: Bool -> String -> a -> a
traceCond False _  = id
traceCond True msg = trace msg
--}

-----------------------
-- * Enhanced states --
-----------------------
class (DState ds, Ord e) => Enhancement ds e where
  listEnh :: (Alternative f) => Int -> ArcGraph -> ds -> f e
  diffEnh :: ArcGraph -> ds -> e -> FreeMod Int (ds,e)

listEStates :: (Enhancement ds e) => ArcGraph -> Int -> Int -> [(ds,e)]
listEStates ag i j = concatMap (\s -> (,) s <$> listEnh j ag s) $ listStates ag i

-----------------
-- * Instances --
-----------------
newtype MapEState pc = MEState (Map.Map pc SL2B)
  deriving (Eq, Ord, Show, Generic, NFData)

instance (DState ds, PComponent pc) => Enhancement ds (MapEState pc) where
  listEnh j ag st = 
    let comps = getComponents (smoothing ag st)
        deg' = j - degree ag st + L.length comps
        subs = filter (\sub -> 2*L.length sub == deg') $ L.subsequences comps
        mapS sub = Map.fromList $! map (\c -> (c, if elem c sub then SLI else SLX)) comps
    in MEState <$> foldr' (\x xs -> pure (mapS x) <|> xs) empty subs

  diffEnh ag st (MEState mp) =
    FM.sumFM $ (diffState ag st :: V.Vector (Int,ds)) >>= \dstp -> do
      let (sign,dst) = dstp
      return $! ((,) dst . MEState . Map.fromList) FM.@$>% sign FM.@*% Frob.tqftZ (doIntersect ag) (Map.toList mp) (getComponents (smoothing ag dst))

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
data KHData a ds e = KHData {
  subject :: ArcGraph,
  rank :: Int,
  tors :: [a],
  cycleV :: [FreeMod a (ds, e)],
  bndryV :: Maybe [FreeMod a (ds, e)] }
  deriving (Show,Eq,Generic, NFData)

vecToSum :: (Coefficient a, Ord b) => [b] -> Vec a -> FreeMod a b
vecToSum bs v = sumFM $! zipWith (@*@%) (vecToList v) bs

cohomologyToKH :: (Coefficient a, DState ds, Enhancement ds e) => ArcGraph -> [(ds, e)] -> Bool -> Homology a -> KHData a ds e
cohomologyToKH ag basis hasBndry hdata =
  let basisAGE = uncurry (AGraphE ag) <$> basis
  in KHData {
    subject = ag,
    rank = L.length (freeCycs hdata),
    tors = fmap fst (torsions hdata),
    cycleV = fmap (vecToSum basis) (freeCycs hdata ++ fmap snd (torsions hdata)),
    bndryV = if hasBndry
             then Just (fmap (vecToSum basis) (bndries hdata))
             else Nothing }

-- | Compute Khovanov homology for given range of cohomological degrees and a given quantum-degree
computeKhovanov :: (ChainEliminable a, NFData ds, NFData e, Show ds, Show e, Eq ds, Eq e, Enhancement ds e) => ArcGraph -> Int -> [ds] -> Bool -> IntMap (KHData a ds e)
computeKhovanov ag qdeg states hasBndry =
  let numCrs = countCross ag
      hdegs = [0..numCrs]
      slimAG = slimCross ag
      basis i = L.concatMap (\st -> (,) st <$> listEnh qdeg slimAG st) (filter ((==i). degree slimAG) states)
      !basisMap = force (IMap.fromAscList (map (\i -> (i,basis i)) hdegs))
      !diffs = flip (parMap rdeepseq) [0..numCrs-1] $ \i ->
        let sbasis = basisMap IMap.! i
            tbasis = basisMap IMap.! (i+1)
        in force $! present (uncurry (diffEnh slimAG)) sbasis tbasis
  in if numCrs == 0
     then -- The case where there is no crossing point;
       IMap.map (\v -> KHData slimAG (L.length v) [] (fmap (1@*@%) v) Nothing) basisMap
     else -- The case where there is at least one crossing point;
       let !hdata = force (homology diffs)
           !hdataMap = IMap.fromAscList $ filter (\hdti -> not (null (freeCycs (snd hdti)) && null (torsions (snd hdti))) ) $ zip hdegs hdata
       in flip IMap.mapWithKey hdataMap $ \i hdt -> cohomologyToKH slimAG (basisMap IMap.! i) hasBndry hdt
