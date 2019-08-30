{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DeriveAnyClass #-}

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

import Control.Monad
import Control.Monad.ST

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
--{- WORK IN PROGRESS
newtype BArray = BArray BA.BitArray
  deriving (Eq, Show, Generic)

instance NFData BArray where
  rnf (BArray _) = ()

data Enhancement t a = Enh (t a) BArray
  deriving (Eq, Show, Generic)

instance (NFData (t a)) => NFData (Enhancement t a)

--}

data ArcGraphE ds = AGraphE {
  arcGraph :: ArcGraph,
  state :: ds,
  enhLabel :: Map.Map [Int] SL2B
  } deriving (Eq,Show,Generic,NFData)

-- | Lexicographical order on ArcGraphE
instance (DState ds) => Ord (ArcGraphE ds) where
  compare (AGraphE _ st coeff) (AGraphE _ st' coeff')
    = case compare st st' of
        LT -> LT
        GT -> GT
        EQ -> compare coeff coeff'

enhancements :: (DState ds) => Int -> ArcGraph -> ds -> V.Vector (Map.Map [Int] SL2B)
enhancements deg ag st = runST $ do
  let comps = components (smoothing ag st)
      deg' = deg + L.length comps
      subis = zip [0..] $ filter (\sub -> 2*L.length sub == deg') $ L.subsequences comps
  mapMV <- MV.unsafeNew (length subis) -- Do not need to initialize the memory
  for_ subis $ \subi -> do
    let (i,sub) = subi
    MV.write mapMV i $ L.foldl' (\xs y -> Map.insert y (if elem y sub then SLI else SLX) xs) Map.empty comps
  V.unsafeFreeze mapMV

enhancedStates :: (DState ds) => Int -> ArcGraph -> ds -> V.Vector (ArcGraphE ds)
enhancedStates deg ag st
  = V.map (AGraphE ag st) $ enhancements deg ag st

enhancedStatesL :: (DState ds) => Int -> ArcGraph -> ds -> [ArcGraphE ds]
enhancedStatesL deg ag st
  = V.foldr' (\x xs -> (AGraphE ag st x):xs) [] $ enhancements deg ag st

hasIntersection :: Eq a => [a] -> [a] -> Bool
hasIntersection x y = not $ L.null (x `L.intersect` y)

differential :: (DState ds) => ArcGraphE ds -> FreeMod Int (ArcGraphE ds)
differential (AGraphE ag st coeffMap) =
  let dStVec = diffState ag st
  in FM.sumFM $ runST $ do
    imageMV <- MV.unsafeNew (V.length dStVec)
    for_ [0..(V.length dStVec -1)] $ \i -> do
      let (sign,dSt) = dStVec V.! i
          imageAg = AGraphE ag dSt FM.@$>% Map.fromList FM.@$>% Frob.tqftZ hasIntersection (Map.toList coeffMap) (components (smoothing ag dSt))
      MV.write imageMV i (sign FM.@*% imageAg)
    V.unsafeFreeze imageMV

----------------------------------------------------------
-- The computation of (unnormalized) Khovanov homology
-----------------------------------------------------------
-- | The type to carry the data of Khovanov homologies
data KHData ds = KHData {
  rank :: Int,
  tors :: [Int],
  cycleV :: [FreeMod Int (ArcGraphE ds)],
  bndryV :: Maybe [FreeMod Int (ArcGraphE ds)] }
  deriving (Show,Eq,Generic, NFData)

vecToSum :: (Integral a, Num a, LA.Element a, Ord b) => [b] -> LA.Vector a -> FreeMod Int b
vecToSum bs v = sumFM $! zipWith (@*@%) (fromIntegral <$!> LA.toList v) bs

cohomologyToKH :: (DState ds) => [ArcGraphE ds] -> Bool -> IntHomology -> KHData ds
cohomologyToKH basis hasBndry hdata =
  KHData {
    rank = L.length (freeCycs hdata),
    tors = fmap (fromIntegral . snd) (torCycs hdata),
    cycleV = fmap (vecToSum basis) (freeCycs hdata ++ fmap fst (torCycs hdata)),
    bndryV = if hasBndry
             then Just (fmap (vecToSum basis) (bndries hdata))
             else Nothing }

-- | Compute Khovanov homology for given range of cohomological degrees and a given quantum-degree
computeKhovanov :: (DState ds) => ArcGraph -> Int -> Int -> Int -> [ds] -> Bool -> Map.Map (Int,Int) (KHData ds)
computeKhovanov ag minhdeg maxhdeg qdeg states hasBndry =
  let numCrs = countCross ag
      minhdeg' = max 0 (minhdeg-1)
      maxhdeg' = min numCrs (maxhdeg+1)
      hdegsNHD = [minhdeg' .. maxhdeg']
      slimAG = slimCross ag
      basis i = L.concatMap (enhancedStatesL (qdeg-i) slimAG) (filter ((==i). degree slimAG) states)
      basisMap = Map.fromSet basis (Set.fromList hdegsNHD)
      diffs = flip (parMap rdeepseq) [minhdeg'..maxhdeg'-1] $ \i ->
        let sbasis = basisMap Map.! i
            tbasis = basisMap Map.! (i+1)
        in force $ if null sbasis || null tbasis
                   then (length tbasis LA.>< length sbasis) []
                   else matiDataToLA $ genMatrix differential sbasis tbasis
  in if maxhdeg' <= minhdeg'
     then -- The case where there is no crossing point;
       Map.mapKeysMonotonic (\i -> (i,qdeg)) $ Map.map (\v -> KHData (L.length v) [] (fmap (1@*@%) v) Nothing) basisMap
     else -- The case where there is at least one crossing point;
       let hdata = force $ intHomology (L.head diffs) (L.tail diffs)
           hdataMap = Map.fromList $ filter (\hdti -> fst hdti >= minhdeg && fst hdti <= maxhdeg && not (null (freeCycs (snd hdti)) && null (torCycs (snd hdti))) ) $ zip [minhdeg'..maxhdeg'] hdata
           khMap = flip Map.mapWithKey hdataMap $ \i hdt -> cohomologyToKH (basisMap Map.! i) hasBndry hdt
       in Map.mapKeysMonotonic (\i -> (i,qdeg {- -2*i -})) khMap

-- | Compute Khovanov homology for given range of cohomological degrees and a given quantum-degree
computeKhovanov' :: (DState ds) => ArcGraph -> Int -> Int -> Int -> [ds] -> Bool -> Map.Map (Int,Int) (KHData ds)
computeKhovanov' ag minhdeg maxhdeg qdeg states hasBndry  =
  let numCrs = countCross ag
      hdegsNHD = [max 0 (minhdeg-1) .. min numCrs (maxhdeg+1)]
      slimAG = slimCross ag
  in runST $ do
    -- Define STRef to write the result
    resultMapRef <- newSTRef Map.empty
    -- Compute basis
    baseMVec <- MV.replicate (numCrs + 2) []
    forM_ [0..numCrs] $ \i -> do
      -- Get selected states at coh.degree i
      let statesi = filter ((==i). degree slimAG) states
      MV.write baseMVec i $ L.concatMap (enhancedStatesL (qdeg-i) slimAG) statesi
    -- Compute kernels and images of differentials
    cycleMV <- MV.replicate (numCrs+1) ([] :: [LA.Vector LA.Z])
    bndryMV <- MV.replicate (numCrs+1) ([] :: [LA.Vector LA.Z])
    diagcMV <- MV.replicate (numCrs+1) ([] :: [Int])
    forM_ hdegsNHD $ \i -> do
      base0 <- MV.read baseMVec i
      base1 <- MV.read baseMVec (i+1)
      if L.null base1
        then MV.write cycleMV i $ LA.toColumns (LA.ident (L.length base0))
        else unless (L.null base0) $ do
          let diffMat = genMatrix differential base0 base1
          let (dvec,ker,im) = kerImOf $ matiDataToLA diffMat
          MV.write cycleMV i ker
          MV.write bndryMV (i+1) im
          MV.write diagcMV (i+1) $ map fromIntegral $ LA.toList dvec
    -- Compute Homology at degree (i,qdeg)
    forM_ [minhdeg..maxhdeg] $ \i ->
      when (i <= numCrs) $ do
        base <- MV.read baseMVec i
        cycleL <- MV.read cycleMV i
        bndryL <- MV.read bndryMV i
        let freeRk = length cycleL - length bndryL
        tor <- L.filter (/=1) <$> MV.read diagcMV i
        when (freeRk > 0 || not (L.null tor)) $ do
          let cycs' = map (map fromIntegral . LA.toList) cycleL
          let bnds' = map (map fromIntegral . LA.toList) bndryL
          let (cycs'',bnds'') = (cycs' L.\\ bnds', bnds' L.\\ cycs')
          modifySTRef' resultMapRef $ Map.insert (i,qdeg) $
            KHData freeRk tor (map (flip zipSum base) cycs'') (if hasBndry then (Just $ map (flip zipSum base) bnds'') else Nothing)
    readSTRef resultMapRef
