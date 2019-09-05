{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DeriveAnyClass #-}

------------------------------------------------
-- |
-- Module    :  ArcGraph.Component
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Compute the set of path-components of link diagrams
--
------------------------------------------------

module ArcGraph.Component where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData)

import Control.Applicative
import Control.Monad
import Control.Monad.ST

import Data.Bifunctor

import Data.Bits
import Data.Word

import qualified Data.List as L

import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IntMap

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as MVU
import Data.STRef
import Data.Foldable (for_)

import ArcGraph
import Data.ChunkedBits

---------------------
-- * Type classes
---------------------
-- | Type class to represent path-components in arc graphs.
class (Ord a) => PComponent a where
  componentAt :: ArcGraph -> a -> [ArcPath]
  getComponents :: ArcGraph -> [a]
  doIntersect :: ArcGraph -> a -> a -> Bool


----------------
-- * Instances
----------------

-- | Naive instance for PComponent
newtype ArcList = AList [Int]
  deriving (Eq,Show,Ord,Generic,NFData)

instance PComponent ArcList where
  componentAt (AGraph ps _) (AList is)
    = fst <$> L.filter ((`elem` is) . snd) (zip ps [0..])
  getComponents ag
    = let (compMap,singles) = componentMap ag
      in fmap AList $ IntMap.elems compMap ++ fmap (:[]) singles
  doIntersect _ (AList x) (AList y) = not $ L.null (x `L.intersect` y)

-- | Path component in terms of bits
newtype ArcBits b = ABits b
  deriving (Eq,Show,Ord,Generic,NFData)

instance (Ord b, BitPool b) => PComponent (ArcBits b) where
  componentAt (AGraph ps _) (ABits bs)
    = fst <$> L.filter ((testBit bs) . snd) (zip ps [0..])
  getComponents ag@(AGraph ps _)
    = let (compMap,singles) = componentMap ag :: (IntMap [Int],[Int])
      in ABits . newBits (length ps) <$> IntMap.elems compMap ++ fmap (:[]) singles
  doIntersect _ (ABits xs) (ABits ys)
    = popCount (xs .&. ys) > 0


----------------------
-- * Implementaions
----------------------

-- | Connection of ArcPath around a crossing according to its CrsState
localCross :: Cross -> (Segment,Segment)
localCross (Crs sega segb Crossing) = (sega,segb)
localCross (Crs (Sgmt v0 v1) (Sgmt w0 w1) Smooth0) = (Sgmt v0 w1, Sgmt v1 w0)
localCross (Crs (Sgmt v0 v1) (Sgmt w0 w1) Smooth1) = (Sgmt v0 w0, Sgmt v1 w1)

-- | Attach a label on elements.
reattachLabel :: (a -> Bool) -> [a] -> MVU.STVector s Int -> ST s Int -> ST s ()
reattachLabel p xs labelMV genLabel = do
  let indxs = L.findIndices p xs
  oldLabs <- filter (>0) <$> (MVU.read labelMV `traverse` indxs)
  newLab <- genLabel
  for_ (L.findIndices p xs) $ \j ->
    MVU.write labelMV j newLab
  labelV <- VU.freeze labelMV
  flip VU.imapM_ labelV $ \j lab ->
    when (lab `elem` oldLabs) $ MVU.write labelMV j newLab

-- | Step up the STRef counter and return the old value
stepST :: (Enum a) => STRef s a -> ST s a
stepST cntRef= do
  cnt <- readSTRef cntRef
  modifySTRef' cntRef succ
  return cnt

-- | Label ArcPath so that two paths are with the same label provided they belong to the same path-components.
-- Return the dictionary of path-components with the indices of unlabled paths.
componentMap :: (Alternative t) => ArcGraph -> (IntMap (t Int), t Int)
componentMap (AGraph ps cs) = runST $ do
  let n = length ps
  edgeMV <- MVU.replicate n (-1)
  cntRef <- newSTRef 0
  forM_ cs $ \c -> do
    let (segA, segB) = localCross c
    -- Connect ends of segA and of segB respectively.
    reattachLabel ((segA `connectable`) . getEndSgmt) ps edgeMV (stepST cntRef)
    reattachLabel ((segB `connectable`) . getEndSgmt) ps edgeMV (stepST cntRef)
  VU.ifoldl' ibin (IntMap.empty ,empty) <$> VU.unsafeFreeze edgeMV
    where
      ibin xs i y = if y >= 0
                    then first (IntMap.insertWith (<|>) y (pure i)) xs
                    else second (pure i <|>) xs
