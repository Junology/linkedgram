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

import Control.Monad
import Control.Monad.ST

import qualified Data.Map.Strict as Map
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import Data.STRef
import Data.Foldable (for_)

import ArcGraph

-- Connection of ArcPath around a crossing according to its CrsState
localCross :: Cross -> (Segment,Segment)
localCross (Crs sega segb Crossing) = (sega,segb)
localCross (Crs (Sgmt v0 v1) (Sgmt w0 w1) Smooth0) = (Sgmt v0 w1, Sgmt v1 w0)
localCross (Crs (Sgmt v0 v1) (Sgmt w0 w1) Smooth1) = (Sgmt v0 w0, Sgmt v1 w1)

-- Attach a label on elements.
mkConnection :: (a -> Bool) -> MV.STVector s (Maybe Int,a) -> STRef s Int -> ST s ()
mkConnection p edgeMV cntRef = do
  edgeV <- V.freeze edgeMV
  cnt <- stepST cntRef
  oldLabRef <- newSTRef []
  for_ (V.findIndices (p . snd) edgeV) $ \j ->
    case fst (edgeV V.! j) of
      Just k -> modifySTRef' oldLabRef (k:)
      Nothing -> MV.modify edgeMV (\x -> (Just cnt,snd x)) j
  oldLab <- readSTRef oldLabRef
  let isOldLab (mayl,_) = maybe False (`elem` oldLab) mayl
  forM_ (V.findIndices isOldLab edgeV) $ \j ->
    MV.modify edgeMV (\x -> (Just cnt,snd x)) j
  where
    -- Step up the STRef counter and return the old value
    stepST :: (Enum a) => STRef s a -> ST s a
    stepST cntRef= do
      cnt <- readSTRef cntRef
      modifySTRef' cntRef succ
      return cnt

-- Divide the indices of ArcPath into components
components :: ArcGraph -> [[Int]]
components (AGraph ps cs) = runST $ do
  let n = length ps
  edgeMV <- MV.new n
  -- Initialize buffer
  forM_ [0..(n-1)] $ \i ->
    MV.unsafeWrite edgeMV i (Nothing, getEndVrtx (ps!!i))
  cntRef <- newSTRef 0
  forM_ cs $ \c -> do
    let (Sgmt v0 v1, Sgmt w0 w1) = localCross c
    -- Connect v0 <--> v1
    mkConnection (\x -> elem v0 x || elem v1 x) edgeMV cntRef
    -- Connect w0 <--> w1
    mkConnection (\x -> elem w0 x || elem w1 x) edgeMV cntRef
  (justcls,nocls) <- V.ifoldl' (\xs i y -> case fst y of {Just k -> ((k,[i]):fst xs,snd xs); Nothing -> (fst xs,[i]:snd xs);}) ([],[]) <$> V.unsafeFreeze edgeMV
  return $ Map.elems (Map.fromListWith (++) justcls) ++ nocls

