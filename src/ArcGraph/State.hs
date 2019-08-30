{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE DeriveGeneric #-}

------------------------------------------------
-- |
-- Module    :  ArcGraph.State
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- States on link diagrams
--
------------------------------------------------

module ArcGraph.State where

import Control.Monad
import Control.Applicative

import Data.Foldable
import qualified Data.List as L

import ArcGraph

type DiagramState = [Int]

class Ord a => DState a where
  degree :: ArcGraph -> a -> Int
  smoothing :: ArcGraph -> a -> ArcGraph
  listStates :: ArcGraph -> Int -> [a]
  -- | Compute next states with sign in differential.
  -- On determining signs, we use the ascending order on corssings.
  diffState :: Alternative t => ArcGraph -> a -> t (Int,a)

instance DState [Int] where
  degree _ = L.length

  smoothing (AGraph ps cs) st
    = AGraph ps $ zipWith mkCrs [0..] cs
    where
      mkCrs i crs@(Crs sega segb crst)
        | crst == Crossing = Crs sega segb (if i `elem` st then Smooth1 else Smooth0)
        | otherwise = crs

  listStates (AGraph _ cs) deg =
    let crsIndices = L.findIndices isCross cs
    in filter ((==deg) . length) $ L.subsequences crsIndices
    where
      isCross (Crs _ _ Crossing) = True
      isCross _ = False

  diffState (AGraph _ cs) st
    = snd $! foldr' bin (1,empty) [0..length cs - 1]
    where
      bin i (sign,sts) = if i `elem` st
                         then (negate sign, sts)
                         else (sign, pure (sign,i:st) <|> sts)

listSmoothing :: [Int] -> ArcGraph -> [ArcGraph]
listSmoothing degrees ag
  = smoothing ag <$!> ((listStates ag :: Int -> [DiagramState]) =<< degrees)
