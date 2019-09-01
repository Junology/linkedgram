{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DeriveAnyClass #-}

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

import GHC.Generics (Generic)
import Control.DeepSeq (NFData)

import Control.Monad
import Control.Applicative

import Data.Foldable
import Data.Bifunctor
import qualified Data.List as L

import ArcGraph

newtype IListState = IListState [Int]
  deriving (Eq,Ord,Show,Generic,NFData)

class Ord a => DState a where
  degree :: ArcGraph -> a -> Int
  getLabel :: ArcGraph -> a -> String
  smoothing :: ArcGraph -> a -> ArcGraph
  listStates :: ArcGraph -> Int -> [a]
  -- | Compute next states with sign in differential.
  -- On determining signs, we use the ascending order on corssings.
  diffState :: Alternative t => ArcGraph -> a -> t (Int,a)

instance DState IListState where
  degree _ (IListState is) = L.length is

  getLabel (AGraph _ cs) (IListState is)
    = map (\c -> if fst c `elem` is then '1' else '0') $ zip [0..] cs

  smoothing (AGraph ps cs) (IListState is)
    = AGraph ps $ zipWith mkCrs [0..] cs
    where
      mkCrs i crs@(Crs sega segb crst)
        | crst == Crossing = Crs sega segb (if i `elem` is then Smooth1 else Smooth0)
        | otherwise = crs

  listStates (AGraph _ cs) deg =
    let crsIndices = L.findIndices isCross cs
    in IListState <$> filter ((==deg) . length) (L.subsequences crsIndices)
    where
      isCross (Crs _ _ Crossing) = True
      isCross _ = False

  diffState (AGraph _ cs) (IListState is)
    = second IListState <$> (snd $! foldr' bin (1,empty) [0..length cs - 1])
    where
      bin i (sign,sts) = if i `elem` is
                         then (negate sign, sts)
                         else (sign, pure (sign,L.sort (i:is)) <|> sts)

listSmoothing :: [Int] -> ArcGraph -> [ArcGraph]
listSmoothing degrees ag
  = smoothing ag <$!> ((listStates ag :: Int -> [IListState]) =<< degrees)
