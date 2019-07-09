{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric #-}

------------------------------------------------
-- |
-- Module    : Numeric.Algebra.Frobenius
-- Copyright : (c) Jun Yoshida 2019
-- License   : BSD3
--
-- Defining Frobenius algebra with free basis
--
------------------------------------------------

module Numeric.Algebra.Frobenius where

import GHC.Generics

import Control.Applicative

import Control.Parallel
import Control.Parallel.Strategies

import Data.List as L (transpose, foldl', splitAt, replicate)

import Numeric.Algebra.FreeModule

data SL2B = SLI | SLX
  deriving (Show, Eq, Ord, Generic)

instance NFData SL2B

-- | Frobenius algebra structure over integers
class Ord b => BFrobeniusI b where
  {-# MINIMAL (mult, unit | foldMult), (diag, counit | foldDiag) #-}
  -- Algebra structure --
  mult :: b -> b -> FreeMod Int b
  -- ^ binary multiplication

  unit :: FreeMod Int b
  -- ^ unit

  foldMult :: [b] -> FreeMod Int b
  -- ^ fold multiplication

  -- Coalgebra structure --
  diag :: b -> FreeMod Int (b,b)
  -- ^ binary comultiplication

  counit :: b -> Int
  -- ^ counit

  foldDiag :: Int -> b -> FreeMod Int [b]
  -- ^ fold comultiplication

  endoaug :: b -> FreeMod Int b
  endoaug x = counit x @*% unit

  -- Default implementations
  mult x y = foldMult [x,y]
  unit = foldMult []
  diag x = (\y -> (head y, head (tail y))) @$> foldDiag 2 x
  counit = getCoeff [] . foldDiag 0
  foldMult = foldl (\x y -> flip mult y @=<< x) unit
  foldDiag n x
    | n < 0      = undefined
    | n == 0     = counit x @*@ empty
    | n == 1     = 1@*@ pure x
    | otherwise  = diag x @>>= \y -> tensorWith (:) (1@*@% fst y) (foldDiag (n-1) (snd y))

-- | Frobenius algebra structure over integers
instance BFrobeniusI SL2B where
  mult SLI SLI = 1@*@SLI
  mult SLI SLX = 1@*@SLX
  mult SLX SLI = 1@*@SLX
  mult SLX SLX = zeroVec

  foldMult = maybe zeroVec (1@*@) . L.foldl' mult' (Just SLI)
    where
      mult' mel SLI = mel
      mult' Nothing _ = Nothing
      mult' (Just SLI) SLX = Just SLX
      mult' (Just SLX) SLX = Nothing

  diag SLI = 1@*@(SLI,SLX) @+ 1@*@(SLX,SLI)
  diag SLX = 1@*@(SLX,SLX)

  foldDiag n SLI = sumFM' $ flip map [0..(n-1)] $ \i ->
    let (xs,ys) = L.splitAt i (replicate (n-1) SLX)
    in 1@*@% (xs ++ (SLI:ys))
  foldDiag n SLX = 1@*@% L.replicate n SLX

  unit = 1@*@SLI
  counit SLI = 0
  counit SLX = 1

{-- -- Abandoned since it is much slower.
tqftZ' :: (Ord a2, BFrobeniusI b)
  => (a1 -> a2 -> Bool) -> [(a1,b)] -> [a2] -> FreeMod Int [(a2,b)]
tqftZ' p doms cods =
  let (!ds,!bs) = unzip doms
                `using` evalTuple2 (parList rpar) (parList rpar)
      !med1 = parMap rpar (foldDiag (length cods)) bs
      !cupcaps = parMap rpar (\x -> map (\y -> if p x y then (1@*@%) else endoaug) cods) ds
      !med2 = zipWith tensorMapFM cupcaps med1 `using` parList rpar
      !med = L.transpose @$>% tensorFM med2
  in zip cods @$>% mapFM foldMult med
--}

-- | TQFT operation of genus Zero represented by coincidence
tqftZ :: (Ord a2, Ord b, BFrobeniusI b)
  => (a1 -> a2 -> Bool) -> [(a1,b)] -> [a2] -> FreeMod Int [(a2,b)]
tqftZ p doms cods
  = let mkU cs i = if i < length cs
                   then if cs!!i
                     then Right (1@*@)
                     else Left unit
                   else
                     Right (1@*@)
        bin (v,x)
          = let cs = map (p v) cods
            in insertMapFM (mkU cs) $ foldDiag (sum $ map fromEnum cs) x
        med = L.transpose @$> tensorFM (parMap rpar bin doms)
    in zip cods @$> mapFM foldMult med
