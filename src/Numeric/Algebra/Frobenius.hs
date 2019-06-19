{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}

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

import Data.List as L (transpose)

import Numeric.Algebra.FreeModule

data SL2B = SLI | SLX
  deriving (Show, Eq)

class BFrobeniusI b where
  mult :: b -> b -> FreeMod Int b
  diag :: b -> FreeMod Int (b,b)
  unit :: FreeMod Int b
  counit :: b -> Int

-- | compare by degrees
instance Ord SL2B where
  compare SLI SLI = EQ
  compare SLI SLX = GT
  compare SLX SLI = LT
  compare SLX SLX = EQ

-- | Frobenius algebra structure over integers
instance BFrobeniusI SL2B where
  mult SLI SLI = 1@*@SLI
  mult SLI SLX = 1@*@SLX
  mult SLX SLI = 1@*@SLX
  mult SLX SLX = zeroVec

  diag SLI = 1@*@(SLI,SLX) @+ 1@*@(SLX,SLI)
  diag SLX = 1@*@(SLX,SLX)

  unit = 1@*@SLI
  counit SLI = 0
  counit SLX = 1

foldMult :: (Eq b, Ord b, BFrobeniusI b) => [b] -> FreeMod Int b
foldMult = foldl (\x y -> ((flip mult y) @=<< x)) unit

diagL :: (Eq b, Ord b, BFrobeniusI b) => b -> FreeMod Int [b]
diagL = ((\x->[fst x,snd x]) @$>) . diag

foldDiag :: (Eq b, Ord b, BFrobeniusI b) => Int -> b -> FreeMod Int [b]
foldDiag n x
  | n < 0  = undefined
  | n == 0 = counit x @*@[]
  | n == 1 = 1@*@[x]
  | n > 1  = headConsMapFM diagL $ foldDiag (n-1) x

-- | TQFT operation of genus Zero represented by coincidence
tqftZ :: (Ord a2, Ord b, BFrobeniusI b)
  => (a1 -> a2 -> Bool) -> [(a1,b)] -> [a2] -> FreeMod Int [(a2,b)]
tqftZ p doms cods
  = let mkU cs i = if i < length cs
                   then case (cs!!i) of
                     True -> Right (1@*@)
                     False -> Left unit
                   else
                     Right (1@*@)
        bin (v,x)
          = let cs = map (p v) cods
            in insertMapFM (mkU cs) $ foldDiag (sum $ map fromEnum cs) x
        med = L.transpose @$> (tensorFM $ fmap bin doms)
    in zip cods @$> mapFM foldMult med
