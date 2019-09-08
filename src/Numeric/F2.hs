{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE StandaloneDeriving, GeneralizedNewtypeDeriving #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE DefaultSignatures #-}

------------------------------------------------
-- |
-- Module    :  Numeric.F2
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Matrix presentations of morphisms
--
------------------------------------------------

module Numeric.F2 (F2) where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData)

import Data.Bits (Bits,FiniteBits,testBit)
import Data.Ratio (numerator,denominator)
import Data.Bifunctor (first)
import Foreign.Storable (Storable)

newtype F2 = F2 Bool
  deriving (Bounded,Enum,Eq,Ord,Generic,NFData,Bits,FiniteBits,Storable)

---------------------------------------
-- * Instances for Numeric classes
---------------------------------------
instance Num F2 where
  (F2 x) + (F2 y) = F2 (x /= y)
  (F2 x) * (F2 y) = F2 (x && y)
  (-) = (+)
  abs = id
  signum = id
  fromInteger x = F2 (testBit x 0)
  negate = id

instance Fractional F2 where
  fromRational x =
    let num = fromInteger (numerator x) :: F2
        den = fromInteger (denominator x) :: F2
    in num*den
  x / (F2 True) = x
  x / (F2 False) = error "Division by zero"
  recip x@(F2 True) = x
  recip (F2 False) = error "Division by zero"

---------------------
-- * Show and Read
---------------------
instance Show F2 where
  show (F2 True)  = "1"
  show (F2 False) = "0"

instance Read F2 where
  readsPrec prior = map (first fromInteger) . (readsPrec prior)
