{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}

------------------------------------------------
-- |
-- Module    :  Data.ChunkedBits
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Array of bits in terms of Word64
--
------------------------------------------------

module Data.ChunkedBits (
  BitPool,
  ChunkedBits,
  newBits
  )where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData, force)

import Control.Monad
import Control.Monad.ST

import Data.Bits
import Data.Word

import Numeric (showIntAtBase)
import Data.Char (intToDigit)

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as MVU

----------------
-- * Classes
---------------
class (Bits b) => BitPool b where
  newBits :: (Foldable t) => Int -> t Int -> b

-- * Types
-------------
type ChunkType = Word64

-- | Array of bits in terms of ChunkType == Word64
newtype ChunkedBits = ChBits (VU.Vector ChunkType)
  deriving (Eq,Ord,Generic,NFData)

----------------------------
-- * Supporting functions
----------------------------
chunkSize :: Int
chunkSize = finiteBitSize (undefined :: ChunkType)
divideInt i = (i `shiftR` 6, i .&. 0x3F) -- cf. 2^6 == 0x3F+1 == chunkSize
showBinary x = showIntAtBase 2 intToDigit x ""

liftJ2 :: (a->a->a) -> Maybe a -> Maybe a -> Maybe a
liftJ2 f Nothing ym = ym
liftJ2 f xm Nothing = xm
liftJ2 f (Just x) (Just y) = Just $! f x y

unionWith :: (Num a, VU.Unbox a) => (a->a->a) -> VU.Vector a -> VU.Vector a -> VU.Vector a
unionWith f xs ys = runST $ do
  let len = max (VU.length xs) (VU.length ys)
  stVec <- MVU.unsafeNew len
  forM_ [0..len-1] $ \i -> do
    MVU.unsafeWrite stVec i $ maybe 0 id (liftJ2 f (xs VU.!? i) (ys VU.!? i))
  VU.unsafeFreeze stVec

-----------------
-- * Instances
-----------------
instance Show ChunkedBits where
  show (ChBits xs)
    = VU.foldl (\s x -> s++showBinary x) "" xs

instance Bits ChunkedBits where
  (ChBits xs) .&. (ChBits ys)
    = ChBits $! unionWith (.&.) xs ys
  (ChBits xs) .|. (ChBits ys)
    = ChBits $! unionWith (.|.) xs ys
  (ChBits xs) `xor` (ChBits ys)
    = ChBits $! unionWith xor xs ys
  complement (ChBits xs)
    = ChBits $! VU.map complement xs
  shiftL (ChBits xs) i
    = let len = VU.length xs
          (ihigh,ilow) = divideInt i
      in ChBits $! force $ runST $ do
    stVec <- MVU.new len
    forM_ [0..(len-ihigh-1)] $ \j -> do
      let xhigh = (xs VU.! (j+ihigh)) `shiftL` ilow
          xlow  = maybe 0 (`shiftR` (chunkSize-ilow)) (xs VU.!? (j+ihigh+1))
      MVU.write stVec j (xhigh .|. xlow)
    VU.unsafeFreeze stVec
  shiftR (ChBits xs) i
    = let len = VU.length xs
          (ihigh,ilow) = divideInt i
      in ChBits $! force $ runST $ do
    stVec <- MVU.new len
    forM_ [0..(len-ihigh-1)] $ \j -> do
      let xhigh = maybe 0 (`shiftL` (64-ilow)) (xs VU.!? (j-1))
          xlow  = (xs VU.! j) `shiftR` ilow
      MVU.write stVec (ihigh+j) (xhigh .|. xlow)
    VU.unsafeFreeze stVec
  rotateL (ChBits xs) i
    = let len = VU.length xs
          (ihigh,ilow) = divideInt i
      in ChBits $! force $ runST $ do
    stVec <- MVU.new len
    forM_ [0..(len-1)] $ \j -> do
      let isrc = (j+ihigh) `mod` len
          isrcR = (j+ihigh+1) `mod` len
          xhigh = (xs VU.! isrc) `shiftL` ilow
          xlow  = (xs VU.! isrcR) `shiftR` (chunkSize-ilow)
      MVU.write stVec j (xhigh .|. xlow)
    VU.unsafeFreeze stVec
  rotateR (ChBits xs) i
    = let len = VU.length xs
          (ihigh,ilow) = divideInt i
      in ChBits $! force $ runST $ do
    stVec <- MVU.new len
    forM_ [0..(len-1)] $ \j -> do
      let isrc  = (j-ihigh) `mod` len
          isrcL = (j-ihigh-1) `mod` len
          xhigh = (xs VU.! isrc) `shiftL` ilow
          xlow  = (xs VU.! isrcL) `shiftR` (chunkSize-ilow)
      MVU.write stVec j (xhigh .|. xlow)
    VU.unsafeFreeze stVec
  bitSize _ = undefined
  -- ^ bitSize must ignore the argument, which prevent us to implement the function.
  bitSizeMaybe (ChBits xs) = Nothing
  isSigned _ = False
  testBit (ChBits xs) i
    = let (ihigh,ilow) = divideInt i
      in maybe False (\x -> testBit x ilow) (xs VU.!? ihigh)
  bit i
    = let (ihigh,ilow) = divideInt i
      in ChBits $! force $! runST $ do
    stVec <- MVU.new (ihigh+1)
    MVU.write stVec ihigh (bit ilow :: ChunkType)
    VU.unsafeFreeze stVec
  popCount (ChBits xs) = VU.foldl' (\n x -> n+popCount x) 0 xs

instance BitPool ChunkedBits where
  newBits n is = ChBits $! force $ runST $ do
    stVec <- MVU.new (n `shiftR` 6 + 1)
    forM_ is $ \i -> do
      let (ihigh,ilow) = divideInt i
      MVU.modify stVec (.|. shiftL 1 ilow) ihigh
    VU.unsafeFreeze stVec
