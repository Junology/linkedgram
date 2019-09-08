{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
{-# LANGUAGE TypeOperators #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Matrix.F2.Mutable
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Mutable matrix with entries in F2; the prime field of characteristic 2
--
------------------------------------------------

module Numeric.Matrix.F2.Mutable where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData)

import Control.Monad.Primitive
import Data.Primitive.MutVar

import Control.Monad (when, forM_)
import Data.Bits
import Data.Word (Word8)
import Foreign.C.Types (CUInt(..), CInt(..))
import Foreign.Ptr (Ptr)
import Foreign.ForeignPtr.Unsafe (unsafeForeignPtrToPtr)
import Foreign.Marshal (fromBool)

import qualified Data.Vector.Storable  as VS
import qualified Data.Vector.Storable.Mutable as MVS

------------------
-- * Data types
------------------

-- | Mutable matrices
data MMatrixF2 s = MMatrixF2 {
  mrows :: {-# UNPACK #-} !(MutVar s Int),
  mcols :: {-# UNPACK #-} !(MutVar s Int),
  mmatData :: {-# UNPACK #-} !(MVS.MVector s Word8)
  } deriving (Generic)

type STMatrixF2 s = MMatrixF2 s
type IOMatrixF2 = MMatrixF2 RealWorld

-- | Mutable vectors
newtype MVectorF2 s = MVectorF2 { mvecData :: MVS.MVector s Word8 }
  deriving (Generic,NFData)

type STVectorF2 s = MVectorF2 s
type IOVectorF2 = MVectorF2 RealWorld

-- | The size of a mutable vector.
msizeF2 :: MVectorF2 m -> Int
msizeF2 = MVS.length . mvecData


------------------------------
-- * Construction of vectors
------------------------------
-- | Create a mutable vector of the given length.
--   All the entries are initialized to the given value.
newVectorF2 :: (Bits a, PrimMonad m) => Int -> a -> m (MVectorF2 (PrimState m))
newVectorF2 l x = MVectorF2 <$> MVS.replicate l (fromBool (testBit x 0))

-- | Create a mutable vector of the given length.
--   Entries will not be initalized.
unsafeNewVectorF2 :: (PrimMonad m) => Int -> m (MVectorF2 (PrimState m))
unsafeNewVectorF2 l = MVectorF2 <$> MVS.unsafeNew l

-- | Clone a mutable vector.
cloneVectorF2 :: (PrimMonad m) => MVectorF2 (PrimState m) -> m (MVectorF2 (PrimState m))
cloneVectorF2 mvec = MVectorF2 <$> MVS.clone (mvecData mvec)


--------------------------------
-- * Construction of matrices
--------------------------------
-- | Create a mutable matrix of the given size.
--   All the entries are initialized to the given value.
newMatrixF2 :: (Bits a, PrimMonad m) => Int -> Int -> a -> m (MMatrixF2 (PrimState m))
newMatrixF2 r c x = do
  mr <- newMutVar r
  mc <- newMutVar c
  dt <- MVS.replicate (r*c) (fromBool (testBit x 0))
  return (MMatrixF2 mr mc dt)

newIdentF2 :: (PrimMonad m) => Int -> m (MMatrixF2 (PrimState m))
newIdentF2 n = do
  mMat@(MMatrixF2 _ _ mdt) <- newMatrixF2 n n False
  forM_ [0..(n-1)] $ \i -> do
    MVS.write mdt (i*(n+1)) 1
  return mMat

-- | Create a mutable matrix of the given size.
--   Entries will not be initalized.
unsafeNewMatrixF2 :: (PrimMonad m) => Int -> Int -> m (MMatrixF2 (PrimState m))
unsafeNewMatrixF2 r c = do
  mr <- newMutVar r
  mc <- newMutVar c
  dt <- MVS.unsafeNew (r*c)
  return (MMatrixF2 mr mc dt)

-- | Clone a mutable matrix.
cloneMatrixF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> m (MMatrixF2 (PrimState m))
cloneMatrixF2 mmat = do
  mr <- newMutVar =<< readMutVar (mrows mmat)
  mc <- newMutVar =<< readMutVar (mcols mmat)
  dt <- MVS.clone (mmatData mmat)
  return (MMatrixF2 mr mc dt)


-----------------------------------------------
-- * Accessing individual entries in vectors
-----------------------------------------------
-- | Read the entry at the given position.
readVectorF2 :: (Num a, PrimMonad m) => MVectorF2 (PrimState m) -> Int -> m a
readVectorF2 mvec i = fromIntegral <$> MVS.read (mvecData mvec) i

-- | Replace the entry at the given position.
writeVectorF2 :: (Bits a, PrimMonad m) => MVectorF2 (PrimState m) -> Int -> a -> m ()
writeVectorF2 mvec i x = MVS.write (mvecData mvec) i (fromBool (testBit x 0))

-- | Modify the entry at the given position.
modifyVectorF2 :: (PrimMonad m) => MVectorF2 (PrimState m) -> (Bool -> Bool) -> Int -> m ()
modifyVectorF2 mvec f i =
  MVS.modify (mvecData mvec) (fromBool . f . flip testBit 0) i

-- | Read the entry at the given position; No bounds check are performed.
unsafeReadVectorF2 :: (Num a, PrimMonad m) => MVectorF2 (PrimState m) -> Int -> m a
unsafeReadVectorF2 mvec i = fromIntegral <$> MVS.read (mvecData mvec) i

-- | Replace the entry at the given position; No bounds check are performed.
unsafeWriteVectorF2 :: (Bits a, PrimMonad m) => MVectorF2 (PrimState m) -> Int -> a -> m ()
unsafeWriteVectorF2 mvec i x =
  MVS.write (mvecData mvec) i (fromBool (testBit x 0))

-- | Modify the entry at the given position; No bounds check are performed.
unsafeModifyVectorF2 :: (PrimMonad m) => MVectorF2 (PrimState m) -> (Bool -> Bool) -> Int -> m ()
unsafeModifyVectorF2 mvec f i =
  MVS.modify (mvecData mvec) (fromBool . f . flip testBit 0) i


-----------------------------------------------
-- * Accessing individual entries in matrices
-----------------------------------------------
-- | Read the entry at the given position.
readMatrixF2 :: (Num a, PrimMonad m) => MMatrixF2 (PrimState m) -> Int -> Int -> m a
readMatrixF2 mmat i j = do
  c <- readMutVar (mcols mmat)
  fromIntegral <$> MVS.read (mmatData mmat) (i*c + j)

-- | Replace the entry at the given position.
writeMatrixF2 :: (Bits a, PrimMonad m) => MMatrixF2 (PrimState m) -> Int -> Int -> a -> m ()
writeMatrixF2 mmat i j x = do
  c <- readMutVar (mcols mmat)
  MVS.write (mmatData mmat) (i*c + j) (fromBool (testBit x 0))

-- | Modify the entry at the given position.
modifyMatrixF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> (Bool -> Bool) -> Int -> Int -> m ()
modifyMatrixF2 mmat f i j = do
  c <- readMutVar (mcols mmat)
  MVS.modify (mmatData mmat) (fromBool . f . flip testBit 0) (i*c + j)

-- | Read the entry at the given position; No bounds check are performed.
unsafeReadMatrixF2 :: (Num a, PrimMonad m) => MMatrixF2 (PrimState m) -> Int -> Int -> m a
unsafeReadMatrixF2 mmat i j = do
  c <- readMutVar (mcols mmat)
  fromIntegral <$> MVS.read (mmatData mmat) (i*c + j)

-- | Replace the entry at the given position; No bounds check are performed.
unsafeWriteMatrixF2 :: (Bits a, PrimMonad m) => MMatrixF2 (PrimState m) -> Int -> Int -> a -> m ()
unsafeWriteMatrixF2 mmat i j x = do
  c <- readMutVar (mcols mmat)
  MVS.write (mmatData mmat) (i*c + j) (fromBool (testBit x 0))

-- | Modify the entry at the given position; No bounds check are performed.
unsafeModifyMatrixF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> (Word8 -> Word8) -> Int -> Int -> m ()
unsafeModifyMatrixF2 mmat f i j = do
  c <- readMutVar (mcols mmat)
  MVS.modify (mmatData mmat) f (i*c + j)


----------------------------
-- * Imported C functions
----------------------------
-- Copied from the source of hmatrix package.
infixr 5 ::>
type (::>) t r =  CUInt -> CUInt -> Ptr t -> r

foreign import ccall unsafe "wrap_copy_transpose" c_copyTranspose :: Word8 ::> Word8 ::> IO CInt

foreign import ccall unsafe "wrap_matmul_bin" c_matmul :: Word8 ::> Word8 ::> Word8 ::> IO CInt

foreign import ccall unsafe "wrap_elim_rows" c_elimRows :: Word8 ::> Word8 ::> Word8 ::> IO CUInt

foreign import ccall unsafe "wrap_diag_rep" c_diagRep :: Word8 ::> Word8 ::> Word8 ::> IO CUInt

-- | Apply an IO action on mutable matrices.
apMMatF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> (a -> m b) -> (CUInt -> CUInt -> Ptr Word8 -> a) -> m b
apMMatF2 x g f = do
  let ptr = unsafeForeignPtrToPtr (fst (MVS.unsafeToForeignPtr0 (mmatData x)))
  r <- fromIntegral <$> readMutVar (mrows x)
  c <- fromIntegral <$> readMutVar (mcols x)
  g (f r c ptr)


(#|) :: (PrimMonad m) => m CInt -> String -> m ()
f #| msg = do
  err <- f
  when (err /= 0) $ fail msg

-- | Transpose a mutable matrix.
mcopyTrF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> m (MMatrixF2 (PrimState m))
mcopyTrF2 mMat = do
  r <- readMutVar (mrows mMat)
  c <- readMutVar (mcols mMat)
  mDest <- unsafeNewMatrixF2 c r
  apMMatF2 mMat (apMMatF2 mDest unsafeIOToPrim) c_copyTranspose #| "c_transpose"
  return mDest

-- | Compute a multiplication of two matrices.
mmatMulF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> MMatrixF2 (PrimState m) -> m (MMatrixF2 (PrimState m))
mmatMulF2 mmatA mmatB = do
  r <- readMutVar (mrows mmatA)
  c <- readMutVar (mcols mmatB)
  mCr <- newMutVar r
  mCc <- newMutVar c
  mmatC <- MMatrixF2 mCr mCc <$> MVS.unsafeNew (r*c)
  apMMatF2 mmatA (apMMatF2 mmatB (apMMatF2 mmatC unsafeIOToPrim)) c_matmul #| "c_matmul"
  return mmatC

-- | Row elimination
melimRowsF2 :: (Num a, PrimMonad m) => MMatrixF2 (PrimState m) -> MMatrixF2 (PrimState m) -> MMatrixF2 (PrimState m) -> m a
melimRowsF2 mmatU mmatUinv mmatA =
  fromIntegral <$> apMMatF2 mmatU (apMMatF2 mmatUinv (apMMatF2 mmatA unsafeIOToPrim)) c_elimRows

-- | Diagonal form
mdiagRepF2 :: (Num a, PrimMonad m) => MMatrixF2 (PrimState m) -> MMatrixF2 (PrimState m) -> MMatrixF2 (PrimState m) -> m a
mdiagRepF2 mmatU mmatA mmatV =
  fromIntegral <$> apMMatF2 mmatU (apMMatF2 mmatA (apMMatF2 mmatV unsafeIOToPrim)) c_diagRep
