{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
{-# LANGUAGE TypeOperators #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Matrix.F2
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Matrix with entries in F2; the prime field of characteristic 2
--
------------------------------------------------

module Numeric.Matrix.F2 where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData)

import Control.Monad.Primitive
import Data.Primitive.MutVar

import Control.Monad (join, forM_)
import Control.Monad.ST (runST)
import Data.Bits
import Data.Word (Word8)

import Foreign.Marshal (fromBool)

import Data.Bifunctor (bimap)
import qualified Data.List as L

import qualified Data.Vector.Storable  as VS
import qualified Data.Vector.Storable.Mutable as MVS

import Numeric.Matrix.F2.Mutable

------------------
-- * Data types
------------------

-- | Matrix with entries in F2
--   The row-major convention is always used.
data MatrixF2 = MatrixF2 {
  rows :: {-# UNPACK #-} !Int,
  cols :: {-# UNPACK #-} !Int,
  matData :: {-# UNPACK #-} !(VS.Vector Word8)
  } deriving (Eq,Ord,Generic,NFData)

instance Show MatrixF2 where
  show (MatrixF2 r c dt) =
    "(" ++ show r ++ "><" ++ show c ++ ")[\n"
    ++ L.intercalate "\n" (fmap show $ mkChunks c dt)
    ++ "]\n"

-- | Vectors with entries in F2
newtype VectorF2 = VectorF2 { vecData :: VS.Vector Word8 }
  deriving (Eq,Ord,Show,Generic,NFData)

-- | The size of a vector.
sizeF2 :: VectorF2 -> Int
sizeF2 = VS.length . vecData

-- | Divide a vector into a list of chunks of a fixed length.
mkChunks :: (VS.Storable a) => Int -> VS.Vector a -> [VS.Vector a]
mkChunks n v =
  L.unfoldr (\x -> if VS.null x then Nothing else Just (VS.splitAt n x)) v

------------------------------------------
-- * Conversion to/from mutables vectors
------------------------------------------
-- | Immutable copy of a mutable vector.
freezeVectorF2 :: (PrimMonad m) => MVectorF2 (PrimState m) -> m VectorF2
freezeVectorF2 mvec = VectorF2 <$> VS.freeze (mvecData mvec)

-- | Mutable copy of immutable vector.
thawVectorF2 :: (PrimMonad m) => VectorF2 -> m (MVectorF2 (PrimState m))
thawVectorF2 vec = MVectorF2 <$> VS.thaw (vecData vec)

-- | Unsafe conversion of a mutable vector to an immutable one
unsafeFreezeVectorF2 :: (PrimMonad m) => MVectorF2 (PrimState m) -> m VectorF2
unsafeFreezeVectorF2 mvec = VectorF2 <$> VS.unsafeFreeze (mvecData mvec)

-- | Unsafe conversion of an immutable vector to a mutable one.
unsafeThawVectorF2 :: (PrimMonad m) => VectorF2 -> m (MVectorF2 (PrimState m))
unsafeThawVectorF2 vec = MVectorF2 <$> VS.unsafeThaw (vecData vec)


-------------------------------------
-- * Construction of Vector/Matrix
-------------------------------------

-- | Vector from a list.
fromList :: (Bits a) => [a] -> VectorF2
fromList = VectorF2 . VS.fromList . fmap (fromBool . flip testBit 0)

-- | Matrix from a list.
--   The interface is same as that of hmatrix package.
(><) :: (Bits a) => Int -> Int -> [a] -> MatrixF2
(><) r c =
  MatrixF2 r c . VS.fromList . fmap (fromBool . flip testBit 0) . take (r*c)

-- | Create a rectangular matrix with diagonals from a list
diagRectL :: (Bits a) => Int -> Int -> [a] -> MatrixF2
diagRectL r c xs = runST $ do
  stVec <- MVS.replicate (r*c) 0
  forM_ (zip [0..] xs) $ \ix -> do
    let (i,x) = ix
    MVS.write stVec (i*(r+1)) (fromBool (testBit x 0))
  MatrixF2 r c <$> VS.unsafeFreeze stVec

-- | Create a diagonal matrix from a list
diagL :: (Bits a) => [a] -> MatrixF2
diagL xs = let n = length xs
           in diagRectL n n xs

-- | The identity matrix
ident :: Int -> MatrixF2
ident n = diagL (replicate n True)

-- | The scalar matrix
scalar :: (Bits a) => Int -> a -> MatrixF2
scalar n x = diagL (replicate n x)

-- | Transpose
tr :: MatrixF2 -> MatrixF2
tr mat = runST $ do
  stMat <- thawMatrixF2 mat
  unsafeFreezeMatrixF2 =<< mcopyTrF2 stMat


------------------------------------------
-- * Conversion to/from mutables matrices
------------------------------------------
-- | Immutable copy of a mutable matrix.
freezeMatrixF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> m MatrixF2
freezeMatrixF2 mmat = do
  r <- readMutVar (mrows mmat)
  c <- readMutVar (mcols mmat)
  MatrixF2 r c <$> VS.freeze (mmatData mmat)

-- | Mutable copy of immutable matrix.
thawMatrixF2 :: (PrimMonad m) => MatrixF2 -> m (MMatrixF2 (PrimState m))
thawMatrixF2 mat = do
  mr <- newMutVar (rows mat)
  mc <- newMutVar (cols mat)
  MMatrixF2 mr mc <$> VS.thaw (matData mat)

-- | Unsafe conversion of a mutable matrix to an immutable one
unsafeFreezeMatrixF2 :: (PrimMonad m) => MMatrixF2 (PrimState m) -> m MatrixF2
unsafeFreezeMatrixF2 mmat = do
  r <- readMutVar (mrows mmat)
  c <- readMutVar (mcols mmat)
  MatrixF2 r c <$> VS.unsafeFreeze (mmatData mmat)

-- | Unsafe conversion of an immutable matrix to a mutable one.
unsafeThawMatrixF2 :: (PrimMonad m) => MatrixF2 -> m (MMatrixF2 (PrimState m))
unsafeThawMatrixF2 mat = do
  mr <- newMutVar (rows mat)
  mc <- newMutVar (cols mat)
  MMatrixF2 mr mc <$> VS.unsafeThaw (matData mat)


---------------------------------------
-- * Vector to/from Matrix conversion
---------------------------------------

-- | Convert a vector to a matrix with it as a row vector.
asRow :: VectorF2 -> MatrixF2
asRow (VectorF2 dt) = MatrixF2 1 (VS.length dt) dt

-- | Convert a vector to a matrix with it as a column vector.
asColumn :: VectorF2 -> MatrixF2
asColumn (VectorF2 dt) = MatrixF2 (VS.length dt) 1 dt

-- | Create a matrix from row vectors.
fromRows :: [VectorF2] -> MatrixF2
fromRows [] = MatrixF2 0 0 VS.empty
fromRows cols@((VectorF2 v):_) =
  let dts = fmap vecData cols
  in MatrixF2 (L.length cols) (VS.length v) (VS.concat dts)

-- | Extract row vectors.
toRows :: MatrixF2 -> [VectorF2]
toRows (MatrixF2 _ c dt)
  = VectorF2 <$> mkChunks c dt

-- | Create a matrix from row vectors.
fromColumns :: [VectorF2] -> MatrixF2
fromColumns = tr . fromRows

-- | Extract column vectors.
toColumns :: MatrixF2 -> [VectorF2]
toColumns = toRows . tr


----------------------------------------
-- * Extraction of submatrices
----------------------------------------
splitRowsAt :: Int -> MatrixF2 -> (MatrixF2,MatrixF2)
splitRowsAt n (MatrixF2 r c dt) =
  let (dt1,dt2) = VS.splitAt (n*c) dt
  in (MatrixF2 n c dt1, MatrixF2 (r-n) c dt2)

splitColumnsAt :: Int -> MatrixF2 -> (MatrixF2,MatrixF2)
splitColumnsAt n = bimap tr tr . splitRowsAt n . tr

-------------------
-- * Elimination
-------------------
-- | Decompose a matrix A so that
-- @
--   let (q,q',u) = elimRows a
--   in a == q <> u && q <> q' = ident (rows A)
-- @
-- yields True.
elimRows :: MatrixF2 -> (MatrixF2, MatrixF2, MatrixF2)
elimRows mat = runST $ do
  stMat <- thawMatrixF2 mat
  stQ <- newIdentF2 (rows mat)
  stQ' <- newIdentF2 (rows mat)
  melimRowsF2 stQ stQ' stMat
  (,,) <$> unsafeFreezeMatrixF2 stQ
    <*> unsafeFreezeMatrixF2 stQ'
    <*> unsafeFreezeMatrixF2 stMat

-- | Compute a presentation by a diagonal matrix; i.e.
-- @
--   let n = cols a
--       (q,rk,p) = diagRep a
--   in a <> p == q <> (diagL (replicate rk True ++ replicate (n-rk) False))
-- @
-- yields True.
diagRep :: MatrixF2 -> (MatrixF2, Int, MatrixF2)
diagRep mat = runST $ do
  stMat <- thawMatrixF2 mat
  stQ <- newIdentF2 (rows mat)
  stP <- newIdentF2 (cols mat)
  rk <- mdiagRepF2 stQ stMat stP
  (\x y -> (x,rk,y)) <$> unsafeFreezeMatrixF2 stQ <*> unsafeFreezeMatrixF2 stP

------------
-- * Misc
------------
toList :: (Bits a) => VectorF2 -> [a]
toList (VectorF2 vec) = map f (VS.toList vec)
  where
    f x = if testBit x 0
          then bit 0
          else zeroBits

-----------------
-- * Instances
-----------------
instance Semigroup MatrixF2 where
  matA <> matB = runST $ unsafeFreezeMatrixF2 =<<
    join (mmatMulF2 <$> unsafeThawMatrixF2 matA <*> unsafeThawMatrixF2 matB)
