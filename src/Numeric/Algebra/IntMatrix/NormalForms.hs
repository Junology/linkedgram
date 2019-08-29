{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntMatrix.NormalForms
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Hermite Normal form and Smith normal form.
-- For the algorithms, see c-source files in src/Internal/C.
--
------------------------------------------------

module Numeric.Algebra.IntMatrix.NormalForms (
  hermiteNFST,
  hermiteNF,
  smithNFST,
  smithNF,
  smithRepST,
  smithRep
  ) where

import Control.Monad.ST
import Control.Monad.ST.Unsafe

import Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra.Devel

import System.IO.Unsafe (unsafePerformIO)
import Foreign.C.Types (CLong(..), CInt(..))
import Foreign.Ptr (Ptr)

-------------------------------------------------------------
-- Copied from the source of hmatrix package.
infixr 5 :>, ::>
type (:>)  t r = CInt -> Ptr t -> r
type (::>) t r =  CInt -> CInt -> CInt -> CInt -> Ptr t -> r
type Ok = IO CInt

--------------------------------------------------------------

-- | Hermite normal form
foreign import ccall unsafe "c_hermiteNF_LLL_partial" cHermiteNFP :: CInt -> Z ::> Z ::> Ok
foreign import ccall unsafe "c_hermiteNF_LLL" cHermiteNF :: Z ::> Z ::> Ok

hermiteNFPST :: Int -> STMatrix s Z -> STMatrix s Z -> ST s ()
hermiteNFPST k stU stH = do
  u <- unsafeFreezeMatrix stU
  h <- unsafeFreezeMatrix stH
  unsafeIOToST (apply u (apply h id) (cHermiteNFP (fromIntegral k)) #| "cHermiteNFP")

hermiteNFP :: Int -> Matrix Z -> (Matrix Z, Matrix Z)
hermiteNFP k mat = runST $ do
  stMat <- thawMatrix mat
  stU <- thawMatrix $ ident (rows mat)
  hermiteNFPST k stU stMat
  (,) <$> unsafeFreezeMatrix stU <*> unsafeFreezeMatrix stMat

hermiteNFST :: STMatrix s Z -> STMatrix s Z -> ST s ()
hermiteNFST stU stH = do
  u <- unsafeFreezeMatrix stU
  h <- unsafeFreezeMatrix stH
  unsafeIOToST (apply u (apply h id) cHermiteNF #| "cHermiteNF")

hermiteNF :: Matrix Z -> (Matrix Z, Matrix Z)
hermiteNF mat = runST $ do
  stMat <- thawMatrix mat
  stU <- thawMatrix $ ident (rows mat)
  hermiteNFST stU stMat
  (,) <$> unsafeFreezeMatrix stU <*> unsafeFreezeMatrix stMat

-- | Smith normal form
foreign import ccall unsafe "c_smithNF" cSmithNF :: Z ::> Z ::> Z ::> Ok

smithNFST :: STMatrix s LA.Z -> STMatrix s LA.Z ->  STMatrix s LA.Z -> ST s ()
smithNFST stU stM stV = do
  u <- unsafeFreezeMatrix stU
  m <- unsafeFreezeMatrix stM
  v <- unsafeFreezeMatrix stV
  unsafeIOToST (apply u (apply m (apply v id)) cSmithNF #| "cSmithNF")

smithNF :: Matrix Z -> (Matrix Z, Matrix Z, Matrix Z)
smithNF mat = runST $ do
  stMat <- thawMatrix mat
  stU <- thawMatrix $ ident (rows mat)
  stV <- thawMatrix $ ident (cols mat)
  smithNFST stU stMat stV
  (\x y z -> (x,y,z))
    <$> unsafeFreezeMatrix stU
    <*> unsafeFreezeMatrix stMat
    <*> unsafeFreezeMatrix stV

-- | Smith representation
foreign import ccall unsafe "c_smithRep" cSmithRep :: Z ::> Z ::> Z::> Ok

smithRepST :: STMatrix s LA.Z -> STMatrix s LA.Z -> STMatrix s LA.Z -> ST s ()
smithRepST stA stM stB = do
  a <- unsafeFreezeMatrix stA
  m <- unsafeFreezeMatrix stM
  b <- unsafeFreezeMatrix stB
  unsafeIOToST (apply a (apply m (apply b id)) cSmithRep #| "cSmithRep")

smithRep :: Matrix Z -> (Matrix Z, Matrix Z, Matrix Z)
smithRep mat = runST $ do
  stMat <- thawMatrix mat
  stU <- thawMatrix $ ident (rows mat)
  stV <- thawMatrix $ ident (cols mat)
  smithRepST stU stMat stV
  (\x y z -> (x,y,z))
    <$> unsafeFreezeMatrix stU
    <*> unsafeFreezeMatrix stMat
    <*> unsafeFreezeMatrix stV
