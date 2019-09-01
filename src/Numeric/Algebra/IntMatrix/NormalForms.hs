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

import Control.Monad
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
foreign import ccall unsafe "c_hermiteNF_LLL" cHermiteNF :: Z ::> Z ::> Z ::> Ok

hermiteNFST :: STMatrix s Z -> STMatrix s Z -> STMatrix s Z -> ST s ()
hermiteNFST stU stUi stH = do
  u <- unsafeFreezeMatrix stU
  ui <- unsafeFreezeMatrix stUi
  h <- unsafeFreezeMatrix stH
  unsafeIOToST (apply u (apply ui (apply h id)) cHermiteNF #| "cHermiteNF")

hermiteNFDecomp :: Matrix Z -> (Matrix Z, Matrix Z, Matrix Z)
hermiteNFDecomp mat = runST $ do
  stMat <- thawMatrix mat
  stU <- thawMatrix $ ident (rows mat)
  stUi <- thawMatrix $ ident (rows mat)
  hermiteNFST stU stUi stMat
  (\x y z -> (x,y,z)) <$!> unsafeFreezeMatrix stU <*> unsafeFreezeMatrix stUi <*> unsafeFreezeMatrix stMat

hermiteNF :: Matrix Z -> (Matrix Z, Matrix Z)
hermiteNF = p . hermiteNFDecomp
  where
    p (x,y,z) = (y,z)

hermiteDecomp :: Matrix Z -> (Matrix Z, Matrix Z)
hermiteDecomp = p . hermiteNFDecomp
  where
    p (x,y,z) = (x,z)

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
