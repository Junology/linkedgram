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
  smithNF
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

foreign import ccall unsafe "c_hermiteNF_LLL" cHermiteNF :: Z ::> Z ::> Ok

hermiteNFST :: STMatrix s Z -> STMatrix s Z -> ST s ()
hermiteNFST stU stH = do
  u <- unsafeFreezeMatrix stU
  h <- unsafeFreezeMatrix stH
  unsafeIOToST (apply u (apply h id) cHermiteNF #| "cHermiteNF")

hermiteNF :: Matrix Z -> (Matrix Z, Matrix Z)
hermiteNF mat = runST $ do
  let (r,c) = size mat
  stMat <- thawMatrix mat
  stU <- thawMatrix $ ident r
  hermiteNFST stU stMat
  (,) <$> unsafeFreezeMatrix stU <*> unsafeFreezeMatrix stMat

-- | Smith normal form
foreign import ccall unsafe "c_smithNF" cSmithNF :: Z ::> Z ::> Z ::> Ok

smithNFST :: STMatrix s LA.Z -> STMatrix s LA.Z ->  STMatrix s LA.Z -> ST s ()
smithNFST stU stM stV = do
  u <- unsafeFreezeMatrix stU
  m <- unsafeFreezeMatrix stM
  v <- unsafeFreezeMatrix stV
  unsafeIOToST $ do
    apply u (apply m (apply v id)) cSmithNF #| "cSmithNF"

smithNF :: Matrix Z -> (Matrix Z, Matrix Z, Matrix Z)
smithNF mat = runST $ do
  let (r,c) = size mat
  stMat <- thawMatrix mat
  stU <- thawMatrix $ ident r
  stV <- thawMatrix $ ident c
  smithNFST stU stMat stV
  (\x y z -> (x,y,z))
    <$> unsafeFreezeMatrix stU
    <*> unsafeFreezeMatrix stMat
    <*> unsafeFreezeMatrix stV
