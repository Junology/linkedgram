{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntMatrix.HermiteNFLLL
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Compute Hermite normal forms in the LLL-based method.
-- The original "pseudo-code" is found in the paper
--   George Havas, Bohdan S. Majewski & Keith R. Matthews (1998) Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction, Experimental Mathematics, 7:2, 125-136, DOI: 10.1080/10586458.1998.10504362 
--
------------------------------------------------

module Numeric.Algebra.IntMatrix.HermiteNFLLL (
  hermiteNFST,
  hermiteNF
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

foreign import ccall unsafe "c_hermiteNFUD_LLL" cHermiteNF :: Z ::> Z ::> Ok

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
