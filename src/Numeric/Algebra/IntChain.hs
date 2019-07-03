{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntChain
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Chain complexes over integers
--
------------------------------------------------

module Numeric.Algebra.IntChain where

import qualified Numeric.LinearAlgebra as LA

import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV

import Numeric.Algebra.FreeModule
import Numeric.Algebra.IntMatrix

-- |
-- The "lengths" must be consistent; i.e.
-- "length bases" == "length diffs" + 1
data IntChain t a = IntChain {
  bases :: t [a],
  diffs :: t (a -> FreeMod Int a) }

data Homology t a = Homology {
  cycles :: t [FreeMod Int a],
  bndrys :: t [FreeMod Int a],
  ranks :: t Int,
  torsions :: t [Int] }

--------------------------------------
-- Chain level Gaussian elimination --
--------------------------------------
-- | Execute the following algorithm on the given matrix, say A:
-- * Compute the Hermite normal form H=UA^t for the transpose;
-- * Extract the vectors spanning the kernel of A from rows of U;
-- * Put H' the matrix obtained from H by omitting rows which have no non-zero entry;
-- * Put H'' the matrix obtained from H' by omitting rows which have 1 on the bottom columns;
-- * return (kernel, index of (colvec H' \\ colvec H''), transpose of H'')
squeezeMat :: LA.Matrix LA.Z -> ([LA.Vector LA.Z], [Int], LA.Matrix LA.Z)
squeezeMat matA = undefined matA -- TO BE DEFINED
