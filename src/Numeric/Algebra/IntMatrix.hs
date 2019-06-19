------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntMatrix
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Operations on integer matrices
--
------------------------------------------------

module Numeric.Algebra.IntMatrix (
  -- Hermite normal form
  hnfLLL,
  -- Smith normal form
  smithNF,
  preSmithNF, -- for Debug
  -- Translation
  matiDataToLA,
  matdDataToLA,
  matiLAToData,
  matdLAToData
  ) where

import Numeric.Algebra.IntMatrix.HNFLLL
import Numeric.Algebra.IntMatrix.SmithNF

import qualified Numeric.LinearAlgebra as LA

import qualified Data.Matrix as Mat

matiDataToLA :: Integral a => Mat.Matrix a -> LA.Matrix LA.Z
matiDataToLA mt =
  (LA.><) (Mat.nrows mt) (Mat.ncols mt) (map fromIntegral $ Mat.toList mt)

matdDataToLA :: Mat.Matrix Double -> LA.Matrix Double
matdDataToLA mt =
  (LA.><) (Mat.nrows mt) (Mat.ncols mt) (Mat.toList mt)

matiLAToData :: Num a => LA.Matrix LA.Z -> Mat.Matrix a
matiLAToData mt =
  Mat.fromList (LA.rows mt) (LA.cols mt) $ map fromIntegral $ LA.toList $ LA.flatten mt

matdLAToData :: LA.Matrix Double -> Mat.Matrix Double
matdLAToData mt =
  Mat.fromList (LA.rows mt) (LA.cols mt) $ LA.toList $ LA.flatten mt
