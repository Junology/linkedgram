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
  hermiteNF,
  -- Smith normal form
  smithNF,
  -- Compute Image and Kernel
  kerImOf,
  -- Translation
  vectiLAToData,
  matiDataToLA,
  matdDataToLA,
  matiLAToData,
  matdLAToData
  ) where

import Control.Monad
import Control.Monad.ST (ST, runST)
import Control.Monad.Loops (whileM_)
import Data.STRef

--import Numeric.Algebra.IntMatrix.HNFLLL
--import Numeric.Algebra.IntMatrix.SmithNF
import Numeric.Algebra.IntMatrix.NormalForms

import qualified Numeric.LinearAlgebra as LA

import qualified Data.Vector as V
import qualified Data.Matrix as Mat

-- Compute the submatrix containing all the non-zero diagonal entries
extractMaxNZDiag :: LA.Matrix LA.Z -> (Int,LA.Matrix LA.Z)
extractMaxNZDiag mx = runST $ do
  dRef <- newSTRef 0
  let p d = d < uncurry min (LA.size mx) && (mx LA.! d LA.! d) /= 0
  whileM_ (p <$> readSTRef dRef) $ modifySTRef' dRef (+1)
  d <- readSTRef dRef
  return (d,LA.subMatrix (0,0) (d,d) mx)

-- D = P <> A <> Q
-- | Return (d,ker,im) where
-- d: the diagonal entries of Smith normal form
-- ker: basis for the kernel
-- im: basis for the image
kerImOf :: LA.Matrix LA.Z -> (LA.Vector LA.Z,[LA.Vector LA.Z],[LA.Vector LA.Z])
kerImOf mt =
  let (ul,h,ur) = smithNF mt
      (ull,_,ulr) = smithNF ul
      (rk,dmt) = extractMaxNZDiag h
      ker' = ur LA.¿ [rk..(LA.cols ur-1)]
      (_,kert) = hermiteNF (LA.tr' ker')
      im' = (ulr LA.<> ull LA.<> h) LA.¿ [0..rk-1]
      (_,imt) = hermiteNF (LA.tr' im')
  in (LA.takeDiag dmt, LA.toRows kert, LA.toRows imt)

-------------------------
-- Translation of data --
-------------------------
vectiLAToData :: (Num a) => LA.Vector LA.Z -> V.Vector a
vectiLAToData = V.fromList . map fromIntegral . LA.toList

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
