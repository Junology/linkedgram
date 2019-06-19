------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntMatrix.SmithNF
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Compute Smith normal form
--
------------------------------------------------

module Numeric.Algebra.IntMatrix.SmithNF where

import Control.Monad (when, forM_, guard)
import Control.Monad.ST
import Control.Monad.Loops (whileM_, whileJust_)

import Data.List as L
import Data.Maybe
import Data.STRef

import Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra.Devel

import Numeric.Algebra.IntMatrix.HNFLLL (hnfLLLST)

import Debug.Trace

-- | Find off-diagonal entries and calculate min of min{row,col} among them
findDiagSz :: LA.Matrix LA.Z -> Maybe Int
findDiagSz ma =
  let offdiag = filter (uncurry (/=)) $ LA.find (/=0) ma
  in L.foldl' bin Nothing offdiag
  where
    bin Nothing (i,j) = Just (min i j)
    bin (Just k) (i,j) = Just $ L.foldl' min k [i,j]

-- | Compute "pre-Smith" form with a function calculating Hermite normal forms
-- "Pre-Smith" means diag(c_1,c_2,..) with not necessarily c_i|c_{i+1}
preSmithNFST :: STMatrix s LA.Z -> STMatrix s LA.Z ->  STMatrix s LA.Z -> ST s ()
preSmithNFST stMatUL stMatA stMatURt = do
  whileJust_ (liftSTMatrix findDiagSz stMatA) $ \diagSz -> do
    -- Extract a submatrix containing all the off-diagonal entries.
    stMatAex <- thawMatrix =<< extractMatrix stMatA (FromRow diagSz) (FromCol diagSz)
    stMatULex <- thawMatrix =<< extractMatrix stMatUL (FromRow diagSz) AllCols
    stMatURtex <- thawMatrix =<< extractMatrix stMatURt (FromRow diagSz) AllCols
    -- Compute the Hermite normal forms of the extracted submatrix and of the transpose of the result.
    hnfLLLST stMatULex stMatAex
    stMatAext <- thawMatrix =<< (LA.tr' <$> freezeMatrix stMatAex)
    hnfLLLST stMatURtex stMatAext
    -- Write out the result
    setMatrix stMatA diagSz diagSz =<< (LA.tr' <$> freezeMatrix stMatAext)
    setMatrix stMatUL diagSz 0 =<< freezeMatrix stMatULex
    setMatrix stMatURt diagSz 0 =<< freezeMatrix stMatURtex

isHeadGCD :: LA.Vector LA.Z -> Bool
isHeadGCD = isHeadGCD' . LA.toList
  where
    isHeadGCD' [] = True
    isHeadGCD' (x:xs) = L.foldl' (\b y -> b && rem y x == 0) True xs

extractMaxNZDiag :: LA.Matrix LA.Z -> (Int,LA.Matrix LA.Z)
extractMaxNZDiag mx = runST $ do
  dRef <- newSTRef 0
  let p d = d < uncurry min (LA.size mx) && (mx!d!d) /= 0
  whileM_ (p <$> readSTRef dRef) $ modifySTRef' dRef (+1)
  d <- readSTRef dRef
  return (d,LA.subMatrix (0,0) (d,d) mx)

-- | Normalize diagonals to form a divisible chain; i.e. c_1 | c_2 | c_3 | ...
normalizeST :: STMatrix s LA.Z -> STMatrix s LA.Z ->  STMatrix s LA.Z -> ST s ()
normalizeST stMatUL stMatA stMatURt = do
  (dlen,matD) <- liftSTMatrix extractMaxNZDiag stMatA
  stMatD <- thawMatrix matD
  stMatULex <- thawMatrix =<< extractMatrix stMatUL (RowRange 0 (dlen-1)) AllCols
  stMatURtex <- thawMatrix =<< extractMatrix stMatURt (RowRange 0 (dlen-1)) AllCols
  normalizeST' stMatULex stMatD stMatURtex
  setMatrix stMatA 0 0 =<< freezeMatrix stMatD
  setMatrix stMatUL 0 0 =<< freezeMatrix stMatULex
  setMatrix stMatURt 0 0 =<< freezeMatrix stMatURtex
  where
    normalizeST' stMatULex stMatD stMatURtex = do
      b <- (isHeadGCD . LA.takeDiag) <$> freezeMatrix stMatD
      n <- liftSTMatrix LA.rows stMatD
      when (not b) $ do
        forM_ [1..(n-1)] $ \i -> do
          readMatrix stMatD i i >>= writeMatrix stMatD i 0
          rowOper (AXPY 1 i 0 AllCols) stMatURtex
        preSmithNFST stMatULex stMatD stMatURtex
      when (n > 1) $ do
        stMatDex <- thawMatrix =<< extractMatrix stMatD (FromRow 1) (FromCol 1)
        stMatULexex <- thawMatrix =<< extractMatrix stMatULex (FromRow 1) AllCols
        stMatURtexex <- thawMatrix =<< extractMatrix stMatURtex (FromRow 1) AllCols
        normalizeST' stMatULexex stMatDex stMatURtexex
        setMatrix stMatD 1 1 =<< freezeMatrix stMatDex
        setMatrix stMatULex 1 0 =<< freezeMatrix stMatULexex
        setMatrix stMatURtex 1 0 =<< freezeMatrix stMatURtexex

-- | Compute Smith normal form together with transform unimodular matrices.
smithNF :: LA.Matrix LA.Z -> (LA.Matrix LA.Z, LA.Matrix LA.Z, LA.Matrix LA.Z)
smithNF matA = runST $ do
  stMatA <- thawMatrix matA
  stMatUL <- thawMatrix $ LA.ident (LA.rows matA)
  stMatURt <- thawMatrix $ LA.ident (LA.cols matA)
  preSmithNFST stMatUL stMatA stMatURt
  normalizeST stMatUL stMatA stMatURt
  resMatS <- freezeMatrix stMatA
  resMatUL <- freezeMatrix stMatUL
  resMatURt <- freezeMatrix stMatURt
  return (resMatUL, resMatS, LA.tr' resMatURt)

-- | For debug.
-- Compute "pre-Smith form."
preSmithNF :: LA.Matrix LA.Z -> (LA.Matrix LA.Z, LA.Matrix LA.Z, LA.Matrix LA.Z)
preSmithNF matA = runST $ do
  stMatA <- thawMatrix matA
  stMatUL <- thawMatrix $ LA.ident (LA.rows matA)
  stMatURt <- thawMatrix $ LA.ident (LA.cols matA)
  preSmithNFST stMatUL stMatA stMatURt
  resMatS <- freezeMatrix stMatA
  resMatUL <- freezeMatrix stMatUL
  resMatURt <- freezeMatrix stMatURt
  return (resMatUL, resMatS, LA.tr' resMatURt)
