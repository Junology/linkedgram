{-# LANGUAGE FlexibleContexts #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntMatrix.HNFLLL
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Compute Hermite normal forms in the LLL-based method.
-- The original "pseudo-code" is found in the paper
--   George Havas, Bohdan S. Majewski & Keith R. Matthews (1998) Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction, Experimental Mathematics, 7:2, 125-136, DOI: 10.1080/10586458.1998.10504362 
--
------------------------------------------------

module Numeric.Algebra.IntMatrix.HNFLLL (
  hnfLLLST,
  hnfLLL,
  ) where

import Control.Monad (when, forM_)
import Control.Monad.ST
import Control.Monad.Loops (whileM_)

import Data.STRef

import qualified Data.Maybe as M

import Data.Int

import Numeric.LinearAlgebra ((!))
import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra.Devel

-- | Type class which guarantees the type can be embedded in the real number field
class (Num t,Eq t,Ord t) => Realizable t where
  toReal :: t -> Double

instance Realizable Double where
  toReal = id

instance Realizable Int64 where
  toReal = fromIntegral

-- | Type class to qualify a type to be used in LLL-based algorithm
class (Realizable t, LA.Indexable (LA.Vector t) t, LA.Container LA.Vector t) => LLLQualified t

instance LLLQualified Double

instance LLLQualified Int64

-- | Data Type to carry resources in the computation of LLL-based algorithm
data HNFLLLResource s t = HNFLLL {
  matA :: STMatrix s t,
  matU :: STMatrix s t,
  lambda :: STMatrix s Double,
  focusC1 :: STRef s (Maybe Int),
  focusC2 :: STRef s (Maybe Int) }

-- | The function which is called Reduce2 in the original "pseudo-code."
reduceLLL :: (LLLQualified t) => Int -> Int -> HNFLLLResource s t -> ST s ()
reduceLLL k i resRef = do
  ma <- freezeMatrix (matA resRef)
  case M.listToMaybe (LA.find (/=0) (ma!i)) of
    (Just j) -> do
      writeSTRef (focusC1 resRef) (Just j)
      when ( (ma!i!j) < 0 ) $ do
        rowOper (SCAL (-1) (Row i) AllCols) (matA resRef)
        rowOper (SCAL (-1) (Row i) AllCols) (matU resRef)
        rowOper (SCAL (-1) (Row (j+1)) (ColRange 0 j)) (lambda resRef)
        rowOper (SCAL (-1) (FromRow (j+2)) (Col (j+1))) (lambda resRef)
    Nothing -> writeSTRef (focusC1 resRef) Nothing
  case M.listToMaybe (LA.find (/=0) (ma!k)) of
    (Just j) -> writeSTRef (focusC2 resRef) (Just j)
    Nothing -> writeSTRef (focusC2 resRef) Nothing
  c1 <- readSTRef (focusC1 resRef)
  q <- fromIntegral.floor <$> do
    case c1 of
      (Just j) ->
        (/) <$> (toReal <$> readMatrix (matA resRef) k j) <*> (toReal <$> readMatrix (matA resRef) i j)
      Nothing -> do
        lki <- readMatrix (lambda resRef) (k+1) (i+1)
        di <- readMatrix (lambda resRef) (i+1) (i+1)
        if 2 * abs lki > di
          then return $ fromIntegral (round (lki / di))
          else return 0
  when (q /= 0) $ do
    rowOper (AXPY (negate q) i k AllCols) (matA resRef)
    rowOper (AXPY (negate q) i k AllCols) (matU resRef)
    rowOper (AXPY (toReal $ negate q) (i+1) (k+1) AllCols) (lambda resRef)

-- | The function which is called Swap1 in the original "pseudo-code."
swapLLL :: (Num t, LA.Element t) => Int -> HNFLLLResource s t -> ST s ()
swapLLL k resRef = do
  rowOper (SWAP k (k-1) AllCols) (matA resRef)
  rowOper (SWAP k (k-1) AllCols) (matU resRef)
  rowOper (SWAP k (k+1) (ColRange 0 (k-1))) (lambda resRef)
  ml <- freezeMatrix (lambda resRef)
  forM_ [(k+2)..(LA.cols ml-1)] $ \i -> do
    writeMatrix (lambda resRef) i k $
      ((ml!i!k)*(ml!(k+1)!k) + (ml!i!(k+1))*(ml!(k-1)!(k-1)))/(ml!k!k)
    writeMatrix (lambda resRef) i (k+1) $
      ((ml!i!k)*(ml!(k+1)!(k+1)) + (ml!i!(k+1))*(ml!(k+1)!k))/(ml!k!k)
  writeMatrix (lambda resRef) k k $
    ((ml!(k+1)!k)*(ml!(k+1)!k) + (ml!(k+1)!(k+1))*(ml!(k-1)!(k-1)))/(ml!k!k)

-- | Recognize whether a swap should occur or not.
shouldSwapST :: (Num t, LA.Element t) => Int -> HNFLLLResource s t -> ST s Bool
shouldSwapST k resRef = do
  ml <- freezeMatrix $ lambda resRef
  mc1 <- readSTRef $ focusC1 resRef
  mc2 <- readSTRef $ focusC2 resRef
  case (mc1,mc2) of
    (Just c1,Just c2) -> return $ c1 <= c2
    (Just _,Nothing) -> return True
    (Nothing,Just _) -> return False
    (Nothing,Nothing) ->
      return $ (ml!(k+1)!k)*(ml!(k+1)!k) + (ml!(k+1)!(k+1))*(ml!(k-1)!(k-1)) < (ml!k!k)*(ml!k!k)

-- | ST to compute the Hermite normal form for the matrix in stMatA.
-- The first argument is also subject to the operation that is executed on stMatA, so one gets the transformation unimodular matrix by passing the identity.
-- See the source of hnfLLL function for an example of usage.
hnfLLLST :: (LLLQualified t) => STMatrix s t -> STMatrix s t -> ST s ()
hnfLLLST stMatU stMatA = do
  m <- LA.rows <$> freezeMatrix stMatA
  newLambda <- thawMatrix $ LA.ident (m+1)
  newC1 <- newSTRef Nothing
  newC2 <- newSTRef Nothing
  let resRef = HNFLLL stMatA stMatU newLambda newC1 newC2
  kRef <- newSTRef (1 :: Int)
  whileM_ ((<m) <$> readSTRef kRef) $ do
    k <- readSTRef kRef
    reduceLLL k (k-1) resRef
    shouldSwap <- shouldSwapST k resRef
    if shouldSwap
      then swapLLL k resRef >> when (k>1) (writeSTRef kRef (k-1))
      else (forM_ [0..k-2] $ \i -> reduceLLL k i resRef) >> writeSTRef kRef (k+1)
  setMatrix stMatA 0 0 =<< (LA.flipud <$> freezeMatrix stMatA)
  setMatrix stMatU 0 0 =<< (LA.flipud <$> freezeMatrix stMatU)

-- | Execute the LLL-based algorithm on a given matrix.
hnfLLL :: LLLQualified t => LA.Matrix t -> (LA.Matrix t, LA.Matrix t)
hnfLLL mtx = runST $ do
  stMatA <- thawMatrix mtx
  stMatU <- thawMatrix $ LA.ident (LA.rows mtx)
  hnfLLLST stMatU stMatA
  h <- freezeMatrix stMatA
  u <- freezeMatrix stMatU
  return (u,h)

