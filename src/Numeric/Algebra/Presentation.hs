{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
{-# LANGUAGE TypeFamilies, TypeFamilyDependencies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE DefaultSignatures #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntChain
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Chain complexes over integers
--
------------------------------------------------

module Numeric.Algebra.Presentation where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData, force)

import Control.Monad
import Control.Monad.ST

import Data.Foldable

import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Devel as LA

import Numeric.Algebra.FreeModule
import Numeric.Algebra.IntMatrix


---------------
-- * Classes
---------------
-- | Class which provides matrix presentations for homomorphisms with coefficients in a.
class (Eq a, Num a, NFData a, Semigroup (Mat a), NFData (Vec a), NFData (Mat a)) => Coefficient a where
  type Vec a = v | v -> a
  -- ^ The type of vectors with coefficients in a.
  type Mat a = m | m -> a
  -- ^ The type of matrices with coefficients in a.
  ident :: Int -> Mat a
  -- ^ Presentation of the identity of the given rank.
  present :: (Integral a', NFData a', Ord c, NFData c) => (b->FreeMod a' c) -> [b] -> [c] -> Mat a
  -- ^ Present a linear map in a matrix form with respect to a given bases
  rankVec :: Vec a -> Int
  -- ^ The rank of the module where the given vector lies.
  rankDom :: Mat a -> Int
  -- ^ The rank of the domain of the given matrix.
  rankCod :: Mat a -> Int
  -- ^ The rank of the codomain of the given matrix.
  fromVecs :: (Foldable t) => t (Vec a) -> Mat a
  -- ^ Give a matrix whose image is spanned by given vectors.
  -- Typically, it is the matrix with given vectors as columns.
  toVecs :: Mat a -> [Vec a]
  -- ^ Give a set of vectors which span the image of the matrix.
  -- Typically, it is the set of column vectors.
  -- ^ The scalar matrix of the given size.
  (@|>) :: Mat a -> Vec a -> Vec a
  -- ^ Multiplication of matrices and vectors.
  (@|>) x = head . toVecs . (x<>) . fromVecs . (:[])
  scale :: a -> Vec a -> Vec a
  -- ^ Scalar multiplication on vectors.

-- | Matrix presentations with normal forms.
class (Coefficient a, NFData (NormalForm a)) => Normalized a where
  data NormalForm a
  -- ^ The type of normal forms.
  fromNF :: NormalForm a -> Mat a
  -- ^ Convert into a matrix form.
  -- >> let xnf = toNF x
  -- Then
  -- >> x <> basisDom xnf == basisCod xnf <> fromNF xnf
  basisDom :: NormalForm a -> Mat a
  -- ^ The basis matrix on the domain.
  -- It should be invertible.
  basisCod :: NormalForm a -> Mat a
  -- ^ The basis matrix on the codomain.
  -- It should be invertible.
  toNF :: Mat a -> NormalForm a
  -- ^ Compute the normal form.
  invariants :: NormalForm a -> [a]
  -- ^ Each normal form should contain the data of "invariant factors;" e.g. those for Smith normal forms of matrices with coefficients in a PID.


---------------------------------------
-- * Instances for LA.Z using hmatrix
---------------------------------------
instance Coefficient (LA.Z) where
  type Vec LA.Z = LA.Vector LA.Z
  type Mat LA.Z = LA.Matrix LA.Z
  ident = LA.ident
  present f dom cod
    = let domRk = length dom
          codRk = length cod
      in runST $ do
    stMat <- LA.newUndefinedMatrix LA.RowMajor codRk domRk
    forM_ [0..domRk-1] $ \j -> do
      let !vals = force $ f (dom !!j)
      forM_ [0..codRk-1] $ \i -> do
        let coeff = getCoeff (cod !! i) vals
        LA.unsafeWriteMatrix stMat i j (fromIntegral coeff)
    LA.unsafeFreezeMatrix stMat
  rankVec = LA.size
  rankDom = LA.cols
  rankCod = LA.rows
  fromVecs = LA.fromColumns . toList
  toVecs = LA.toColumns
  (@|>) = (LA.#>)
  scale = LA.scale

instance Normalized (LA.Z) where
  -- | Presentation in Smith normal forms.
  data NormalForm LA.Z = SmithNF {
    _invFactor :: LA.Vector LA.Z, -- Invariant factor
    _basisS :: LA.Matrix LA.Z,    -- basis for source
    _basisT :: LA.Matrix LA.Z     -- basis for target
    } deriving (Show, Generic, NFData)
  fromNF (SmithNF invs bss bts)
    = runST $ do
    stMat <- LA.newMatrix 0 (LA.cols bts) (LA.cols bss)
    forM_ [0..(LA.size invs)-1] $ \i ->
      LA.writeMatrix stMat i i  1
    LA.unsafeFreezeMatrix stMat
  basisDom = _basisS
  basisCod = _basisT
  toNF x
    = let (u,s,v) = smithRep x
      in SmithNF (LA.takeDiag s) v u
  invariants = LA.toList . _invFactor
