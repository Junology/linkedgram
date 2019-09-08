{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
{-# LANGUAGE TypeFamilies, TypeFamilyDependencies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE DefaultSignatures #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.Presentation
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Matrix presentations of morphisms
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
import Numeric.Matrix.Integral
import Numeric.F2 (F2)
import qualified Numeric.Matrix.F2 as F2
import qualified Numeric.Matrix.F2.Mutable as F2


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
  vecToList :: Vec a -> [a]
  -- ^ Vec a is convertible to a list.
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
  scale :: a -> Vec a -> Vec a
  -- ^ Scalar multiplication on vectors.
  scalar :: Int -> a -> Mat a
  -- ^ Scalar * identity matrix.
  -- | Default implementations
  ident n = scalar n 1
  rankVec = length . vecToList
  (@|>) x = head . toVecs . (x<>) . fromVecs . (:[])
  scale x v = scalar (rankVec v) x @|> v
  scalar n x = fromVecs (map (scale x) (toVecs (ident n)))

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


------------------------------------------------
-- * Instances in terms of LA.Z using hmatrix
------------------------------------------------
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
  vecToList = LA.toList
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
      LA.writeMatrix stMat i i  (invs `LA.atIndex` i)
    LA.unsafeFreezeMatrix stMat
  basisDom = _basisS
  basisCod = _basisT
  toNF x
    = let (u,s,v) = smithRep x
      in SmithNF (LA.takeDiag s) v u
  invariants = LA.toList . _invFactor

---------------------------------
-- * Instances in terms of F2
---------------------------------
instance Coefficient F2 where
  type Vec F2 = F2.VectorF2
  type Mat F2 = F2.MatrixF2
  ident = F2.ident
  present f dom cod
    = let domRk = length dom
          codRk = length cod
      in runST $ do
    stMat <- F2.unsafeNewMatrixF2 codRk domRk
    forM_ [0..domRk-1] $ \j -> do
      let !vals = force $ f (dom !!j)
      forM_ [0..codRk-1] $ \i -> do
        let coeff = getCoeff (cod !! i) vals
        F2.unsafeWriteMatrixF2 stMat i j (fromIntegral coeff :: F2)
    F2.unsafeFreezeMatrixF2 stMat
  vecToList = F2.toList
  rankVec = F2.sizeF2
  rankDom = F2.cols
  rankCod = F2.rows
  fromVecs = F2.fromColumns . toList
  toVecs = F2.toColumns
  scalar = F2.scalar

instance Normalized F2 where
  data NormalForm F2 =
    DiagNF {-# UNPACK #-} !Int {-# UNPACK #-} !F2.MatrixF2 {-# UNPACK #-} !F2.MatrixF2
    deriving (Show, Generic,NFData)
  invariants (DiagNF rk _ _) = replicate rk 1
  basisDom (DiagNF _ bss _) = bss
  basisCod (DiagNF _ _ bst) = bst
  fromNF (DiagNF rk bss bst) = F2.diagRectL (rankCod bst) (rankCod bss) (repeat 1 :: [Int])
  toNF x
    = let (u,rk,v) = F2.diagRep x
      in DiagNF rk v u
