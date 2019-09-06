{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE TypeFamilies, TypeFamilyDependencies #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.Homology
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Chain complexes over integers
--
------------------------------------------------

module Numeric.Algebra.Homology
  (Homology(..), homology) where

import GHC.Generics (Generic, Generic1)
import Control.DeepSeq (NFData, force)

import Control.Parallel.Strategies

import Control.Monad
import Control.Monad.ST

import Control.Applicative

import Data.STRef

import Data.Bifunctor
import Data.Traversable
import Data.Tuple (swap)
import Data.Maybe
import Data.List

import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Devel as LA

import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Presentation
import Numeric.Algebra.IntMatrix
import Numeric.Algebra.IntMatrix.NormalForms

{-- for Debug
import Debug.Trace

traceCond :: Bool -> String -> a -> a
traceCond False _  = id
traceCond True msg = trace msg
--}

--------------
-- * Classes
--------------
class (Normalized a) => ChainEliminable a where
  elimChainHead :: Mat a -> Mat a -> (NormalForm a, Mat a)
  -- ^ If
  -- @ let (a',v) = elimChainHead u a @
  -- Then
  -- @ a' == toNF (a <> u^{-1}) @
  -- and the vectors of v (i.e. toVecs v) induces a basis for the cokernel of a⊗ℚ.


-------------------------------------
-- * Computation of homology groups
-------------------------------------
data Homology a = Homology {
  freeCycs :: [Vec a],
  torsions :: [(a,Vec a)],
  bndries :: [Vec a]
  } deriving (Generic)

deriving instance (Coefficient a, Show a, Show (Vec a)) => Show (Homology a)
deriving instance (Coefficient a) => NFData (Homology a)

-- | Compute the kernel and the image of a matrix from its normal form.
homologyOfNF :: (Normalized a) => NormalForm a -> (Homology a, Homology a)
homologyOfNF nf
  = let invs = invariants nf ++ repeat 0
        cycs = map snd $ filter ((==0) . fst) $ zip invs (toVecs (basisDom nf))
        img = filter ((/=0) . fst) $ zip invs (toVecs (basisCod nf))
        (tors,bnds') = partition ((/=1) . fst) img
        bnds = map (uncurry scale) bnds'
    in (Homology [] tors bnds, Homology cycs [] [])

-- | Decompose complexes into two-term complexes with differential in a Smith normal form.
chainToNormals :: (ChainEliminable a, Traversable t) => t (Mat a) -> (Mat a,t (NormalForm a))
chainToNormals = first (maybe (ident 0) id) . mapAccumL chainToNormals' Nothing
  where
    chainToNormals' mayU mat
      = let (u,mat') = case mayU of
                         Just u -> (u, mat <> u)
                         Nothing -> (ident (rankDom mat), mat)
            (matNF, coker) = elimChainHead u mat'
        in (Just coker, matNF)

-- | Compute homology group with integral coefficient.
homology :: (ChainEliminable a, Traversable t, Alternative t, NFData (t (NormalForm a)), NFData (t (Homology a, Homology a))) => t (Mat a) -> t (Homology a)
homology diffs =
  let (!coker, !nfs) = force (chainToNormals diffs)
      !hs = force (fmap homologyOfNF nfs `using` parTraversable (rparWith rdeepseq))
  in (uncurry ((<|>) . pure)) (mapAccumR pile (Homology (toVecs coker) [] []) hs)
  where
    pile (Homology cycs _ _) (thisH, nextCycH) =
      (nextCycH, thisH {freeCycs = cycs})

intHomology :: [LA.Matrix LA.Z] -> [Homology LA.Z]
intHomology = homology

-----------------
-- * Instances
-----------------
-- | Enables chain-level Gaussian eliminations on integer matrices.
instance ChainEliminable LA.Z where
  elimChainHead basis matA = runST $ do
    stMatA <- LA.thawMatrix matA
    stBasis <- LA.thawMatrix basis
    stU <- LA.thawMatrix (LA.ident (LA.rows matA))
    -- Compute the Smith normal form
    smithRepST stU stMatA stBasis
    -- Write the new basis of the source
    resBasis <- LA.unsafeFreezeMatrix stBasis
    -- Take the invariant factor.
    invFs <- LA.takeDiag <$!> LA.unsafeFreezeMatrix stMatA
    let rkQ = length (takeWhile (/=0) (LA.toList invFs)) -- rank over the field of fractions; e.g. Q for Z
    -- Free closure of the image of matA & its orthogonal.
    fimage <- if LA.rows matA > 0 && rkQ > 0
              then LA.extractMatrix stU LA.AllRows (LA.ColRange 0 (rkQ-1))
              else return $ (LA.rows matA LA.>< 0) []
    fcoker <- if LA.rows matA > 0 && rkQ < LA.rows matA
              then LA.extractMatrix stU LA.AllRows (LA.FromCol rkQ)
              else return $ (LA.rows matA LA.>< 0) []
    return (SmithNF invFs resBasis fimage, fcoker)

