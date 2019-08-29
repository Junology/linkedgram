{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DeriveAnyClass #-}

------------------------------------------------
-- |
-- Module    :  Numeric.Algebra.IntChain
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Chain complexes over integers
--
------------------------------------------------

module Numeric.Algebra.IntChain
  (IntHomology(..), intHomology) where

import GHC.Generics (Generic, Generic1)
import Control.DeepSeq (NFData, NFData1, force)

import Control.Monad
import Control.Monad.ST

import Control.Applicative

import Data.STRef

import Data.Traversable
import Data.Maybe
import Data.List

import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Devel as LA

import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV

import Numeric.Algebra.FreeModule
import Numeric.Algebra.IntMatrix
import Numeric.Algebra.IntMatrix.NormalForms

{-- for Debug
import Debug.Trace
--}

-- | Matrix representation in a Smith normal form.
data RepSmith a = RepSmith {
  invFactor :: [a],      -- Invariant factor
  basisS :: LA.Matrix a, -- basis for source
  basisT :: LA.Matrix a  -- basis for target
  }
  deriving (Show, Generic, Generic1, NFData)

data IntHomology = IntHomology {
  freeCycs :: [LA.Vector LA.Z],
  torCycs :: [(LA.Vector LA.Z, LA.Z)],
  bndries :: [LA.Vector LA.Z] }
  deriving (Show, Generic, NFData)

-- | @ squeezeMat basis a b = (RepSmith u s v, q, b') @ implies
-- * @ a == v LA.<> LA.diag (LA.fromList s) LA.<> u @ is the Smith normal form of a;
-- * b' is the restriction of b to the "cokernel" of a, for which q is the basis; where "cokernel" is a maximal free submodule orthogonal to the image.
squeezeMat :: LA.Matrix LA.Z -> LA.Matrix LA.Z -> LA.Matrix LA.Z -> (RepSmith LA.Z, LA.Matrix LA.Z, LA.Matrix LA.Z)
squeezeMat basis matA matB = runST $ do
  stMatA <- LA.thawMatrix matA
  stBasis <- LA.thawMatrix basis
  stU <- LA.thawMatrix (LA.ident (LA.rows matA))
  -- Compute the Smith normal form
  smithRepST stU stMatA stBasis
  -- Write the new basis of the source
  resBasis <- LA.unsafeFreezeMatrix stBasis
  -- Take the invariant factor.
  invFs <- takeWhile (/=0) <$!> LA.toList <$!> LA.takeDiag <$!> LA.unsafeFreezeMatrix stMatA
  let rkQ = length invFs -- rank over the field of fractions; e.g. Q for Z
  -- Free closure of the image of matA & its orthogonal.
  fimage <- if LA.rows matA > 0
            then LA.extractMatrix stU LA.AllRows (LA.ColRange 0 (rkQ-1))
            else return $ (LA.rows matA LA.>< 0) []
  fcoker <- if LA.rows matA > 0 && rkQ < LA.rows matA
            then LA.extractMatrix stU LA.AllRows (LA.FromCol rkQ)
            else return $ (LA.rows matA LA.>< 0) []
  {-- Debug
  when (not (null invFs)) $ do
    trace (show matA) $ return ()
    trace (show invFs) $ return ()
    trace (show fimage ++ show fcoker) $ return ()
  --}
  return (RepSmith invFs resBasis fimage, fcoker, matB LA.<> fcoker)

-- | Decompose complexes into two-term complexes with differential in a Smith normal form.
smithize :: (Traversable t, Alternative t) => LA.Matrix LA.Z -> t (LA.Matrix LA.Z) -> (t (RepSmith LA.Z), LA.Matrix LA.Z)
smithize headD tailDs =
  uncurry combine $ mapAccumL (uncurry step) (LA.ident (LA.cols headD), headD) tailDs
  where
    step bs hd td = let (res, basis, newhead) = force (squeezeMat bs hd td)
                    in ((basis,newhead), res)
    combine (basis, lastD) ts =
      let (res, coker, _) = squeezeMat basis lastD (LA.ident (LA.rows lastD))
      in (ts <|> pure (force res), coker)

-- | Compute homology group with integral coefficient.
intHomology :: (Traversable t, Alternative t) => LA.Matrix LA.Z -> t (LA.Matrix LA.Z) -> t IntHomology
intHomology headD tailDs =
  let (smiths,coker) = smithize headD tailDs
  in uncurry (combine coker) $ mapAccumL step (IntHomology [] [] []) smiths
  where
    step hdata smith =
      let newcyc = drop (length (invFactor smith)) $ LA.toColumns (basisS smith)
          (ones,tors) = span (==1) (invFactor smith)
          torcycs = LA.toColumns $ LA.dropColumns (length ones) (basisT smith)
          bnds = zipWith LA.scale (invFactor smith) (LA.toColumns (basisT smith))
      in (IntHomology [] (zip torcycs tors) bnds, hdata {freeCycs = newcyc})
    combine coker hdata ts = ts <|> pure hdata {freeCycs = LA.toColumns coker}
