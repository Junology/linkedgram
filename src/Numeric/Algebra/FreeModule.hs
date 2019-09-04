{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}

------------------------------------------------
-- |
-- Module    : Numeric.Algebra.FreeModule
-- Copyright : (c) Jun Yoshida 2019
-- License   : BSD3
--
-- Free modules over Num type.
--
------------------------------------------------

module Numeric.Algebra.FreeModule where

import GHC.Generics

import Control.DeepSeq

import Control.Monad
import Control.Monad.ST

import qualified Data.List as L
import qualified Data.Maybe as M

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map

import Data.Matrix (Matrix)
import qualified Data.Matrix as Mat

import qualified Numeric.LinearAlgebra as LA
import qualified Numeric.LinearAlgebra.Devel as LA

import Data.Foldable

------------------
-- * Data types --
------------------
newtype FreeMod a b = FMod { termMap :: Map b a }
  deriving (Eq, Ord, Generic, NFData)

instance (Show a, Show b) => Show (FreeMod a b) where
  show (FMod mp)
    = let ts = Map.toList mp
      in if L.null ts
         then "@0"
         else L.intercalate " @+ " $ show' <$> ts
    where
      show' (x,r) = show r ++ "@*@" ++ show x

--------------------------------
-- * Miscellaneous operations --
--------------------------------
getCoeff :: (Num a, Ord b) => b -> FreeMod a b -> a
getCoeff x = M.fromMaybe 0 . Map.lookup x . termMap

spanning :: FreeMod a b -> [b]
spanning = Map.keys . termMap

removeZero :: (Num a, Eq a) => FreeMod a b -> FreeMod a b
removeZero = FMod . Map.filter (/=0) . termMap

forEachWithInterM :: Monad m => (a -> b -> m ()) -> m () -> FreeMod a b -> m ()
forEachWithInterM f g (FMod mp) =
  let mpl = Map.toList mp
  in forM_ (L.intersperse Nothing (map Just mpl)) $ maybe g (uncurry (flip f))

-----------------------------------
-- * Basic Arithmetic operations --
-----------------------------------
infixr 7 @*@
infixr 7 @*@%
infixr 7 @*
infixr 7 @*%
infixl 6 @+
infixl 6 @+%

zeroVec :: FreeMod a b
zeroVec = FMod Map.empty

(@*@%) :: a -> b -> FreeMod a b
(@*@%) r x = FMod (Map.singleton x r)

(@*@) :: (Eq a, Num a) => a -> b -> FreeMod a b
(@*@) r x = removeZero (r @*@% x)

(@*%) :: (Num a) => a -> FreeMod a b -> FreeMod a b
(@*%) r (FMod mp) = FMod $ Map.map (r*) mp

(@*) :: (Eq a, Num a) => a -> FreeMod a b -> FreeMod a b
(@*) r x = removeZero (r @* x)

(@+%) :: (Num a, Ord b) => FreeMod a b -> FreeMod a b -> FreeMod a b
(@+%) (FMod mp) (FMod mq) = FMod $ Map.unionWith (+) mp mq

(@+) :: (Eq a, Num a, Ord b) => FreeMod a b -> FreeMod a b -> FreeMod a b
(@+) x y = removeZero (x @+% y)

------------------------
-- * Monad operations --
------------------------
infixr 1 @=<<
infixr 1 @=<<%

(@=<<%) :: (Eq a, Num a, Ord b, Ord c)
  => (b -> FreeMod a c) -> FreeMod a b -> FreeMod a c
(@=<<%) f (FMod mp)
  = Map.foldlWithKey bin zeroVec mp
  where
    bin x b r = x @+% r@*%f b

(@>>=%) :: (Eq a, Num a, Ord b, Ord c)
  => FreeMod a b -> (b -> FreeMod a c) -> FreeMod a c
(@>>=%) = flip (@=<<%)

(@=<<) :: (Eq a, Num a, Ord b, Ord c)
  => (b -> FreeMod a c) -> FreeMod a b -> FreeMod a c
(@=<<) f fmp
  = removeZero $ f @=<<% fmp

(@>>=) :: (Eq a, Num a, Ord b, Ord c)
  => FreeMod a b -> (b -> FreeMod a c) -> FreeMod a c
(@>>=) = flip (@=<<)

infixr 1 @>=>
infixr 1 @>=>%

(@>=>%) :: (Eq a, Num a, Ord b, Ord c, Ord d)
  => (b -> FreeMod a c) -> (c -> FreeMod a d) -> b -> FreeMod a d
(@>=>%) f g = (g@=<<%) . f

(@>=>) :: (Eq a, Num a, Ord b, Ord c, Ord d)
  => (b -> FreeMod a c) -> (c -> FreeMod a d) -> b -> FreeMod a d
(@>=>) f g = (g@=<<) . f

infixr 1 @$>
infixr 1 @$>%

(@$>%) :: (Eq a, Num a, Ord b, Ord c)
  => (b -> c) -> FreeMod a b -> FreeMod a c
(@$>%) f = (@=<<%) ((1@*@%) . f)

(@$>) :: (Eq a, Num a, Ord b, Ord c)
  => (b -> c) -> FreeMod a b -> FreeMod a c
(@$>) f = removeZero . (f@$>%)

-------------------------
-- * Tensor operations --
-------------------------
sumFM' :: (Num a, Ord b, Foldable t) => t (FreeMod a b) -> FreeMod a b
sumFM' = foldl' (@+%) zeroVec

sumFM :: (Eq a, Num a, Ord b, Foldable t) => t (FreeMod a b) -> FreeMod a b
sumFM = removeZero . sumFM'

zipSum' :: (Foldable t, Applicative t, Num a, Ord b)
  => t a -> t b -> FreeMod a b
zipSum' rs es = foldl' (@+%) zeroVec $ zipWith (@*@%) (toList rs) (toList es)

zipSum :: (Eq a, Foldable t, Applicative t, Num a, Ord b)
  => t a -> t b -> FreeMod a b
zipSum rs es = removeZero $ zipSum' rs es

tensorWith :: (Eq a, Num a, Ord b, Ord c, Ord d)
  => (b -> c -> d) -> FreeMod a b -> FreeMod a c -> FreeMod a d
tensorWith f (FMod mp) y
  = removeZero $ Map.foldlWithKey bin zeroVec mp
  where
    bin x b r = x @+% r@*%( f b @$> y)

tensorFM :: (Eq a, Num a, Ord b, Foldable t) => t (FreeMod a b) -> FreeMod a [b]
tensorFM = foldr' (tensorWith (:)) (1@*@%[])

tensorMapFM :: (Eq a, Num a, Ord b, Ord c)
  => [b -> FreeMod a c] -> FreeMod a [b] -> FreeMod a [c]
tensorMapFM fs = (@=<<) (tensorFM . zipWith ($) fs)

insertMap :: (Int -> Either b (a -> b)) -> [a] -> [b]
insertMap f xs = L.unfoldr cmb (xs,[f i|i<-[0..]])
  where
    cmb (_,[]) = undefined
    cmb (ys, Left b:ebs) = Just (b,(ys,ebs))
    cmb ([],Right _:_) = Nothing
    cmb (y:ys,Right g:ebs) = Just (g y,(ys,ebs))

insertMapFM :: (Eq a, Num a, Ord b, Ord c)
  => (Int -> Either (FreeMod a c) (b -> FreeMod a c)) -> FreeMod a [b] -> FreeMod a [c]
insertMapFM f (FMod mp)
  = removeZero $ Map.foldlWithKey bin zeroVec mp
  where
    bin x bs r = x @+% r@*% tensorFM (insertMap f bs)

mapFM :: (Eq a, Num a, Ord b, Ord c)
  => (b -> FreeMod a c) -> FreeMod a [b] -> FreeMod a [c]
mapFM f = insertMapFM (const (Right f))

-----------------------------
-- * Matrix representation --
-----------------------------
genMatrix :: (Num a, Ord c) => (b -> FreeMod a c) -> [b] -> [c] -> Matrix a
genMatrix _ [] _ = Mat.fromList 0 0 []
genMatrix f dom cod
  = let domRk = length dom
        codRk = length cod
        el j i = if i <= domRk && j <= codRk
                 then getCoeff (cod !! (j-1)) (f (dom !! (i-1)))
                 else 0
    in Mat.matrix codRk domRk $ uncurry el

genMatrixAR :: (Num a, Ord c) => (b -> FreeMod a c) -> [b] -> Matrix a
genMatrixAR f dom
  = let cod = L.nub $ L.concatMap spanning $ fmap f dom
    in genMatrix f dom cod

genMatrixILA :: (Integral a, Num a', LA.Element a', Ord c) => (b -> FreeMod a c) -> [b] -> [c] -> LA.Matrix a'
genMatrixILA f dom cod
  = let domRk = length dom
        codRk = length cod
    in runST $ do
  stMat <- LA.newUndefinedMatrix LA.RowMajor codRk domRk
  forM_ [0..codRk-1] $ \i ->
    forM_ [0..domRk-1] $ \j -> do
      let coeff = getCoeff (cod !! i) (f (dom !! j))
      LA.unsafeWriteMatrix stMat i j (fromIntegral coeff)
  LA.unsafeFreezeMatrix stMat
