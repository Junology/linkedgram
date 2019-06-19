{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}

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

import qualified Data.List as L

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map

import Data.Matrix (Matrix)
import qualified Data.Matrix as Mat

import Data.Foldable

----------------
-- Data types --
----------------
data FreeMod a b = FMod (Map b a)
  deriving Eq

instance (Show a, Show b) => Show (FreeMod a b) where
  show (FMod mp)
    = let ts = Map.toList mp
      in if L.null ts
         then "@0"
         else L.intercalate " @+ " $ show' <$> ts
    where
      show' (x,r) = show r ++ "@*@" ++ show x

------------------------------
-- Miscellaneous operations --
------------------------------
getCoeff :: (Num a, Ord b) => b -> FreeMod a b -> a
getCoeff x (FMod mp) =
  case Map.lookup x mp of
    Just coeff -> coeff
    Nothing -> fromIntegral 0

spanning :: FreeMod a b -> [b]
spanning (FMod mp) = Map.keys mp

removeZero :: (Num a, Eq a) => FreeMod a b -> FreeMod a b
removeZero (FMod x) = FMod $ Map.filter (/=fromIntegral 0) x

---------------------------------
-- Basic Arithmetic operations --
---------------------------------
infixr 7 @*@
infixr 7 @*@%
infixr 7 @*
infixr 7 @*%
infixl 6 @+
infixl 6 @+%

zeroVec :: FreeMod a b
zeroVec = FMod (Map.empty)

(@*@%) :: a -> b -> FreeMod a b
(@*@%) r x = FMod (Map.singleton x r)

(@*@) r x = removeZero (r @*@% x)

(@*%) :: (Eq a, Num a) => a -> FreeMod a b -> FreeMod a b
(@*%) r (FMod mp) = FMod $ Map.map (r*) mp

(@*) r x = removeZero (r @* x)

(@+%) :: (Eq a, Num a, Ord b) => FreeMod a b -> FreeMod a b -> FreeMod a b
(@+%) (FMod mp) (FMod mq) = FMod $ Map.unionWith (+) mp mq

(@+) x y = removeZero (x @+% y)

---------------------------
-- Monad-like operations --
---------------------------
infixr 1 @=<<
infixr 1 @=<<%

(@=<<%) :: (Eq a, Num a, Ord b, Ord c)
  => (b -> FreeMod a c) -> FreeMod a b -> FreeMod a c
(@=<<%) f (FMod mp)
  = Map.foldlWithKey bin zeroVec mp
  where
    bin x b r = x @+% r@*%f b

(@=<<) :: (Eq a, Num a, Ord b, Ord c)
  => (b -> FreeMod a c) -> FreeMod a b -> FreeMod a c
(@=<<) f fmp
  = removeZero $ f @=<<% fmp

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

-----------------------
-- Tensor operations --
-----------------------
sumFM :: (Eq a, Num a, Ord b, Foldable t) => t (FreeMod a b) -> FreeMod a b
sumFM = removeZero . foldl' (@+%) zeroVec

tensorWith :: (Eq a, Num a, Ord b, Ord c, Ord d)
  => (b -> c -> d) -> FreeMod a b -> FreeMod a c -> FreeMod a d
tensorWith f (FMod mp) y
  = removeZero $ Map.foldlWithKey bin zeroVec mp
  where
    bin x b r = x @+% r@*%( f b @$> y)

tensorFM :: (Eq a, Num a, Ord b, Foldable t) => t (FreeMod a b) -> FreeMod a [b]
tensorFM = foldr' (tensorWith (:)) (1@*@%[])

headConsMapFM :: (Eq a, Num a, Ord b)
  => (b -> FreeMod a [b]) -> FreeMod a [b] -> FreeMod a [b]
headConsMapFM f (FMod mp)
  = removeZero $ Map.foldlWithKey bin zeroVec mp
  where
    bin x [] r = x @+% r@*@%[]
    bin x (b:bs) r
      = x @+% r@*%( (++bs) @$> f b)

insertMap :: (Int -> Either b (a -> b)) -> [a] -> [b]
insertMap f xs = L.unfoldr cmb (xs,[f i|i<-[0..]])
  where
    cmb (_,[]) = undefined
    cmb (ys,(Left b:ebs)) = Just (b,(ys,ebs))
    cmb ([],Right _:_) = Nothing
    cmb (y:ys,Right g:ebs) = Just (g y,(ys,ebs))

insertMapFM :: (Eq a, Num a, Ord b, Ord c)
  => (Int -> Either (FreeMod a c) (b -> FreeMod a c)) -> FreeMod a [b] -> FreeMod a [c]
insertMapFM f (FMod mp)
  = removeZero $ Map.foldlWithKey bin zeroVec mp
  where
    bin x bs r = x @+% r@*%( tensorFM (insertMap f bs) )

mapFM :: (Eq a, Num a, Ord b, Ord c)
  => (b -> FreeMod a c) -> FreeMod a [b] -> FreeMod a [c]
mapFM f = insertMapFM (const (Right f))

---------------------------
-- Matrix representation --
---------------------------
genMatrix :: (Num a, Ord c) => (b -> FreeMod a c) -> [b] -> [c] -> Matrix a
genMatrix f dom cod
  = let el i j = getCoeff (cod !! (i-1)) ((fmap f dom) !! (j-1))
    in Mat.matrix (length cod) (length dom) $ uncurry el

genMatrixAR :: (Num a, Ord c) => (b -> FreeMod a c) -> [b] -> Matrix a
genMatrixAR f dom
  = let cod = L.nub $ L.concatMap spanning $ fmap f dom
    in genMatrix f dom cod
