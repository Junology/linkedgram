module Main where

import Criterion.Main

import Data.List as L
import Data.Set as Set

testlist :: [Int]
testlist = [0..9]

naiveSubseq :: Int -> [Int] -> [[Int]]
naiveSubseq n = L.filter ((==n) . length) . subsequences

yanaiveSubseq :: Int -> [Int] -> [[Int]]
yanaiveSubseq n []
  | n == 0    = [[]]
  | otherwise = []
yanaiveSubseq 0 (_:_) = [[]]
yanaiveSubseq n (x:xs) = yanaiveSubseq n xs ++ fmap (x:) (yanaiveSubseq (n-1) xs)

setSubseq :: Int -> [Int] -> [[Int]]
setSubseq n = Set.toList . Set.mapMonotonic Set.toList . Set.filter ((==n) . Set.size) . powerSet . fromList

main :: IO ()
main = defaultMain [
  bgroup "subseq" [
      bench "L.subsequences" $ nf (naiveSubseq 5) testlist,
      bench "yanaiveSubseq" $ nf (yanaiveSubseq 5) testlist,
      bench "Set.powerSet" $ nf (setSubseq 5) testlist] ]

{-- Benchmark for matrix normal forms.
import Data.List as L
import Data.Matrix as Mat

import Numeric.LinearAlgebra

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius
import qualified Numeric.Algebra.IntMatrix.HNFLLL as HASKHNF
import qualified Numeric.Algebra.IntMatrix.NormalForms as FFIHNF

import qualified Numeric.Algebra.IntMatrix.SmithNF as HASKSMITH


testmat1 = (100 >< 100) $ cycle [-2,-1,0,1,2,3,4]

hasIntersection :: Eq a => [a] -> [a] -> Bool
hasIntersection x y = not $ L.null (x `L.intersect` y)

testSL2B :: Int -> [(Int,SL2B)]
testSL2B n = zip [1..n] (cycle [SLI,SLI,SLX])

divides :: Int -> Int -> Bool
divides m n = n `rem` m == 0

main :: IO ()
main = defaultMain [
  bgroup "hermiteNF" [
      bench "pure Haskell" $ nf HASKHNF.hnfLLL testmat1,
      bench "ffi" $ nf FFIHNF.hermiteNF testmat1 ],
  bgroup "smithNF" [
      bench "Haskell" $ nf HASKSMITH.smithNF testmat1,
      bench "ffi" $ nf FFIHNF.smithNF testmat1 ],
  bgroup "TQFT operation" [
      bench "tqftZ" $ nf (genMatrixAR (tqftZ divides (testSL2B 200))) (subsequences [1..11]) ] ]
--}
