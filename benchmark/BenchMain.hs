import Criterion.Main

import Data.List as L
import Data.Matrix as Mat

import Numeric.LinearAlgebra

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius
import qualified Numeric.Algebra.IntMatrix.HNFLLL as HASKHNF
import qualified Numeric.Algebra.IntMatrix.HermiteNFLLL as FFIHNF

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
  bgroup "TQFT operation" [
      bench "tqftZ" $ nf (genMatrixAR (tqftZ divides (testSL2B 200))) (subsequences [1..11]) ] ]

