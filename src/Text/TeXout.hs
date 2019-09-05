{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE FlexibleInstances, UndecidableInstances #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}

------------------------------------------------
-- |
-- Module    :  Text.TeXout
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Export TeX codes.
--
------------------------------------------------

module Text.TeXout where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData)

import Control.Monad
import Control.Monad.ST
import Data.STRef

import Data.Text (Text)
import qualified Data.Text as T

import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map

import Numeric (showFFloat)

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius

-- | Class of types which are exhibited in LaTeX display maths.
class TeXMathShow a where
  texMathShow :: a -> Text

-- * Instances
instance (TeXMathShow a) => TeXMathShow (Maybe a) where
  texMathShow = maybe T.empty texMathShow

instance TeXMathShow String where
  texMathShow = T.pack

instance TeXMathShow Int where
  texMathShow = T.pack . show

instance TeXMathShow Double where
  texMathShow x = T.pack $ showFFloat (Just 4) x ""

instance TeXMathShow SL2B where
  texMathShow SLI = T.singleton '1'
  texMathShow SLX = T.singleton 'X'

instance (TeXMathShow a, TeXMathShow b) => TeXMathShow (a,b) where
  texMathShow (x,y)
    = T.pack "\\left(" <>
      texMathShow x <> T.pack "," <> texMathShow y <>
      T.pack "\\right)"

instance (TeXMathShow a, Ord b, TeXMathShow b) => TeXMathShow (FreeMod a b) where
  texMathShow fm
    = case Map.toList (termMap fm) of
        [] -> texMathShow (0 :: Int)
        ts -> T.intercalate (T.singleton '+') (fmap showProd ts)
    where
      showProd (term,coeff) = texMathShow coeff <> texMathShow term

instance (TeXMathShow a) => TeXMathShow (Map (Int,Int) a) where
  texMathShow mp
    = case foldl rangeFinder Nothing (Map.keys mp) of
        Nothing -> T.empty
        (Just (imin,imax,jmin,jmax))
          -> runST $ do
          stTeX <- newSTRef T.empty
          writeSTRef stTeX $ beginEnvLn "array" [FixArg ("r|" ++ replicate (jmax-jmin+1) 'c')]
          modifySTRef' stTeX (<> mkArrayLn "i\\backslash j" [jmin..jmax])
          modifySTRef' stTeX (<> T.pack "\\\\\\hline\n")
          forM_ [imin..imax] $ \i -> do
            let arrayLn = mkArrayLn i (fmap (\j-> mp Map.!?(i,j)) [jmin..jmax])
            modifySTRef stTeX (<> arrayLn)
            if i < imax
              then (modifySTRef' stTeX (<> T.pack "\\\\\n"))
              else (modifySTRef' stTeX (<> T.singleton '\n'))
          modifySTRef' stTeX (<> endEnv "array")
          readSTRef stTeX
    where
      rangeFinder Nothing (i,j) = Just (i,i,j,j)
      rangeFinder (Just (imin,imax,jmin,jmax)) (i,j) =
        let imin' = if i < imin then i else imin
            imax' = if i > imax then i else imax
            jmin' = if j < jmin then j else jmin
            jmax' = if j > jmax then j else jmax
        in Just (imin',imax',jmin',jmax')

------------------------------
-- * Macros argment control
------------------------------
embrace :: String -> Text
embrace str = '{' `T.cons` T.pack str `T.snoc` '}'

data TeXArg = OptArg String | FixArg String
  deriving (Eq, Ord, Generic, NFData)

encloseArg :: TeXArg -> Text
encloseArg (OptArg arg)
  = case arg of
      []        -> T.empty
      otherwise -> '[' `T.cons` T.pack arg `T.snoc` ']'
encloseArg (FixArg arg) = embrace arg

encloseArgs :: (Foldable t) => t TeXArg -> Text
encloseArgs = foldl (\enc arg -> enc <> encloseArg arg) T.empty


----------------------------
-- * Headers on TeX source
----------------------------
argEnclose :: String -> Text
argEnclose arg = '{' `T.cons` T.pack arg `T.snoc` '}'

optArgEnclose :: String -> Text
optArgEnclose []  = T.empty
optArgEnclose arg = '[' `T.cons` T.pack arg `T.snoc` ']'

documentClass :: String -> String -> Text
documentClass opts cls
  = T.pack "\\documentclass" <> optArgEnclose opts <> argEnclose cls

usePackage :: String -> String -> Text
usePackage opts pkg
  = T.pack "\\usepackage" <> optArgEnclose opts <> argEnclose pkg

------------------------------
-- * Macros and Environments
------------------------------
macro :: (Foldable t) => String -> t TeXArg -> Text
macro csname args
  = '\\' `T.cons` T.pack csname <> encloseArgs args

beginEnv :: (Foldable t) => String -> t TeXArg -> Text
beginEnv env args
  = T.pack "\\begin" <> encloseArg (FixArg env) <> encloseArgs args

beginEnvLn :: (Foldable t) => String -> t TeXArg -> Text
beginEnvLn env args = beginEnv env args `T.snoc` '\n'

endEnv :: String -> Text
endEnv env
  = T.pack "\\end" <> encloseArg (FixArg env)

endEnvLn :: String -> Text
endEnvLn env
  = endEnv env `T.snoc` '\n'

---------------------
-- * Miscellaneous
---------------------
horizontalLine :: Text
horizontalLine
  = T.pack "\\leavevmode\n\\\n\\noindent\\makebox[\\linewidth]{\\rule{\\paperwidth}{0.4pt}}\n"

mkArrayLn :: (TeXMathShow a, TeXMathShow b) => a -> [b] -> Text
mkArrayLn x ys
  = T.intercalate (T.pack " & ") (texMathShow x: fmap texMathShow ys)

mathBB :: String -> Text
mathBB c = macro "mathbb" [FixArg c]

cyclicGrp :: Int -> Text
cyclicGrp n = mathBB "Z" <> T.singleton '/' <> texMathShow n

ofRank :: (Integral a, TeXMathShow a) => Text -> a -> Text
ofRank txt 1 = txt
ofRank txt n = txt <> T.pack "^{\\oplus " <> texMathShow n `T.snoc` '}'

flatZip :: Eq a => [a] -> [(Int,a)]
flatZip [] = []
flatZip xs@(_:_) = uncurry (:) $ foldr bin ((1,last xs),[]) (init xs)
  where
    bin y ((n,z),zs)
      | y==z      = ((n+1,z),zs)
      | otherwise = ((1,y),(n,z):zs)

finAbGroup :: Int -> [Int] -> Text
finAbGroup freeRk torsions =
  let dsum = if freeRk > 0 && not (null torsions)
             then T.pack "\\oplus "
             else T.empty
      freepart = if freeRk > 0
                 then mathBB "Z" `ofRank` freeRk
                 else T.empty
      torGrps = fmap (\tor -> cyclicGrp (snd tor) `ofRank` fst tor) (flatZip torsions)
      torpart = T.intercalate (T.pack "\\oplus ") torGrps
  in freepart <> dsum <> torpart
