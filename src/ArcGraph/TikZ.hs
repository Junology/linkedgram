{-# LANGUAGE ScopedTypeVariables #-}

------------------------------------------------
-- |
-- Module    :  ArcGraph.TikZ
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- LaTeX + TikZ output of ArcGraph
--
------------------------------------------------

module ArcGraph.TikZ where

import Control.Monad
import Control.Monad.State.Strict
import Control.Monad.ST

import Data.STRef

import Data.Maybe
import Data.List
import Data.Text (Text)
import qualified Data.Text as T
import qualified Data.Map.Strict as Map
import Data.Foldable

import Numeric

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius

import Text.TeXout
import ArcGraph
import ArcGraph.Component
import ArcGraph.State
import ArcGraph.EnhancedState
import ArcGraph.Common

texHeader :: String -> String -> String
texHeader option cls
  = "\\documentclass[pdftex"
  ++ (if null option then "" else ",")
  ++ option ++ "]{" ++ cls ++ "}\n"

texPreamble :: String
texPreamble = "\\usepackage{amsmath,amssymb}\n\\usepackage{tikz}\n\n\\allowdisplaybreaks[2]\n\n"

texBegin :: String
texBegin = "\\begin{document}\n"

texEnd :: String
texEnd = "\\end{document}\n"

texParagraph :: String -> String
texParagraph str = "\\paragraph{" ++ str ++ "}\n\\mbox{}\\\\\\leavevmode\n"

texHorizontalLine :: String
texHorizontalLine
  = "\\leavevmode\n\\\n\\noindent\\makebox[\\linewidth]{\\rule{\\paperwidth}{0.4pt}}\n"

tikzFloat :: Double -> String
tikzFloat x = showFFloat (Just 4) x ""

tikzVertex :: Vertex -> String
tikzVertex (x,y) = "(" ++ tikzFloat x ++ "," ++ tikzFloat y ++ ")"

tikzBegin :: Double -> String
tikzBegin scale = "\\begin{tikzpicture}[scale=" ++ tikzFloat scale ++ "]\n"

tikzEnd :: String
tikzEnd = "\\end{tikzpicture}\n"

drawCommand :: Maybe String -> String
drawCommand Nothing = "\\draw[very thick, cap=round]"
drawCommand (Just colname) = "\\draw[" ++ colname ++ ",very thick, cap=round]"

nodeCommand :: Double -> Double -> String -> String
nodeCommand x y str = "\\node at " ++ tikzVertex (x,y) ++ "{" ++ str ++ "};"

tikzColors :: [String]
tikzColors = ["red",
              "green",
              "blue",
              "cyan",
              "magenta",
              "yellow",
              -- "black",
              "gray",
              "darkgray",
              "lightgray",
              "brown",
              "lime",
              "olive",
              "orange",
              "pink",
              "purple",
              "teal",
              "violet"] -- and "white"

lineToCommand :: Maybe String -> Vertex -> Vertex -> String
lineToCommand colname v w
  = drawCommand colname ++ tikzVertex v ++ " -- " ++ tikzVertex w ++ ";\n"

bezierToCommand :: Maybe String -> BezierCube -> String
bezierToCommand colname (v0, v1, v2, v3)
  = drawCommand colname
    ++ " " ++ tikzVertex v0
    ++ " .. controls " ++ tikzVertex v1
    ++ " and " ++ tikzVertex v2
    ++ " .. " ++ tikzVertex v3
    ++ ";\n"

showArcPathTikz :: Maybe String -> ArcPath -> String
showArcPathTikz colname pth
  = let bs = fmap elevateBezier (mkBezier pth)
    in concatMap (bezierToCommand colname) bs

showCrossTikz :: Double -> Cross -> String
showCrossTikz bd (Crs sega segb crst)
  = let (Sgmt v0@(v0x,v0y) v1@(v1x,v1y)) = sega
        (Sgmt w0@(w0x,w0y) w1@(w1x,w1y)) = segb
        (cx,cy) = fromJust $ calcCross sega segb
    in case crst of
         Crossing -> flip execState "" $ do
           modify' (++lineToCommand Nothing v0 v1)
           let d0 = distance (w0x,w0y) (cx,cy)
           modify' (++lineToCommand Nothing w0 (cx+bd*(w0x-cx)/d0, cy+bd*(w0y-cy)/d0))
           let d1 = distance (w1x,w1y) (cx,cy)
           modify' (++lineToCommand Nothing w1 (cx+bd*(w1x-cx)/d1, cy+bd*(w1y-cy)/d1))
         Smooth0 -> flip execState "" $ do
           let p1 = ((v0x+2.0*cx)/3.0, (v0y+2.0*cy)/3.0)
               p2 = ((w1x+2.0*cx)/3.0, (w1y+2.0*cy)/3.0)
               q1 = ((w0x+2.0*cx)/3.0, (w0y+2.0*cy)/3.0)
               q2 = ((v1x+2.0*cx)/3.0, (v1y+2.0*cy)/3.0)
           modify' (++bezierToCommand Nothing (v0,p1,p2,w1))
           modify' (++bezierToCommand Nothing (w0,q1,q2,v1))
         Smooth1 ->
           showCrossTikz bd (Crs (tposeSegment sega) segb Smooth0)

showArcGraphTikz :: Double -> ArcGraph -> String
showArcGraphTikz bd (AGraph ps cs) = flip execState "" $ do
  for_ ps $ \path ->
    modify' (++ showArcPathTikz Nothing path)
  modify' (++ concatMap (showCrossTikz bd) cs)

showArcGraphTikzWithComp :: Double -> ArcGraph -> String
showArcGraphTikzWithComp bd ag@(AGraph _ cs) = flip execState "" $ do
  for_ (zip (getComponents ag :: [ArcList]) (cycle $ map Just tikzColors)) $ \ipsc -> do
    let (ips,c) = ipsc
    for_ (componentAt ag ips) $ \x ->
      modify' (++ showArcPathTikz c x)
  modify' (++ concatMap (showCrossTikz bd) cs)

showArcGraphEnhTikz :: (DState ds, PComponent pc) => Double -> ArcGraph -> ds -> MapEState pc -> String
showArcGraphEnhTikz bd ag st (MEState coeffMap) = flip execState "" $ do
  let normAG = normalize 1.0 ag
  modify' (++ showArcGraphTikzWithComp bd (smoothing normAG st))
  forM_ (Map.toList coeffMap) $ \kv -> do
    let (comp,coeff) = kv
        (x,y) = fromMaybe (0,0) $ findMostVrtx (\v w -> fst v < fst w) $ componentAt normAG comp
    -- Draw label
    modify' (++ (nodeCommand x y $ case coeff of {SLI -> "$1$"; SLX -> "$X$";}))

showStateSumTikzText :: (DState ds, PComponent pc) => ArcGraph -> FreeMod Int (ds, MapEState pc) -> Text
showStateSumTikzText ag vect =
  if vect == zeroVec
  then T.singleton '0'
  else runST $ do
    let normAG = normalize 1.0 ag
    stTeX <- newSTRef T.empty
    let drawAG coeff (st,enh) = do
          modifySTRef' stTeX (<> texMathShow coeff)
          modifySTRef' stTeX (<> macro "tikz" [OptArg "baseline=-.5ex"])
          modifySTRef' stTeX (<> T.pack "{%\n")
          modifySTRef' stTeX (<> T.pack (showArcGraphEnhTikz 0.15 normAG st enh))
          modifySTRef' stTeX (<> T.pack "}")
    forEachWithInterM drawAG (modifySTRef' stTeX (<> T.singleton '+')) vect
    readSTRef stTeX

showStateSumTikz :: (DState ds, PComponent pc) => Double -> ArcGraph -> FreeMod Int (ds, MapEState pc) -> String
showStateSumTikz bd ag vect = flip execState "" $
  if vect == zeroVec
  then modify' (++ "\\[0\\]\n")
  else do
    modify' (++"\\begin{multline*}\n")
    forEachWithInterM drawAG (modify' (++"+")) vect
    modify' (++"\\end{multline*}\n")
      where
        drawAG :: (DState ds, PComponent pc) => Int -> (ds, MapEState pc) -> State String ()
        drawAG coeff (st,enh) = do
          modify' (++ show coeff)
          modify' (++"\\tikz{%\n")
          modify' (++ showArcGraphEnhTikz bd ag st enh)
          modify' (++"}\n")
