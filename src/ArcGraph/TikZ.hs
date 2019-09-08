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
import Control.Monad.ST

import Data.STRef

import Numeric.Natural (Natural)
import Data.Maybe
import Data.List
import Data.Text (Text)
import qualified Data.Text as T
import Data.Map (Map)
import qualified Data.Map.Strict as Map
import Data.Foldable

import Numeric (showFFloat)

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius

import Text.TeXout
import ArcGraph
import ArcGraph.Component
import ArcGraph.State
import ArcGraph.EnhancedState
import ArcGraph.Common

tikzFloat :: Double -> Text
tikzFloat x = T.pack $ showFFloat (Just 4) x ""

tikzVertex :: Vertex -> Text
tikzVertex (x,y)
  = '(' `T.cons` tikzFloat x <> T.singleton ',' <> tikzFloat y `T.snoc` ')'

drawCommand :: Maybe String -> Text
drawCommand Nothing = macro "draw" [OptArg "very thick, cap=round"]
drawCommand (Just colname)
  = macro "draw" [OptArg (colname ++ ",very thick, cap=round")]

nodeCommand :: Double -> Double -> String -> Text
nodeCommand x y str = T.pack "\\node at " <> tikzVertex (x,y) <> embrace str

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

lineToCommand :: Maybe String -> Vertex -> Vertex -> Text
lineToCommand colname v w
  = drawCommand colname
    <> tikzVertex v <> T.pack " -- " <> tikzVertex w <> T.pack ";\n"

bezierToCommand :: Maybe String -> BezierCube -> Text
bezierToCommand colname (v0, v1, v2, v3)
  = drawCommand colname
    <> T.singleton ' ' <> tikzVertex v0
    <> T.pack " .. controls " <> tikzVertex v1
    <> T.pack " and " <> tikzVertex v2
    <> T.pack " .. " <> tikzVertex v3
    <> T.pack ";\n"

showArcPathTikz :: Maybe String -> ArcPath -> Text
showArcPathTikz colname pth
  = let bs = fmap elevateBezier (mkBezier pth)
    in foldl' (<>) T.empty $ bezierToCommand colname <$> bs

showCrossTikz :: Double -> Cross -> Text
showCrossTikz bd (Crs sega segb crst)
  = let (Sgmt v0@(v0x,v0y) v1@(v1x,v1y)) = sega
        (Sgmt w0@(w0x,w0y) w1@(w1x,w1y)) = segb
        (cx,cy) = fromJust $ calcCross sega segb
    in case crst of
         Crossing ->
           let d0 = distance (w0x,w0y) (cx,cy)
               d1 = distance (w1x,w1y) (cx,cy)
           in lineToCommand Nothing v0 v1
              <> lineToCommand Nothing w0 (cx+bd*(w0x-cx)/d0, cy+bd*(w0y-cy)/d0)
              <> lineToCommand Nothing w1 (cx+bd*(w1x-cx)/d1, cy+bd*(w1y-cy)/d1)
         Smooth0 ->
           let p1 = ((v0x+2.0*cx)/3.0, (v0y+2.0*cy)/3.0)
               p2 = ((w1x+2.0*cx)/3.0, (w1y+2.0*cy)/3.0)
               q1 = ((w0x+2.0*cx)/3.0, (w0y+2.0*cy)/3.0)
               q2 = ((v1x+2.0*cx)/3.0, (v1y+2.0*cy)/3.0)
           in bezierToCommand Nothing (v0,p1,p2,w1)
              <> bezierToCommand Nothing (w0,q1,q2,v1)
         Smooth1 ->
           showCrossTikz bd (Crs (tposeSegment sega) segb Smooth0)

showArcGraphTikz :: Double -> ArcGraph -> Text
showArcGraphTikz bd (AGraph ps cs)
  = foldl' (<>) T.empty
    $ fmap (showArcPathTikz Nothing) ps ++ fmap (showCrossTikz bd) cs

showArcGraphTikzWithComp :: Double -> ArcGraph -> Text
showArcGraphTikzWithComp bd ag@(AGraph _ cs) = runST $ do
  stTeX <- newSTRef T.empty
  forM_ (zip (getComponents ag :: [ArcBits Natural]) (cycle $ map Just tikzColors)) $ \ipsc -> do
    let (ips,c) = ipsc
    for_ (componentAt ag ips) $ \x ->
      modifySTRef' stTeX (<> showArcPathTikz c x)
  forM_ cs $ \crs ->
    modifySTRef' stTeX (<> showCrossTikz bd crs)
  readSTRef stTeX

showArcGraphEnhTikz :: (DState ds, PComponent pc) => Double -> ArcGraph -> ds -> MapEState pc -> Text
showArcGraphEnhTikz bd ag st (MEState coeffMap) = runST $ do
  let normAG = normalize 1.0 ag
  stTeX <- newSTRef T.empty
  modifySTRef' stTeX (<> showArcGraphTikzWithComp bd (smoothing normAG st))
  forM_ (Map.toList coeffMap) $ \kv -> do
    let (comp,coeff) = kv
        (x,y) = fromMaybe (0,0) $ findMostVrtx (\v w -> fst v < fst w) $ componentAt normAG comp
    -- Draw label
    modifySTRef' stTeX (<> (nodeCommand x y (showSL2B coeff) <> T.pack ";\n"))
  readSTRef stTeX
    where
      showSL2B SLI = "$1$"
      showSL2B SLX = "$X$"

showStateSumTikz :: (TeXMathShow a, Num a, Eq a, DState ds, PComponent pc) => ArcGraph -> FreeMod a (ds, MapEState pc) -> Text
showStateSumTikz ag vect =
  if vect == zeroVec
  then T.singleton '0'
  else runST $ do
    let normAG = normalize 1.0 ag
    stTeX <- newSTRef T.empty
    let drawAG coeff (st,enh) = do
          modifySTRef' stTeX (<> texMathShow coeff)
          modifySTRef' stTeX (<> macro "tikz" [OptArg "baseline=-.5ex"])
          modifySTRef' stTeX (<> T.pack "{%\n")
          modifySTRef' stTeX (<> showArcGraphEnhTikz 0.15 normAG st enh)
          modifySTRef' stTeX (<> T.pack "}")
    forEachWithInterM drawAG (modifySTRef' stTeX (<> T.singleton '+')) vect
    readSTRef stTeX

chunksOfMap :: (Ord k) => Int -> Map k a -> [Map k a]
chunksOfMap n = unfoldr $ \mp ->
  if Map.null mp
  then Nothing
  else Just (Map.splitAt n mp)

showStateSumTikzMultlined :: (TeXMathShow a, Num a, Eq a, DState ds, PComponent pc) => ArcGraph -> FreeMod a (ds, MapEState pc) -> Int -> (Bool,Text)
showStateSumTikzMultlined ag vect n =
  if vect == zeroVec
  then (False,T.singleton '0')
  else 
    let termLines = chunksOfMap n (termMap vect)
        textLines = fmap (showStateSumTikz ag . FMod) termLines
    in (length textLines > 1, T.intercalate (T.pack "\\\\\n+") textLines)
