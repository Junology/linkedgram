------------------------------------------------
-- |
-- Module    :  ArcGraph.TikZ
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- LaTeX + TikZ output of ArcGraph
--
------------------------------------------------

module ArcGraph.TikZ (
  showArcGraphTikz,
  typesetArcGraphTikz,
  docArcGraphTikz
  ) where

import Control.Monad
import Control.Monad.State.Strict

import Data.Maybe
import Data.List

import qualified Graphics.Rendering.Cairo as Cairo

import ArcGraph
import ArcGraph.Common

texHeader :: String -> String -> String
texHeader option cls
  = "\\documentclass[pdftex"
  ++ (if null option then "" else ",")
  ++ option ++ "]{" ++ cls ++ "}\n"

texPreamble :: String
texPreamble = "\\usepackage{tikz}\n"

texBegin :: String
texBegin = "\\begin{document}\n"

texEnd :: String
texEnd = "\\end{document}\n"

texHorizontalLine :: String
texHorizontalLine
  = "\\noindent\\makebox[\\linewidth]{\\rule{\\paperwidth}{0.4pt}}\n"

tikzBegin :: String
tikzBegin = "\\begin{tikzpicture}\n"

tikzEnd :: String
tikzEnd = "\\end{tikzpicture}\n"

drawCommand :: String
drawCommand = "\\draw[ultra thick, cap=round]"

lineToCommand :: Vertex -> Vertex -> String
lineToCommand v w = drawCommand ++ show v ++ " -- " ++ show w ++ ";\n"

bezierToCommand :: BezierCube -> String
bezierToCommand (v0, v1, v2, v3)
  = drawCommand
    ++ " " ++ show v0
    ++ " .. controls " ++ show v1
    ++ " and " ++ show v2
    ++ " .. " ++ show v3
    ++ ";\n"

showArcPathTikz :: ArcPath -> String
showArcPathTikz pth
  = let bs = fmap elevateBezier (mkBezier pth)
    in foldl (++) "" $ fmap bezierToCommand bs

showCrossTikz :: Double -> Cross -> String
showCrossTikz bd (Crs sega segb crst)
  = let (Sgmt v0@(v0x,v0y) v1@(v1x,v1y)) = sega
        (Sgmt w0@(w0x,w0y) w1@(w1x,w1y)) = segb
        (cx,cy) = fromJust $ calcCross sega segb
    in case crst of
         Crossing -> flip execState "" $ do
           modify' (++lineToCommand v0 v1)
           let d0 = distance (w0x,w0y) (cx,cy)
           modify' (++lineToCommand w0 (cx+bd*(w0x-cx)/d0, cy+bd*(w0y-cy)/d0))
           let d1 = distance (w1x,w1y) (cx,cy)
           modify' (++lineToCommand w1 (cx+bd*(w1x-cx)/d1, cy+bd*(w1y-cy)/d1))
         Smooth0 -> flip execState "" $ do
           let p1 = ((v0x+2.0*cx)/3.0, (v0y+2.0*cy)/3.0)
               p2 = ((w1x+2.0*cx)/3.0, (w1y+2.0*cy)/3.0)
               q1 = ((w0x+2.0*cx)/3.0, (w0y+2.0*cy)/3.0)
               q2 = ((v1x+2.0*cx)/3.0, (v1y+2.0*cy)/3.0)
           modify' (++bezierToCommand (v0,p1,p2,w1))
           modify' (++bezierToCommand (w0,q1,q2,v1))
         Smooth1 -> do
           showCrossTikz bd (Crs (tposeSegment sega) segb Smooth0)

showArcGraphTikz :: Double -> ArcGraph -> String
showArcGraphTikz bd ag = flip execState "" $ do
  let (AGraph ps cs) = normalize ag
  modify' (++tikzBegin)
  modify' (++ concat (map showArcPathTikz ps))
  modify' (++ concat (map (showCrossTikz bd) cs))
  modify' (++tikzEnd)

typesetArcGraphTikz :: String -> String -> [ArcGraph] -> String
typesetArcGraphTikz option cls ags = flip execState "" $ do
  modify' (++texHeader option cls)
  modify' (++texPreamble)
  modify' (++texBegin)
  forM_ ags $ \ag -> do
    modify' (++showArcGraphTikz 0.15 ag)
  modify' (++texEnd)

docArcGraphTikz :: String -> String -> [ArcGraph] -> String
docArcGraphTikz option cls ags = flip execState "" $ do
  modify' (++texHeader option cls)
  modify' (++texPreamble)
  modify' (++texBegin)
  forM_ ags $ \ag -> do
    modify' (++"\n\\section{}\n")
    modify' (++"\\begin{center}\n")
    modify' (++showArcGraphTikz 0.15 ag)
    modify' (++"\\end{center}\n")
    forM_ [0..countCross ag] $ \i -> do
      modify' (++texHorizontalLine)
      modify' (++"\\subsection*{$H^{" ++ show i ++ ",\\star}$}\n")
      -- modify' (++texHorizontalLine)
      modify' (++"\\begin{center}\n")
      forM_ (listSmoothing [i] ag) $ \agsm -> do
        modify' (++showArcGraphTikz 0.15 agsm)
        modify' (++"\\hfill\n")
      modify' (++"\\end{center}\n")
  modify' (++texEnd)