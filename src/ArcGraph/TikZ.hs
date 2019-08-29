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
  showStateSumTikz,
  docKhovanovTikz,
  typesetArcGraphTikz,
  docArcGraphTikz
  ) where

import Control.Monad
import Control.Monad.State.Strict

import Data.Maybe
import Data.List
import qualified Data.Map.Strict as Map
import Data.Foldable

import Numeric

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius

import ArcGraph
import ArcGraph.EnhancedState
import ArcGraph.Common

texHeader :: String -> String -> String
texHeader option cls
  = "\\documentclass[pdftex"
  ++ (if null option then "" else ",")
  ++ option ++ "]{" ++ cls ++ "}\n"

texPreamble :: String
texPreamble = "\\usepackage{amsmath,amssymb}\n\\usepackage{tikz}\n\\usepackage{breqn}\n"

texBegin :: String
texBegin = "\\begin{document}\n"

texEnd :: String
texEnd = "\\end{document}\n"

texParagraph :: String -> String
texParagraph str = "\\paragraph{" ++ str ++ "}\n\\mbox\\\\\\leavevmode\n"

texHorizontalLine :: String
texHorizontalLine
  = "\\leavevmode\n\\\\\n\\noindent\\makebox[\\linewidth]{\\rule{\\paperwidth}{0.4pt}}\n"

tikzFloat :: Double -> String
tikzFloat x = showFFloat (Just 4) x ""

tikzVertex :: Vertex -> String
tikzVertex (x,y) = "(" ++ tikzFloat x ++ "," ++ tikzFloat y ++ ")"

tikzBegin :: Double -> String
tikzBegin scale = "\\begin{tikzpicture}[scale=" ++ tikzFloat scale ++ "]\n"

tikzEnd :: String
tikzEnd = "\\end{tikzpicture}\n"

drawCommand :: Maybe String -> String
drawCommand Nothing = "\\draw[ultra thick, cap=round]"
drawCommand (Just colname) = "\\draw[" ++ colname ++ ",ultra thick, cap=round]"

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
showArcGraphTikz bd ag = flip execState "" $ do
  let (AGraph ps cs) = slimCross $ normalize 1.0 ag
  for_ ps $ \path ->
    modify' (++ showArcPathTikz Nothing path)
  modify' (++ concatMap (showCrossTikz bd) cs)

showArcGraphTikzWithComp :: Double -> ArcGraph -> DiagramState -> String
showArcGraphTikzWithComp bd ag st = flip execState "" $ do
  let ag'@(AGraph ps cs) = smoothing st $ normalize 1.0 ag
  for_ (zip (components ag') (cycle $ map Just tikzColors)) $ \ipsc -> do
    let (ips,c) = ipsc
    for_ ips $ \i ->
      modify' (++ showArcPathTikz c (ps!!i))
  modify' (++ concatMap (showCrossTikz bd) cs)

showArcGraphEnhTikz :: Double -> ArcGraphE -> String
showArcGraphEnhTikz bd (AGraphE ag st coeffMap) = flip execState "" $ do
  let normAG@(AGraph ps _) = normalize 1.0 ag
  modify' (++ showArcGraphTikzWithComp bd normAG st)
  forM_ (Map.toList coeffMap) $ \kv -> do
    let (comp,coeff) = kv
        (x,y) = fromMaybe (0,0) $ findMostVrtx (\v w -> fst v < fst w) $ map (ps!!) comp
    -- Draw label
    modify' (++ (nodeCommand x y $ case coeff of {SLI -> "$1$"; SLX -> "$X$";}))

showStateSumTikz :: Double -> FreeMod Int ArcGraphE -> String
showStateSumTikz bd vect = flip execState "" $ do
  modify' (++"\\begin{dmath*}\n")
  forEachWithInterM drawAG (modify' (++"+")) vect
  modify' (++"\\end{dmath*}\n")
  where
    drawAG :: Int -> ArcGraphE -> State String ()
    drawAG coeff agE = do
      modify' (++ show coeff)
      modify' (++"\\tikz{%\n")
      modify' (++ showArcGraphEnhTikz bd agE)
      modify' (++"}\n")

flatZip :: Eq a => [a] -> [(Int,a)]
flatZip [] = []
flatZip xs@(_:_) = uncurry (:) $ foldr bin ((1,last xs),[]) (init xs)
  where
    bin y ((n,z),zs)
      | y==z      = ((n+1,z),zs)
      | otherwise = ((1,y),(n,z):zs)

showAbGroupTeX :: Int -> [Int] -> String
showAbGroupTeX freeRk torsions =
  case torsions of
    []{- no torsion -}
      | freeRk > 0 -> freepart freeRk
      | otherwise  -> "0"
    _ {- has torsions -}
      | freeRk > 0 -> freepart freeRk ++ "\\oplus " ++ torpart torsions
      | otherwise  -> torpart torsions
  where
    freepart rk
      | rk == 1   = "\\mathbb Z"
      | otherwise = "\\mathbb Z^{\\oplus " ++ show rk ++ "}"
    markupTor (r,t)
      | r <= 0 = "0"
      | r == 1 = "\\mathbb Z/" ++ show t
      | r >= 2 = "\\left(\\mathbb Z/" ++ show t ++ "\\right)^{" ++ show r ++ "}"
    torpart trs = intercalate "\\oplus " $ map markupTor (flatZip trs)

showHomologyTableTeX :: Map.Map (Int,Int) KHData -> String
showHomologyTableTeX khMap =
  let mayRange = foldl' rangeFinder Nothing (Map.keys khMap)
  in case mayRange of
       Nothing -> ""
       (Just (imin,imax,jmin,jmax)) -> flip execState "" $ do
         modify' (++("\\begin{array}{r|" ++ replicate (jmax-jmin+1) 'c' ++ "}\n"))
         modify' (++"i\\backslash j & ")
         forM_ [jmin..jmax] $ \j ->
           modify'(++(show j ++ if j<jmax then "&" else "\\\\\\hline\n"))
         forM_ [imin..imax] $ \i -> do
           modify' (++(show i ++ "&"))
           forM_ [jmin..jmax] $ \j -> do
             case khMap Map.!? (i,j) of
               Just khdata ->
                 modify' (++ showAbGroupTeX (rank khdata) (tors khdata))
               Nothing ->
                 return ()
             when (j < jmax) $ modify' (++"&")
           when (i < imax) $ modify' (++"\\\\\n")
         modify' (++ "\\end{array}")
  where
    rangeFinder Nothing (i,j) = Just (i,i,j,j)
    rangeFinder (Just (imin,imax,jmin,jmax)) (i,j) =
      let imin' = if i < imin then i else imin
          imax' = if i > imax then i else imax
          jmin' = if j < jmin then j else jmin
          jmax' = if j > jmax then j else jmax
      in Just (imin',imax',jmin',jmax')

showKHDataTikz :: Int -> Int -> KHData -> String
showKHDataTikz i j khData = flip execState "" $ do
  modify' (++ texParagraph "Homology group")
  modify' (++ "\\[\n")
  modify' (++("\\overline{\\mathit{Kh}}^{" ++ show i ++ "," ++ show j ++ "}\n"))
  modify' (++("\\cong" ++ showAbGroupTeX (rank khData) (tors khData)))
  modify' (++ "\\]\n")
  modify' (++ texParagraph "Generating cycles")
  forM_ (cycleV khData) $ \cyc -> modify' (++ showStateSumTikz 0.15 cyc)
  modify' (++ texParagraph "Killing boundaries")
  forM_ (bndryV khData) $ \bnd -> modify' (++ showStateSumTikz 0.15 bnd)

docKhovanovTikz :: ArcGraph -> String -> String -> Map.Map (Int,Int) KHData -> String
docKhovanovTikz ag opts cls khMap = flip execState "" $ do
  modify' (++ texHeader opts cls)
  modify' (++ texPreamble)
  modify' (++ texBegin)
  modify' (++"\n\\section{}\n")
  modify' (++"\\begin{center}\n")
  modify' (++tikzBegin 2.0)
  modify' (++showArcGraphTikz 0.15 ag)
  modify' (++tikzEnd)
  modify' (++"\\end{center}\n")
  modify' (++"\\section*{Table of homology groups}\n")
  modify' (++"\\[")
  modify' (++showHomologyTableTeX khMap)
  modify' (++"\\]")
  forM_ (Map.keys khMap) $ \ind -> do
    let (i,j) = ind
    modify' (++texHorizontalLine)
    modify' (++ showKHDataTikz i j (khMap Map.! (i,j)))
  modify' (++ texEnd)

typesetArcGraphTikz :: String -> String -> [ArcGraph] -> String
typesetArcGraphTikz option cls ags = flip execState "" $ do
  modify' (++texHeader option cls)
  modify' (++texPreamble)
  modify' (++texBegin)
  forM_ ags $ \ag -> do
    modify' (++tikzBegin 1.0)
    modify' (++showArcGraphTikz 0.15 ag)
    modify' (++tikzEnd)
  modify' (++texEnd)

docArcGraphTikz :: String -> String -> [ArcGraph] -> String
docArcGraphTikz option cls ags = flip execState "" $ do
  modify' (++texHeader option cls)
  modify' (++texPreamble)
  modify' (++texBegin)
  forM_ ags $ \ag -> do
    modify' (++"\n\\section{}\n")
    modify' (++"\\begin{center}\n")
    modify' (++tikzBegin 1.0)
    modify' (++showArcGraphTikz 0.15 ag)
    modify' (++tikzEnd)
    modify' (++"\\end{center}\n")
    forM_ [0..countCross ag] $ \i -> do
      modify' (++texHorizontalLine)
      modify' (++"\\subsection*{$H^{" ++ show i ++ ",\\star}$}\n")
      -- modify' (++texHorizontalLine)
      modify' (++"\\begin{center}\n")
      forM_ (listSmoothing [i] ag) $ \agsm -> do
        modify' (++tikzBegin 1.0)
        modify' (++showArcGraphTikz 0.15 agsm)
        modify' (++tikzEnd)
        modify' (++"\\hfill\n")
      modify' (++"\\end{center}\n")
  modify' (++texEnd)
