------------------------------------------------
-- |
-- Module    :  ArcGraph.Cairo
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Cairo output of ArcGraph
--
------------------------------------------------

module ArcGraph.Cairo (
  drawVertex,
  drawArcPath,
  drawArcVertex,
  drawArcGraph,
  drawArcGraphVrtx
  ) where

import Control.Monad
import Control.Monad.ST

import qualified Data.Map.Strict as Map
import Data.Maybe
import Data.List
import Data.STRef

import qualified Graphics.Rendering.Cairo as Cairo

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius

import ArcGraph
import ArcGraph.Common
import ArcGraph.EnhancedState

drawVertex :: Double -> Double -> Double -> Double -> Vertex -> Cairo.Render ()
drawVertex r g b rad (x,y) = do
  Cairo.save
  Cairo.newPath
  Cairo.setLineWidth 3.0
  Cairo.arc x y rad 0.0 (2.0*pi)
  Cairo.setSourceRGB 1.0 1.0 1.0
  Cairo.fillPreserve
  Cairo.setSourceRGB r g b
  Cairo.strokePreserve
  Cairo.restore

drawArcPath :: ArcPath -> Cairo.Render ()
drawArcPath pth =
  case fmap elevateBezier (mkBezier pth) of
    [] -> return ()
    bs@(((x0,y0),_,_,_):_) -> do
      Cairo.newPath
      Cairo.moveTo x0 y0
      forM_ bs $ \b -> do
        let (_, (x1,y1), (x2,y2), (x3,y3)) = b
        Cairo.curveTo x1 y1 x2 y2 x3 y3
      Cairo.stroke

drawArcVertex :: Double -> Double -> Double -> ArcPath -> Cairo.Render ()
drawArcVertex r g b (APath _ vs)
  = forM_ vs (drawVertex r g b 5.0)

drawCross :: Double -> Cross -> Cairo.Render ()
drawCross bd (Crs sega segb crst)
  = let (Sgmt (v0x,v0y) (v1x,v1y)) = sega
        (Sgmt (w0x,w0y) (w1x,w1y)) = segb
        (cx,cy) = fromJust $ calcCross sega segb
    in case crst of
         Crossing -> do
           Cairo.moveTo v0x v0y
           Cairo.lineTo v1x v1y
           Cairo.stroke
           Cairo.moveTo w0x w0y
           let d0 = distance (w0x,w0y) (cx,cy)
           Cairo.lineTo (cx + bd*(w0x-cx)/d0) (cy + bd*(w0y-cy)/d0)
           let d1 = distance (w1x,w1y) (cx,cy)
           Cairo.moveTo (cx + bd*(w1x-cx)/d1) (cy + bd*(w1y-cy)/d1)
           Cairo.lineTo w1x w1y
           Cairo.stroke
         Smooth0 -> do
           let (p1x,p1y) = ((v0x+2.0*cx)/3.0, (v0y+2.0*cy)/3.0)
               (p2x,p2y) = ((w1x+2.0*cx)/3.0, (w1y+2.0*cy)/3.0)
               (q1x,q1y) = ((w0x+2.0*cx)/3.0, (w0y+2.0*cy)/3.0)
               (q2x,q2y) = ((v1x+2.0*cx)/3.0, (v1y+2.0*cy)/3.0)
           Cairo.moveTo v0x v0y
           Cairo.curveTo p1x p1y p2x p2y w1x w1y
           Cairo.moveTo w0x w0y
           Cairo.curveTo q1x q1y q2x q2y v1x v1y
           Cairo.stroke
         Smooth1 -> do
           drawCross bd (Crs (tposeSegment sega) segb Smooth0)

drawArcGraphVrtx :: Double -> Double -> Double -> ArcGraph -> Cairo.Render ()
drawArcGraphVrtx r g b ag = do
  forM_ (arcGraphVrtx ag) (drawVertex r g b 5.0)

drawArcGraph :: ArcGraph -> Maybe (Double,Double,Double) -> Cairo.Render ()
drawArcGraph ag@(AGraph ps cs) mayrgb = do
  forM_ ps drawArcPath
  forM_ cs (drawCross 10.0)
  case mayrgb of
    Just (r,g,b) -> drawArcGraphVrtx r g b ag
    Nothing -> return ()

drawArcGraphEnh :: ArcGraphE -> Double -> Double -> Double -> Double -> Cairo.Render ()
drawArcGraphEnh (AGraphE ag@(AGraph ps _) coeffMap) sz r g b = do
  -- Draw arc graph
  drawArcGraph ag Nothing
  -- Draw labels in the enhanced state
  Cairo.save
  Cairo.setSourceRGB r g b
  Cairo.setFontSize sz
  forM_ (Map.toList coeffMap) $ \kv -> do
    let (comp,coeff) = kv
        (x,y) = fromMaybe (0,0) $ findMostVrtx (\v w -> fst v < fst w) $ map (ps!!) comp
    -- Draw label
    Cairo.moveTo x y
    Cairo.showText $ case coeff of {SLI -> "1"; SLX -> "X";}
  Cairo.restore

drawStateSum :: FreeMod Int ArcGraphE -> Double -> Double -> Double -> Cairo.Render ()
drawStateSum vect r g b = do
  Cairo.save
  Cairo.setFontSize 20
  forEachWithInterM drawAG drawPlus vect
  Cairo.restore
  where
    drawAG coeff agE@(AGraphE ag _) = do
      Cairo.moveTo 0 0
      Cairo.showText (show coeff)
      Cairo.save
      Cairo.moveTo 20 0
      Cairo.scale 0 (-1)
      drawArcGraphEnh (agE {state = normalize 80 ag}) 20 r g b
      Cairo.restore
      Cairo.translate 100 0
    drawPlus = do
      Cairo.showText "+"
      Cairo.translate 20 0
