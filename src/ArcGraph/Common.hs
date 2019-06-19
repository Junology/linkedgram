------------------------------------------------
-- |
-- Module    :  ArcGraph.Common
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Utility tools on drawing ArcGraph
--
------------------------------------------------

module ArcGraph.Common where

import Data.List

import ArcGraph

type BezierQuad = (Vertex, Vertex, Vertex)
type BezierCube = (Vertex, Vertex, Vertex, Vertex)

elevateBezier :: BezierQuad -> BezierCube
elevateBezier (v0@(x0,y0), (x1, y1), v2@(x2, y2))
  = let w1 = ((x0+x1*2.0)/3.0, (y0+y1*2.0)/3.0)
        w2 = ((x1*2.0+x2)/3.0, (y1*2.0+y2)/3.0)
    in (v0, w1, w2, v2)

dupEnds :: [a] -> [a]
dupEnds = dupHead . dupLast
  where dupHead [] = []
        dupHead (x:xs) = x:x:xs
        dupLast [] = []
        dupLast [x] = [x,x]
        dupLast (x:xs@(_:_)) = x:dupLast xs

swapEnds :: [a] -> [a]
swapEnds [] = []
swapEnds [x] = [x]
swapEnds (x:xs@(_:_)) = last xs : init xs ++ [x]

applyAdj3 :: (a -> a -> a -> b) -> [a] -> [b]
applyAdj3 _ [] = []
applyAdj3 _ (_:[]) = []
applyAdj3 f (x:y:ys) = snd $ mapAccumL upgrade (x,y) ys
  where
    upgrade (u,v) w = ((v,w), f u v w)

processBnd :: ArcPath -> [Vertex]
processBnd (APath OpenPath ps) = dupEnds ps
processBnd (APath ClosedPath ps) = (swapEnds . dupEnds) ps

mkBezier :: ArcPath -> [BezierQuad]
mkBezier pth
  = let vs = processBnd pth
        bs = applyAdj3 mkBezier' vs
    in bs
  where
    mkBezier' (x0,y0) (x1,y1) (x2,y2)
      = let x'0 = (x0+x1)/2.0
            y'0 = (y0+y1)/2.0
            x'2 = (x1+x2)/2.0
            y'2 = (y1+y2)/2.0
        in ((x'0,y'0), (x1,y1), (x'2,y'2))

normalize :: Double -> ArcGraph -> ArcGraph
normalize rad ag
  = let vstot = nub $ arcGraphVrtx ag
        c@(cx,cy) = barycenter vstot
        mayag' = moveVrtx vstot (negate cx) (negate cy) ag
        maxrad = maximum $ map (distance c) vstot
        vnormalizer = \x -> (rad * fst x/maxrad, rad * snd x/maxrad)
    in case mayag' of
         Just (AGraph ps' cs')
           -> AGraph (map (mapArcPath vnormalizer) ps') (map (mapCross (mapSgmt vnormalizer)) cs')
         Nothing
           -> AGraph [] []
