{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}

------------------------------------------------
-- |
-- Module    :  ArcGraph
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- The data structure of arcs and crossings of a diagram
--
------------------------------------------------

module ArcGraph where

import Control.Monad
import Control.Monad.Identity
import Control.Monad.ST

import Data.Maybe
import Data.List as L
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import Data.STRef
import Data.Foldable (for_)

{- For debug
import Debug.Trace
--}

-----------
-- Types --
-----------
data Orientation = Positive | Negative | Zero
  deriving (Show, Eq)

instance Semigroup Orientation where
  Positive <> Positive = Positive
  Positive <> Negative = Negative
  Negative <> Positive = Negative
  Negative <> Negative = Positive
  Zero <> _ = Zero
  _ <> Zero = Zero

type Vertex = (Double,Double)
data Segment = Sgmt Vertex Vertex
  deriving (Show, Read)

data PathType = OpenPath | ClosedPath
  deriving (Eq, Read, Show)

data ArcPath = APath PathType [Vertex]
  deriving (Eq,Show, Read)

data CrsState = Crossing | Smooth0 | Smooth1
  deriving (Eq,Show, Read)

instance Ord CrsState where
  compare Smooth0 Smooth0 = EQ
  compare Smooth0 Crossing = LT
  compare Smooth0 Smooth1 = LT
  compare Crossing Smooth0 = GT
  compare Crossing Crossing = EQ
  compare Crossing Smooth1 = LT
  compare Smooth1 Smooth0 = GT
  compare Smooth1 Crossing = GT
  compare Smooth1 Smooth1 = EQ

data Cross = Crs Segment Segment CrsState
  deriving (Show, Read)

-- | Segment is unordered
instance Eq Segment where
  (Sgmt v0 v1) == (Sgmt w0 w1)
    = (v0==w0 && v1==w1) || (v0==w1 && v1==w0)

-- | The order of segments in Cross is irrelevant
instance Eq Cross where
  (Crs sega segb _) == (Crs segc segd _)
    = (sega == segc && segb == segd) || (sega == segd && segb == segc)

-- | Cross are ordered in terms of crossing states
instance Ord Cross where
  compare (Crs _ _ crst) (Crs _ _ crst') = compare crst crst'

data ArcGraph = AGraph [ArcPath] [Cross]
  deriving (Eq,Show, Read)

-- | Order comparison on crossings
instance Ord ArcGraph where
  compare (AGraph _ cs) (AGraph _ cs') = compare cs cs'

-----------------------
-- Utility Functions --
-----------------------
applyAdj :: (a -> a -> b) -> [a] -> [b]
applyAdj _ [] = []
applyAdj _ [_] = []
applyAdj f (x:xs@(y:_)) = f x y : applyAdj f xs

applyCycAdj :: (a -> a -> b) -> [a] -> [b]
applyCycAdj _ [] = []
applyCycAdj _ [_] = [] -- Necessary case
applyCycAdj f (x:xs) = applyCycAdj' x xs
  where
    applyCycAdj' y [] = [f y x]
    applyCycAdj' y (z:zs) = f y z : applyCycAdj' z zs

mapLast :: (a -> a) -> [a] -> [a]
mapLast _ [] = []
mapLast f [x] = [f x]
mapLast f (x:xs) = x : mapLast f xs

groupWith :: (a -> a -> Maybe b) -> [a] -> ([b],[[a]])
groupWith _ [] = ([],[])
groupWith f (x:xs)
  = let (ws, yss, ys,z) = foldl wheel ([],[],[],x) xs
    in (ws, yss ++ [ys++[z]])
  where
    --wheel :: ([b],[[a]],[a],a) -> a -> ([b],[[a]],[a],a)
    wheel (ws,xss,ys,y) z
      = case f y z of
          Just u -> (ws ++ [u], xss++[ys++[y]], [], z)
          Nothing -> (ws, xss, ys++[y], z)

groupCycWith :: (a -> a -> Maybe b) -> [a] -> ([b],[[a]])
groupCycWith _ [] = ([],[])
groupCycWith f xs@(x:_) =
  case groupWith f (xs ++ [x]) of
    (mt, []) -> (mt, [])
    (mt, []:wss) -> (mt,wss)
    (mt, (y:ys):wss) -> (mt, mapLast (++ys) wss)

condMap :: Functor f => (a -> Bool) -> (a -> a) -> f a -> f a
condMap p f = fmap (\x -> if p x then f x else x)

mapHead :: (a -> a) -> [a] -> [a]
mapHead _ [] = []
mapHead f (x:xs) = (f x):xs

mapAt :: Int -> (a -> a) -> [a] -> [a]
mapAt n f xs
  = let (ls,rs) = splitAt n xs
    in ls ++ mapHead f rs

distance :: Vertex -> Vertex -> Double
distance (a0,a1) (b0,b1)
  = let (v0,v1) = (b0-a0,b1-a1)
    in sqrt (v0*v0 + v1*v1)

mapRemoveM ::((a -> Identity a) -> b -> Identity b) -> (a -> a) -> b -> b
mapRemoveM mAp f = runIdentity . mAp (return . f)

--------------
-- Vertices --
--------------
-- recognize if three points are in counter-clockwise, in clockwise, or on a line.
orientation :: Vertex -> Vertex -> Vertex -> Orientation
orientation (ax,ay) (bx,by) (cx,cy)
  = case (bx-ax)*(cy-ay) - (cx-ax)*(by-ay) of
      det
        | det < 0  -> Negative
        | det > 0  -> Positive
        | det == 0 -> Zero

barycenter :: Foldable f => f Vertex -> Vertex
barycenter vs = let l = length vs in foldl' (mapSum l) (0.0,0.0) vs
  where
    mapSum l (x0,y0) (x1,y1) = (x0+(x1/fromIntegral l), y0+(y1/fromIntegral l))

--------------
-- Segments --
--------------
tposeSegment :: Segment -> Segment
tposeSegment (Sgmt v w) = Sgmt w v

mapSgmtM :: (Monad m) => (Vertex -> m Vertex) -> Segment -> m Segment
mapSgmtM f (Sgmt v w) = do {v' <- f v; w' <- f w; return (Sgmt v' w')}

mapSgmt :: (Vertex -> Vertex) -> Segment -> Segment
mapSgmt = mapRemoveM mapSgmtM

--------------
-- Crossing --
--------------
calcCross :: Segment -> Segment -> Maybe Vertex
calcCross (Sgmt (x0,y0) (x1,y1)) (Sgmt (w0,z0) (w1,z1)) = do
  let (dx,dy) = (x1-x0, y1-y0)
      (dw,dz) = (w1-w0, z1-z0)
      det = dy*dw - dx*dz
  when (det == 0.0) $ fail ""
  let t = ((x0-w0)*dz - (y0-z0)*dw) / det
      s = ((x0-w0)*dy - (y0-z0)*dx) / det
  when (t <= 0.0 || t >= 1.0 || s <= 0.0 || s >= 1.0) $ fail ""
  let (cx,cy) = (x0 + dx*t, y0 + dy*t)
  return (cx,cy)

crossing :: CrsState -> Segment -> Segment -> Maybe Cross
crossing crst sega segb = do
  calcCross sega segb
  return $ Crs sega segb crst

crsSgmt :: Segment -> Vertex -> Vertex -> Maybe Segment
crsSgmt seg v w = do
  let seg' = Sgmt v w
  calcCross seg seg'
  return seg'

crsVrtx :: Cross -> [Vertex]
crsVrtx (Crs (Sgmt v0 v1) (Sgmt w0 w1) _) = [v0,v1,w0,w1]

positate :: Cross -> Cross
positate c@(Crs sega@(Sgmt v0 v1) (Sgmt w0 w1) crst)
  = case orientation v0 v1 w1 of
      Negative -> Crs sega (Sgmt w1 w0) crst
      _ -> c

crsChange :: Cross -> Cross
crsChange (Crs sega segb crst)
  = Crs (tposeSegment segb) sega crst

crsModifyState :: (CrsState -> CrsState) -> Cross -> Cross
crsModifyState f (Crs sega segb crst) = Crs sega segb (f crst)

crsSetState :: CrsState -> Cross -> Cross
crsSetState = crsModifyState . const

crsRotateState :: Cross -> Cross
crsRotateState = crsModifyState $ \st ->
  case st of
    Crossing -> Smooth0
    Smooth0 -> Smooth1
    Smooth1 -> Crossing

mapCrossM :: (Monad m) => (Segment -> m Segment) -> Cross -> m Cross
mapCrossM f (Crs sega segb crst)
  = do {sega'<- f sega; segb' <- f segb; return (Crs sega' segb' crst)}

mapCross :: (Segment -> Segment) -> Cross -> Cross
mapCross = mapRemoveM mapCrossM

hasCrossing :: Cross -> Bool
hasCrossing (Crs sega segb _)
  = isJust $ calcCross sega segb

---------------------
-- Path Operations --
---------------------
arcLength :: ArcPath -> Int
arcLength (APath OpenPath vs) = length vs - 1
arcLength (APath ClosedPath []) = 0
arcLength (APath ClosedPath vs@(_:_)) = length vs

hasEdge :: ArcPath -> Bool
hasEdge p = arcLength p >= 1

arcSegments :: ArcPath -> [Segment]
arcSegments (APath OpenPath vs) = applyAdj Sgmt vs
arcSegments (APath ClosedPath vs) = applyCycAdj Sgmt vs

arcVertices :: ArcPath -> [Vertex]
arcVertices (APath _ vs) = vs

cutArc :: CrsState -> Segment -> ArcPath -> ([ArcPath],[Cross])
cutArc crst seg (APath OpenPath vs) =
  let (cs,vss) = groupWith (\v w -> crossing crst seg (Sgmt v w)) vs
  in (filter hasEdge (map (APath OpenPath) vss), map positate cs)
cutArc _ _ (APath ClosedPath []) = ([],[])
cutArc crst seg p@(APath ClosedPath vs) =
  let (cs,vs') = groupCycWith (\v w -> crossing crst seg (Sgmt v w)) vs
  in if null cs
     then ([p],[])
     else (filter hasEdge (fmap (APath OpenPath) vs'), map positate cs)

cutArcs :: CrsState -> Segment -> [ArcPath] -> ([ArcPath],[Cross])
cutArcs crst seg ps
  = let (pss, css) = unzip $ map (cutArc crst seg) ps
    in (concat pss, concat css)

cutArcsByArcs :: CrsState -> [ArcPath] -> ([ArcPath],[Cross])
cutArcsByArcs crst ps
  = let segs = concatMap arcSegments ps
        (qs,css) = mapAccumL (flip (cutArcs crst)) ps segs
        --(css,pss) = unzip $ map (flip (cutArcs crst) ps) segs
    in (qs,nub (concat css))

mapArcPathM :: (Monad m) => (Vertex -> m Vertex) -> ArcPath -> m ArcPath
mapArcPathM f (APath ty vs)
  = do {vs' <- mapM f vs; return (APath ty vs')}

mapArcPath :: (Vertex -> Vertex) -> ArcPath -> ArcPath
mapArcPath = mapRemoveM mapArcPathM

getEndVrtx :: ArcPath -> [Vertex]
getEndVrtx (APath _ []) = []
getEndVrtx (APath ClosedPath (x:_)) = [x]
getEndVrtx (APath OpenPath [x]) = [x]
getEndVrtx (APath OpenPath (x:xs@(_:_))) = [x,last xs]

--------------
-- ArcGraph --
--------------
mkArcGraph :: CrsState -> [ArcPath] -> ArcGraph
mkArcGraph crst = uncurry AGraph . cutArcsByArcs crst

insertArcs :: CrsState -> [ArcPath] -> ArcGraph -> ArcGraph
insertArcs crst ps (AGraph qs cs)
  = let (rs,ds) = cutArcsByArcs crst (ps++qs)
    in AGraph rs (ds++cs)

sumGraph :: CrsState -> ArcGraph -> ArcGraph -> ArcGraph
sumGraph crst (AGraph ps cs) (AGraph qs ds)
  = let (rs,es) = cutArcsByArcs crst (ps++qs)
    in AGraph rs (cs++ds++es)

moveVrtx :: [Vertex] -> Double -> Double -> ArcGraph -> Maybe ArcGraph
moveVrtx vs x y (AGraph ps cs) = do
  let ps' = map updPath ps
  cs' <- mapM updCrs cs
  return $ AGraph ps' cs'
  where
    updVrtx :: Vertex -> Vertex
    updVrtx (x0,y0) = (x0+x,y0+y)
    updPath :: ArcPath -> ArcPath
    updPath (APath ty ws)
      = APath ty (condMap (`elem` vs) updVrtx ws)
    updSeg :: Segment -> Segment
    updSeg (Sgmt a b)
      = let a' = if a `elem` vs then updVrtx a else a
            b' = if b `elem` vs then updVrtx b else b
        in Sgmt a' b'
    updCrs :: Cross -> Maybe Cross
    updCrs (Crs sega segb st) = do
      let sega' = updSeg sega
          segb' = updSeg segb
      calcCross sega' segb'
      return $ Crs sega' segb' st

removeVrtx :: Vertex -> ArcGraph -> ArcGraph
removeVrtx v (AGraph ps cs)
  = let remover :: Vertex -> Maybe Vertex
        (remover,ps'') = mapAccumL rmVrtxFromPath return ps
        ps' = filter hasEdge ps''
        cs' = filter hasCrossing $ mapMaybe (mapCrossM (mapSgmtM remover)) cs
    in AGraph ps' cs'
  where
    rmVrtxFromPath :: (Vertex -> Maybe Vertex) -> ArcPath -> (Vertex -> Maybe Vertex,ArcPath)
    rmVrtxFromPath f p@(APath OpenPath vs)
      = case elemIndex v vs of
          (Just i)
            | length vs >= 3 && i==0
              -> let f' = f >=> \v'-> if v==v' then Just (vs!!1) else Just v'
                 in (f', APath OpenPath (tail vs))
            | length vs >= 3 && i==length vs-1
              -> let f' = f >=> \v'-> if v==v' then Just (vs!!(i-1)) else Just v'
                 in(f', APath OpenPath (init vs))
            | otherwise
              -> (f,APath OpenPath (delete v vs))
          Nothing -> (f,p)
    rmVrtxFromPath f p@(APath ClosedPath vs)
      = case elemIndex v vs of
          (Just i) -> (f,APath ClosedPath (delete v vs))
          Nothing -> (f,p)

arcGraphVrtx :: ArcGraph -> [Vertex]
arcGraphVrtx (AGraph ps cs) = concatMap arcVertices ps ++ concatMap crsVrtx cs

modifyCross :: Int -> (Cross -> Cross) -> ArcGraph -> ArcGraph
modifyCross n f (AGraph ps cs)
  = AGraph ps (mapAt n f cs)

setCrsState :: Int -> CrsState -> ArcGraph -> ArcGraph
setCrsState n crst = modifyCross n (crsSetState crst)

countCross :: ArcGraph -> Int
countCross (AGraph _ cs) = length cs

tryGrowArcM :: (Monad m) => Vertex -> Vertex -> ArcPath -> m ArcPath
tryGrowArcM _ _ (APath ClosedPath _) = fail ""
tryGrowArcM _ _ (APath OpenPath []) = fail ""
tryGrowArcM vrt vnew (APath OpenPath vs@(v:_))
  | vrt == v       = return (APath OpenPath (vnew:vs))
  | vrt == last vs = return (APath OpenPath (vs++[vnew]))
  | otherwise      = fail ""

findIndexWith :: (a -> Maybe b) -> V.Vector a -> Maybe (b,Int)
findIndexWith f = V.ifoldl' bin Nothing
  where
    bin Nothing i x = case f x of {Just y -> Just (y,i); Nothing -> Nothing;}
    bin j@(Just _) _ _  = j

slimCross :: ArcGraph -> ArcGraph
slimCross (AGraph ps cs) = runST $ do
  psMVecRef <- newSTRef =<< V.thaw (V.fromList ps)
  csMVec <- V.thaw $ V.fromList cs
  for_  [0..(length cs-1)] $ \i -> do
    let (Crs sega@(Sgmt v0 v1) segb@(Sgmt w0 w1) crst) = cs!!i
        c = fromJust $ calcCross sega segb
        vs'@(v0':v1':w0':w1':_) = map (\v -> barycenter [v,c]) [v0,v1,w0,w1]
    for_ [0..3] $ \k -> do
      let v = [v0,v1,w0,w1]!!k
          v'= vs'!!k
      psVec <- V.freeze =<< readSTRef psMVecRef
      case findIndexWith (tryGrowArcM v v') psVec of
        Just (p,j) -> do
          psMVec <- readSTRef psMVecRef
          MV.write psMVec j p
        Nothing -> do
          writeSTRef psMVecRef =<< flip MV.unsafeGrow 1 =<< readSTRef psMVecRef
          psMVec <- readSTRef psMVecRef
          MV.write psMVec (MV.length psMVec-1) (APath OpenPath [v,v'])
    MV.write csMVec i $ Crs (Sgmt v0' v1') (Sgmt w0' w1') crst
  ps' <- V.toList <$> (V.freeze =<< readSTRef psMVecRef)
  cs' <- V.toList <$> V.freeze csMVec
  return (AGraph ps' cs')

-------------
-- HitTest --
-------------
hitVertex :: Double -> Vertex -> Vertex -> Bool
hitVertex r v w = distance v w < r

hitCross :: Vertex -> Cross -> Bool
hitCross v (Crs (Sgmt a0 a1) (Sgmt b0 b1) _)
  = orientation a0 b0 v <> orientation a0 b0 a1 == Positive
    && orientation a0 b1 v <> orientation a0 b1 a1 == Positive
    && orientation a1 b0 v <> orientation a1 b0 a0 == Positive
    && orientation a1 b1 v <> orientation a1 b1 a0 == Positive

hitTest :: Double -> Vertex -> ArcGraph -> (Maybe Vertex, Maybe Int)
hitTest rad v ag@(AGraph _ cs)
  = let mv = find (hitVertex rad v) (arcGraphVrtx ag)
        mi = findIndex (hitCross v) cs
    in (mv,mi)
