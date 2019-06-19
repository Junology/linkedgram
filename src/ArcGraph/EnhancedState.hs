{-# LANGUAGE Strict #-}
{-# LANGUAGE StrictData #-}

------------------------------------------------
-- |
-- Module    :  ArcGraph.EnhancedState
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Implementing manipulations on enhanced states
--
------------------------------------------------

module ArcGraph.EnhancedState where

import Control.Applicative
import Control.Monad
import Control.Monad.ST

import Data.STRef
import Data.Maybe as M
import Data.Foldable
import qualified Data.List as L
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import qualified Data.Map.Strict as Map

import ArcGraph

import Numeric.Algebra.FreeModule as FM
import Numeric.Algebra.Frobenius as Frob
import Numeric.Algebra.IntMatrix

---------------
-- Utilities --
---------------
-- Step up the STRef counter and return the old value
stepST :: (Enum a) => STRef s a -> ST s a
stepST cntRef= do
  cnt <- readSTRef cntRef
  modifySTRef' cntRef succ
  return cnt

-----------------------
-- Unenhanced states --
-----------------------
type DiagramState = [Cross]

smoothing :: DiagramState -> ArcGraph -> ArcGraph
smoothing st (AGraph ps cs)
  = AGraph ps $ map (\c -> if elem c st then crsSetState Smooth1 c else crsSetState Smooth0 c) cs

listSmoothing :: [Int] -> ArcGraph -> [ArcGraph]
listSmoothing degrees ag@(AGraph _ cs)
  = let states = filter ((`elem` degrees) . length) $ L.subsequences cs
    in map (flip smoothing ag) states

localCross :: Cross -> (Segment,Segment)
localCross (Crs sega segb Crossing) = (sega,segb)
localCross (Crs (Sgmt v0 v1) (Sgmt w0 w1) Smooth0) = (Sgmt v0 w1, Sgmt v1 w0)
localCross (Crs (Sgmt v0 v1) (Sgmt w0 w1) Smooth1) = (Sgmt v0 w0, Sgmt v1 w1)

---------------------
-- Enhanced states --
---------------------
mkConnection :: (a -> Bool) -> MV.STVector s (Maybe Int,a) -> STRef s Int -> ST s ()
mkConnection p edgeMV cntRef = do
  edgeV <- V.freeze edgeMV
  cnt <- stepST cntRef
  oldLabRef <- newSTRef []
  for_ (V.findIndices (p . snd) edgeV) $ \j -> do
    case fst (edgeV V.! j) of
      Just k -> modifySTRef' oldLabRef (k:)
      Nothing -> MV.modify edgeMV (\x -> (Just cnt,snd x)) j
  oldLab <- readSTRef oldLabRef
  let isOldLab (mayl,_) = M.fromMaybe False $ (`elem` oldLab) <$> mayl
  forM_ (V.findIndices isOldLab edgeV) $ \j -> do
    MV.modify edgeMV (\x -> (Just cnt,snd x)) j

components :: ArcGraph -> [[Int]]
components ag@(AGraph ps cs) = runST $ do
  let n = length ps
  edgeMV <- MV.new n
  -- Initialize buffer
  forM_ [0..(n-1)] $ \i -> do
    MV.unsafeWrite edgeMV i (Nothing, getEndVrtx (ps!!i))
  cntRef <- newSTRef 0
  forM_ cs $ \c -> do
    let (Sgmt v0 v1, Sgmt w0 w1) = localCross c
    -- Connect v0 <--> v1
    mkConnection (\x -> elem v0 x || elem v1 x) edgeMV cntRef
    -- Connect w0 <--> w1
    mkConnection (\x -> elem w0 x || elem w1 x) edgeMV cntRef
  (justcls,nocls) <- V.ifoldl' (\xs i y -> case fst y of {Just k -> ((k,[i]):fst xs,snd xs); Nothing -> (fst xs,[i]:snd xs);}) ([],[]) <$> V.unsafeFreeze edgeMV
  return $ (Map.elems $ Map.fromListWith (++) justcls) ++ nocls

data ArcGraphE = AGraphE ArcGraph (Map.Map [Int] SL2B)
  deriving (Eq,Show)

-- | Lexicographical order on ArcGraphE
instance Ord ArcGraphE where
  compare (AGraphE ag coeff) (AGraphE ag' coeff')
    = case compare ag ag' of
        LT -> LT
        GT -> GT
        EQ -> compare coeff coeff'

enhancements :: ArcGraph -> V.Vector (Map.Map [Int] SL2B)
enhancements ag = runST $ do
  let comps = components ag
      subis = zip [0..] $ L.subsequences comps
  mapMV <- MV.unsafeNew (length subis) -- Do not need to initialize the memory
  for_ subis $ \subi -> do
    let (i,sub) = subi
    MV.write mapMV i $ L.foldl' (\xs y -> Map.insert y (if elem y sub then SLI else SLX) xs) Map.empty comps
  V.unsafeFreeze mapMV

enhancedStates :: ArcGraph -> V.Vector ArcGraphE
enhancedStates ag = V.map (AGraphE ag) $ enhancements ag

enhancedStatesL :: ArcGraph -> [ArcGraphE]
enhancedStatesL ag = V.foldr' (\x xs -> (AGraphE ag x):xs) [] $ enhancements ag

-- | Compute next states with sign in differential.
-- On determining signs, we use the descending order on corssings.
diffState :: ArcGraph -> V.Vector (Int,ArcGraph)
diffState ag@(AGraph _ cs) = runST $ do
  resVecRef <- newSTRef V.empty
  signRef <- newSTRef 1
  for_ [0..(length cs - 1)] $ \i -> do
    case cs!!i of
      (Crs _ _ Smooth0) -> do
        sign <- readSTRef signRef
        modifySTRef' resVecRef $ flip V.snoc (sign,setCrsState i Smooth1 ag)
      (Crs _ _ Smooth1) -> do
        modifySTRef' signRef (*(-1))
      (Crs _ _ Crossing) -> do
        return ()
  readSTRef resVecRef

hasIntersection :: Eq a => [a] -> [a] -> Bool
hasIntersection x y = not $ L.null (x `L.intersect` y)

differential :: ArcGraphE -> FreeMod Int ArcGraphE
differential (AGraphE ag coeffMap) =
  let dStVec = diffState ag
      comps = components ag
  in sumFM $ runST $ do
    imageMV <- MV.unsafeNew (V.length dStVec)
    for_ [0..(V.length dStVec -1)] $ \i -> do
      let (sign,dAg) = (dStVec V.! i)
          imageAg = AGraphE dAg FM.@$>% Map.fromList FM.@$>% Frob.tqftZ hasIntersection (Map.toList coeffMap) (components dAg)
      MV.write imageMV i (sign FM.@*% imageAg)
    V.unsafeFreeze imageMV
