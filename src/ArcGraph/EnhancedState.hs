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

import ArcGraph

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
-- Step up the STRef counter and return the old value
stepST :: (Enum a) => STRef s a -> ST s a
stepST cntRef= do
  cnt <- readSTRef cntRef
  modifySTRef' cntRef succ
  return cnt

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
  indcs <- V.ifoldl' (\xs i y -> (i,fst y):xs) [] <$> V.unsafeFreeze edgeMV
  let justEqSnd (_,Just x) (_,Just y) = x==y
      justEqSnd _ _ = False
  return $ map (map fst) $ L.groupBy justEqSnd indcs

