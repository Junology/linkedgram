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

import Control.Monad
import Control.Monad.ST

import Data.STRef
import Data.Maybe as M
import Data.Foldable
import qualified Data.List as L
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import qualified Data.Map.Strict as Map

import qualified Numeric.LinearAlgebra as LA

import ArcGraph

import Numeric.Algebra.FreeModule as FM
import Numeric.Algebra.Frobenius as Frob
import Numeric.Algebra.IntMatrix

{-- for debug
import Debug.Trace
--}

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
type DiagramState = [Int]

smoothing :: DiagramState -> ArcGraph -> ArcGraph
smoothing st (AGraph ps cs)
  = AGraph ps $ map mkCrs $ zip [0..] cs
  where mkCrs (i,Crs sega segb _) = if i `elem` st
                                    then Crs sega segb Smooth1
                                    else Crs sega segb Smooth0

listStates :: ArcGraph -> Int -> [DiagramState]
listStates (AGraph _ cs) deg =
  filter ((==deg) . length) $ L.subsequences [0..(length cs - 1)]

listSmoothing :: [Int] -> ArcGraph -> [ArcGraph]
listSmoothing degrees ag
  = map (flip smoothing ag) $ concatMap (listStates ag) degrees

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
components (AGraph ps cs) = runST $ do
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

data ArcGraphE = AGraphE {
  arcGraph :: ArcGraph,
  state :: DiagramState,
  enhLabel :: (Map.Map [Int] SL2B)
  } deriving (Eq,Show)

-- | Lexicographical order on ArcGraphE
instance Ord ArcGraphE where
  compare (AGraphE _ st coeff) (AGraphE _ st' coeff')
    = case compare st st' of
        LT -> LT
        GT -> GT
        EQ -> compare coeff coeff'

enhancements :: Int -> ArcGraph -> DiagramState -> V.Vector (Map.Map [Int] SL2B)
enhancements deg ag st = runST $ do
  let comps = components (smoothing st ag)
      deg' = deg + L.length comps
      subis = zip [0..] $ filter (\sub -> 2*L.length sub == deg') $ L.subsequences comps
  mapMV <- MV.unsafeNew (length subis) -- Do not need to initialize the memory
  for_ subis $ \subi -> do
    let (i,sub) = subi
    MV.write mapMV i $ L.foldl' (\xs y -> Map.insert y (if elem y sub then SLI else SLX) xs) Map.empty comps
  V.unsafeFreeze mapMV

enhancedStates :: Int -> ArcGraph -> DiagramState -> V.Vector ArcGraphE
enhancedStates deg ag st
  = V.map (AGraphE ag st) $ enhancements deg ag st

enhancedStatesL :: Int -> ArcGraph -> DiagramState -> [ArcGraphE]
enhancedStatesL deg ag st
  = V.foldr' (\x xs -> (AGraphE ag st x):xs) [] $ enhancements deg ag st

-- | Compute next states with sign in differential.
-- On determining signs, we use the descending order on corssings.
diffState :: ArcGraph -> DiagramState -> V.Vector (Int,DiagramState)
diffState (AGraph _ cs) st = runST $ do
  resVecRef <- newSTRef V.empty
  signRef <- newSTRef 1
  for_ [0..(length cs - 1)] $ \i -> do
    if i `elem` st
      then modifySTRef' signRef negate
      else do
        sign <- readSTRef signRef
        modifySTRef' resVecRef $ flip V.snoc (sign, L.sort (i:st))
  readSTRef resVecRef

hasIntersection :: Eq a => [a] -> [a] -> Bool
hasIntersection x y = not $ L.null (x `L.intersect` y)

differential :: ArcGraphE -> FreeMod Int ArcGraphE
differential (AGraphE ag st coeffMap) =
  let dStVec = diffState ag st
  in FM.sumFM $ runST $ do
    imageMV <- MV.unsafeNew (V.length dStVec)
    for_ [0..(V.length dStVec -1)] $ \i -> do
      let (sign,dSt) = (dStVec V.! i)
          imageAg = AGraphE ag dSt FM.@$>% Map.fromList FM.@$>% Frob.tqftZ hasIntersection (Map.toList coeffMap) (components (smoothing dSt ag))
      MV.write imageMV i (sign FM.@*% imageAg)
    V.unsafeFreeze imageMV

----------------------------------------------------------
-- The computation of (unnormalized) Khovanov homology
-----------------------------------------------------------
-- | The type to carry the data of Khovanov homologies
data KHData = KHData {
  rank :: Int,
  tors :: [Int],
  cycleV :: [FreeMod Int ArcGraphE],
  bndryV :: [FreeMod Int ArcGraphE] }

-- | Compute Khovanov homology for given range of cohomological degrees and a given quantum-degree
computeKhovanov :: ArcGraph -> [Int] -> Int -> [DiagramState] -> Map.Map (Int,Int) KHData
computeKhovanov ag hdegs qdeg states =
  let numCrs = countCross ag
      hdegsNHD = filter isNhd [0..numCrs]
        where isNhd i = elem i hdegs || elem (i+1) hdegs || elem (i-1) hdegs
      slimAG = slimCross ag
  in runST $ do
    -- Define STRef to write the result
    resultMapRef <- newSTRef (Map.empty)
    -- Compute basis
    baseMVec <- MV.replicate (numCrs + 2) []
    forM_ [0..numCrs] $ \i -> do
      -- Get selected states at coh.degree i
      let statesi = filter ((==i). L.length) states
      MV.write baseMVec i $ L.concatMap (enhancedStatesL (qdeg-i) slimAG) statesi
    -- Compute kernels and images of differentials
    cycleMV <- MV.replicate (numCrs+1) ([] :: [LA.Vector LA.Z])
    bndryMV <- MV.replicate (numCrs+1) ([] :: [LA.Vector LA.Z])
    diagcMV <- MV.replicate (numCrs+1) ([] :: [Int])
    forM_ hdegsNHD $ \i -> do
      base0 <- MV.read baseMVec i
      base1 <- MV.read baseMVec (i+1)
      if L.null base1
        then MV.write cycleMV i $ LA.toColumns (LA.ident (L.length base0))
        else unless (L.null base0) $ do
          let diffMat = genMatrix differential base0 base1
          let (dvec,ker,im) = kerImOf $ matiDataToLA diffMat
          MV.write cycleMV i $ ker
          MV.write bndryMV (i+1) $ im
          MV.write diagcMV (i+1) $ map fromIntegral $ LA.toList dvec
    -- Compute Homology at degree (i,qdeg)
    forM_ hdegs $ \i -> do
      when (i <= numCrs) $ do
        base <- MV.read baseMVec i
        cycleL <- MV.read cycleMV i
        bndryL <- MV.read bndryMV i
        let freeRk = (length cycleL) - (length bndryL)
        tor <- L.filter (/=1) <$> MV.read diagcMV i
        when (freeRk > 0 || not (L.null tor)) $ do
          let cycs' = map (map fromIntegral . LA.toList) cycleL
          let bnds' = map (map fromIntegral . LA.toList) bndryL
          let (cycs'',bnds'') = (cycs' L.\\ bnds', bnds' L.\\ cycs')
          modifySTRef' resultMapRef $ Map.insert (i,qdeg) $
            KHData freeRk tor (map (flip zipSum base) cycs'') (map (flip zipSum base) bnds'')
    readSTRef resultMapRef
