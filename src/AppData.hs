------------------------------------------------
-- |
-- Module    :  AppData
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Defining the main data structure that the application is carrying.
--
------------------------------------------------

module AppData where

import Data.List

import ArcGraph

data MoveMode = MoveVertex [Vertex] [Int] | MoveOrigin
  deriving Show

data AppMode
  = ModeView [Vertex] [Int]
  | ModeMove MoveMode (Double,Double)
  | ModeDraw ArcPath
  deriving Show

data AppData = AppData {
  appMode :: AppMode,
  origX :: Double,
  origY :: Double,
  arcGraph :: ArcGraph,
  fileName :: String,
  changed :: Bool }
  deriving Show

initialAppData :: AppData
initialAppData
  = AppData (ModeView [] []) 0.0 0.0 (AGraph [] []) "" True

openAppData :: String -> Double -> Double -> ArcGraph -> AppData
openAppData fname ox oy ag
  = let (cx,cy) = barycenter $ nub $ arcGraphVrtx ag
    in AppData (ModeView [] []) (ox-cx) (oy-cy) ag fname False

moveOrigin :: Double -> Double -> AppData -> AppData
moveOrigin ox oy appData = appData {origX = ox, origY = oy}

confirmChange :: AppData -> AppData
confirmChange appData = appData {changed = False}

setFileName :: String -> AppData -> AppData
setFileName name appData = appData {fileName = name}

hitRadius :: Double
hitRadius = 5.0

changeMode :: AppMode -> AppData -> AppData
changeMode mode appData = appData {appMode = mode}

mouseDown :: Double -> Double -> Bool -> AppData -> AppData
mouseDown mx my modif appData@(AppData mode ox oy graph fname ch)
  = case mode of
      ModeView vs is
        -> let (vs',is') = if modif then (vs,is) else ([],[])
           in case hitTest hitRadius (mx-ox,my-oy) graph of
                (Just v,_) -> appData {appMode = ModeView (v:vs') is'}
                (Nothing,Just i) -> appData {appMode = ModeView vs' (i:is')}
                (Nothing,Nothing) -> appData {appMode = ModeView vs' is'}
      ModeMove _ _
        -> appData
      ModeDraw (APath ty vs)
        -> appData {appMode = ModeDraw (APath ty (vs++[(mx-ox,my-oy)]))}

mouseDouble :: Double -> Double -> AppData -> AppData
mouseDouble mx my appData@(AppData mode ox oy graph fname _)
  = case mode of
      ModeView _ _
        -> case hitTest hitRadius (mx-ox,my-oy) graph of
             (Nothing,Just i)
               -> appData {arcGraph = modifyCross i (crsChange) graph,
                           changed = True}
             _ {-- Anything Else --}
               -> appData
      ModeMove _ _
        -> appData
      ModeDraw (APath ty vs)
        -> case findIndex (hitVertex hitRadius (mx-ox,my-oy)) vs of
             Just i -> appData {appMode = ModeDraw (APath ty (take i vs))}
             Nothing -> appData

mouseRight :: Double -> Double -> CrsState -> AppData -> AppData
mouseRight mx my crst appData@(AppData mode ox oy graph fname _)
  = case mode of
      ModeView _ _
        -> case hitTest hitRadius (mx-ox,my-oy) graph of
             (Nothing, Just i)
               -> appData {arcGraph = modifyCross i (crsRotateState) graph,
                           changed = True}
             _ {-- else --}
               -> appData
      ModeMove _ _
        -> appData
      ModeDraw p
        -> appData {
             appMode = ModeView [] [],
             arcGraph = insertArcs crst [p] graph }

mouseDrag :: Double -> Double -> AppData -> AppData
mouseDrag mx my appData@(AppData mode ox oy graph fname _)
  = case mode of
      ModeView vs is
        -> if null vs && null is
           then appData {appMode = ModeMove MoveOrigin (mx,my)}
           else appData {appMode = ModeMove (MoveVertex vs is) (mx,my)}
      ModeMove _ _
        -> appData
      ModeDraw _
        -> appData

mouseUp :: Double -> Double -> AppData -> AppData
mouseUp mx my appData@(AppData mode ox oy graph fname _)
  = case mode of
      ModeView _ _
        -> appData
      ModeMove MoveOrigin (ox,oy)
        -> appData {
             appMode = ModeView [] [],
             origX = origX appData+mx-ox,
             origY = origY appData+my-oy }
      ModeMove (MoveVertex vs is) (ox,oy)
        -> let (AGraph _ cs) = graph
               vstot = vs ++ concatMap crsVrtx (map (cs!!) is)
           in case moveVrtx vstot (mx-ox) (my-oy) graph of
                Just ag ->
                  let vs' = map (\x -> (fst x + (mx-ox), snd x + (my-oy))) vs
                  in appData {appMode = ModeView vs' is,
                              arcGraph = ag,
                              changed = True }
                Nothing ->
                  appData {appMode = ModeView vs is}
      ModeDraw _
        -> appData
