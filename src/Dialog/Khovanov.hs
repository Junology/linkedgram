{-# LANGUAGE FlexibleInstances, TypeSynonymInstances #-}

------------------------------------------------
-- |
-- Module    :  Dialog.Khovanov
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Dialog to play with Khovanov chain complex
--
------------------------------------------------

module Dialog.Khovanov where

import Control.Monad

import Data.Functor ((<&>))
import Data.Maybe
import Data.List as L
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import Data.IORef

import Graphics.UI.Gtk
import qualified Graphics.Rendering.Cairo as Cairo

import ArcGraph
import ArcGraph.Common
import ArcGraph.Component
import ArcGraph.State
import ArcGraph.EnhancedState
import ArcGraph.Cairo

import Dialog

import Dialog.Khovanov.Common
import Dialog.Khovanov.Integral
import Dialog.Khovanov.F2

{-- for debug
import System.IO
import Debug.Trace
--}

--------------------------
-- * The view for states
--------------------------
-- | Size of Pixbuf in the icon view.
pictSize :: Num a => a
pictSize = 80

-- | Generate Pixbuf of smoothing diagrams
genPixbufArcGraph :: DState ds => ArcGraph -> Int -> IO ([ds], Map.Map ds Pixbuf)
genPixbufArcGraph ag deg = do
  let states = listStates ag deg
      normAg = normalize (pictSize / 2.0) ag
  pixMapRef <- newIORef (Map.empty :: Map.Map ds Pixbuf)
  surface <- Cairo.createImageSurface Cairo.FormatRGB24 pictSize pictSize
  forM_ states $ \cs -> do
    Cairo.renderWith surface $ do
      Cairo.setSourceRGB 1.0 1.0 1.0
      Cairo.paint
      Cairo.translate (pictSize / 2.0) (pictSize / 2.0)
      Cairo.scale 1.0 (-1.0)
      Cairo.setSourceRGB 1.0 0.0 0.0
      drawArcGraph (smoothing normAg cs) Nothing
    pixbuf <- pixbufNewFromSurface surface 0 0 pictSize pictSize
    modifyIORef' pixMapRef $ Map.insert cs pixbuf
  Cairo.surfaceFinish surface
  pixMap <- readIORef pixMapRef
  return (states, pixMap)

-- | Create IconView containing smoothings of a designated degree
createArcGraphView :: (DState ds) => ArcGraph -> Int -> IO (ListStore ds,IconView)
createArcGraphView ag deg = do
  iconview <- iconViewNew
  (diagSt,pixMap) <- genPixbufArcGraph ag deg
  smoothList <- listStoreNew diagSt
  let colAGPixbuf = makeColumnIdPixbuf 1
      colLabel = makeColumnIdString 2
  treeModelSetColumn smoothList colAGPixbuf (pixMap Map.!)
  treeModelSetColumn smoothList colLabel (getLabel ag)
  set iconview [
    iconViewModel := Just smoothList,
    iconViewTextColumn := colLabel,
    iconViewPixbufColumn := colAGPixbuf,
    iconViewItemWidth := pictSize,
    iconViewColumnSpacing := 0,
    iconViewMargin := 0 ]
  return (smoothList, iconview)

-- | Show a dialog window to compute Khovanov homology
showKhovanovDialog :: ArcGraph -> Maybe Window -> IO ()
showKhovanovDialog ag mayparent = do
  -- Create Dialog and set the parent
  let slimAG = slimCross $ normalize 1.0 ag
  khovanovDlg <- dialogNew
  case mayparent of
    Just win -> set khovanovDlg [windowTransientFor := win]
    Nothing  -> return()
  set khovanovDlg [
    windowTitle := "Compute (unnormalized) Khovanov Chain",
    windowWindowPosition := WinPosCenter,
    windowDefaultWidth := 640,
    windowDefaultHeight := 480 ]
  -- Add buttons
  btnExport <- dialogAddButton khovanovDlg "Export" ResponseNone
  btnClose <- dialogAddButton khovanovDlg "Close" ResponseOk

  -- Add HBox to contain smoothing diagrams
  hboxSm <- hBoxNew True 0
  listMVec <- MV.unsafeNew (countCross slimAG+1)
  forM_ [0..countCross slimAG] $ \i -> do
    (smthList,arcGraphView) <- createArcGraphView slimAG i :: IO (ListStore IListState, IconView)
    set arcGraphView [
      iconViewSelectionMode := SelectionMultiple,
      iconViewColumns := 1 ]
    iconViewSelectAll arcGraphView
    MV.write listMVec i (smthList,arcGraphView)
    boxPackStart hboxSm arcGraphView PackNatural 3
  -- Add containers
  vbox <- castToBox <$> dialogGetContentArea khovanovDlg
  set vbox [ boxHomogeneous := False ]

  -- Scroll window containing the list of states
  scroll <- scrolledWindowNew Nothing Nothing
  scrolledWindowSetPolicy scroll PolicyAutomatic PolicyAutomatic
  scrolledWindowSetShadowType scroll ShadowIn
  boxPackStart vbox scroll PackGrow 0
  containerAdd scroll hboxSm

  -- Frame containing widgets to input cohomology data
  frameCohomology <- frameNew
  set frameCohomology [
    frameLabel := "Cohomology data", frameLabelXAlign := 0 ]
  tableCohomology <- tableNew 2 2 False
  containerAdd frameCohomology tableCohomology
  -- Spin button to input the number of positive crossings
  hlabelCrs <- labelNew (Just "Number of positive crossings: ")
  miscSetAlignment hlabelCrs 1.0 0.5
  let nCrs = countCross slimAG
  spinPCrs <-spinButtonNewWithRange 0 (fromIntegral nCrs) 1
  hlabelCoeff <- labelNew (Just "Coefficients: ")
  miscSetAlignment hlabelCoeff 1.0 0.5
  radioCoeffZ <- radioButtonNewWithLabel "Integer"
  radioCoeffF2 <- radioButtonNewWithLabelFromWidget radioCoeffZ "char 2"
  hboxCoeff <- hBoxNew False 2
  boxPackEnd hboxCoeff radioCoeffF2 PackNatural 0
  boxPackEnd hboxCoeff radioCoeffZ PackNatural 0
  tableAttachDefaults tableCohomology hlabelCrs 0 1 0 1
  tableAttach tableCohomology spinPCrs 1 2 0 1 [Fill] [Fill] 3 3
  tableAttachDefaults tableCohomology hlabelCoeff 0 1 1 2
  tableAttach tableCohomology hboxCoeff 1 2 1 2 [Fill] [Fill] 3 3

  -- Frame containing widgets to set export options.
  frameConfig <- frameNew
  set frameConfig [
    frameLabel := "Export setting", frameLabelXAlign := 0 ]
  vboxConfig <- vBoxNew False 2
  containerAdd frameConfig vboxConfig
  -- CheckButton to detemine wheter boundaries are exported or not.
  checkBndry <- checkButtonNewWithLabel "Export boundaries"
  toggleButtonSetActive checkBndry True -- Default: exported.
  -- CheckButton to determine degree convention for quantum-degree.
  checkSlim <- checkButtonNewWithLabel "Slim Table of Cohomology"
  toggleButtonSetActive checkSlim False -- Default: not slim.
  -- CheckButton to determine whether components are colored or not.
  checkColor <- checkButtonNewWithLabel "Color the components"
  toggleButtonSetActive checkColor True -- Default: colored
  -- CheckButton to determine whether TOC appears or not.
  checkTOC <- checkButtonNewWithLabel "Make Table-Of-Contents"
  toggleButtonSetActive checkTOC False -- Default: no TOC
  -- Pack CheckButton to vboxConfig
  boxPackStart vboxConfig checkBndry PackNatural 0
  boxPackStart vboxConfig checkSlim PackNatural 0
  boxPackStart vboxConfig checkColor PackNatural 0
  boxPackStart vboxConfig checkTOC PackNatural 0

  -- Add all child widgets to the main vbox
  boxPackStart vbox frameCohomology PackNatural 3
  boxPackStart vbox frameConfig PackNatural 0
  widgetShowAll khovanovDlg

  btnExport `on` buttonActivated $ do
    -- Write to the file
    mayfname <- showSaveDialog mayparent ""
    when (isJust mayfname) $ do
      -- Get selected states
      agViewListV <- V.unsafeFreeze listMVec
      states <- fmap (V.foldl' (++) []) $ V.forM agViewListV $ \agvlv -> do
        let (smthList,agView) = agvlv
        mapM (listStoreGetValue smthList)
          =<< (map head <$> iconViewGetSelectedItems agView :: IO [Int])
      -- The bound for quantum degrees
      let maxQDeg = let (AGraph ps _) = slimAG in L.length ps
      -- Read settings
      exportKh <- toggleButtonGetActive radioCoeffZ >>= \b ->
        if b then return exportKhovanovZ else return exportKhovanovF2
      nPCrs <- spinButtonGetValueAsInt spinPCrs
      hasBndry <- toggleButtonGetActive checkBndry
      modif <- toggleButtonGetActive checkSlim <&> \ifslim -> \i j ->
        if ifslim
        then (i-(nCrs-nPCrs), j-2*i+nPCrs)
        else (i-(nCrs-nPCrs), j-2*(nCrs-nPCrs)+nPCrs)
      color <- toggleButtonGetActive checkColor
      toc <- toggleButtonGetActive checkTOC
      let cfg = ExpConfig {
            filePath = fromJust mayfname,
            targetArcGraph = ag,
            statesUsed = states,
            qdegBnd = maxQDeg,
            numPCrs = nPCrs,
            hasBndry = hasBndry,
            degModifier = modif,
            hasTOC = toc,
            isColored = color }
      exportKh cfg

  -- Close the dialog when "Close" is pressed
  btnClose `on` buttonActivated $ do
    agViewListV <- V.unsafeFreeze listMVec
    -- Clear the list of states manually to avoid Gtk warnings.
    forM_ agViewListV $ \l ->
      listStoreClear (fst l)
    widgetDestroy khovanovDlg

  -- Run dialog
  dialogRun khovanovDlg
  return ()
