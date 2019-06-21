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

import qualified Data.Text as T
import Data.List as L
import qualified Data.Map.Strict as Map
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import Data.IORef

import System.Directory

import Graphics.UI.Gtk
import qualified Graphics.Rendering.Cairo as Cairo

import Config
import ArcGraph
import ArcGraph.Common
import ArcGraph.EnhancedState
import ArcGraph.Cairo

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius
import Numeric.Algebra.IntMatrix

pictSize :: Num a => a
pictSize = fromIntegral 80

-- | Generate Pixbuf of smoothing diagrams
genPixbufArcGraph :: ArcGraph -> Int -> IO ([ArcGraph], Map.Map [Cross] Pixbuf)
genPixbufArcGraph ag deg = do
  let normAgs = listSmoothing [deg] $ normalize (pictSize / 2.0) ag
  pixMapRef <- newIORef Map.empty
  surface <- Cairo.createImageSurface Cairo.FormatRGB24 pictSize pictSize
  forM_ normAgs $ \ag -> do
    let (AGraph _ cs) = ag
    Cairo.renderWith surface $ do
      Cairo.setSourceRGB 1.0 1.0 1.0
      Cairo.paint
      Cairo.translate (pictSize / 2.0) (pictSize / 2.0)
      Cairo.scale 1.0 (-1.0)
      Cairo.setSourceRGB 1.0 0.0 0.0
      drawArcGraph ag Nothing
    pixbuf <- pixbufNewFromSurface surface 0 0 pictSize pictSize
    modifyIORef' pixMapRef $ Map.insert cs pixbuf
  Cairo.surfaceFinish surface
  pixMap <- readIORef pixMapRef
  return (normAgs, pixMap)

-- | Create IconView containing smoothings of a designated degree
createArcGraphView :: ArcGraph -> Int -> IO (ListStore ArcGraph,IconView)
createArcGraphView ag deg = do
  iconview <- iconViewNew
  (normAgs,pixMap) <- genPixbufArcGraph ag deg
  smoothList <- listStoreNew normAgs
  let colAGPixbuf = makeColumnIdPixbuf 1
      colLabel = makeColumnIdString 2
  treeModelSetColumn smoothList colAGPixbuf $ \ag ->
    let (AGraph _ cs) = ag in pixMap Map.! cs
  treeModelSetColumn smoothList colLabel $ \ag ->
    let (AGraph _ cs) = ag
        crsChar (Crs _ _ Smooth0) = '0'
        crsChar (Crs _ _ Crossing) = '*'
        crsChar (Crs _ _ Smooth1) = '1'
    in map crsChar cs
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
  btnCompute <- dialogAddButton khovanovDlg "Compute" ResponseOk
  btnClose <- dialogAddButton khovanovDlg "Close" ResponseCancel

  -- Add HBox to contain smoothing diagrams
  hboxSm <- hBoxNew True 0
  listMVec <- MV.unsafeNew (countCross ag+1)
  forM_ [0..countCross ag] $ \i -> do
    (smthList,arcGraphView) <- createArcGraphView ag i
    set arcGraphView [
      iconViewSelectionMode := SelectionMultiple,
      iconViewColumns := 1 ]
    MV.write listMVec i (smthList,arcGraphView)
    boxPackStart hboxSm arcGraphView PackNatural 3

  -- Add containers
  vbox <- castToBox <$> dialogGetContentArea khovanovDlg
  scroll <- scrolledWindowNew Nothing Nothing
  scrolledWindowSetPolicy scroll PolicyAutomatic PolicyAutomatic
  scrolledWindowSetShadowType scroll ShadowIn
  boxPackStart vbox scroll PackGrow 0
  containerAdd scroll hboxSm
  widgetShowAll khovanovDlg

  -- Compute (unnormalized) Khovanov homology
  btnCompute `on` buttonActivated $ do
    forM_ [1..(MV.length listMVec - 1)] $ \i -> do
      -- Get selecteds in the previous column
      (smthList0,agView0) <- MV.read listMVec (i-1)
      smth0 <- mapM (listStoreGetValue smthList0) =<< map head <$> iconViewGetSelectedItems agView0
      let smth0Enh = L.concatMap (enhancedStatesL . slimCross) smth0
      -- Get selecteds in the current column
      (smthList1,agView1) <- MV.read listMVec i
      smth1 <- mapM (listStoreGetValue smthList1) =<< map head <$> iconViewGetSelectedItems agView1
      let smth1Enh = L.concatMap (enhancedStatesL . slimCross) smth1
      unless (L.null smth0 || L.null smth1) $ do
        -- Print the diff matrix
        putStrLn $ "diff " ++ show i
        {-
        let (_,h,_) = smithNF $ matiDataToLA $ genMatrix differential smth0Enh smth1Enh
        print $ matiLAToData h
        -}
        let (dvec,_,_) = kerImOf $ matiDataToLA $ genMatrix differential smth0Enh smth1Enh
        print $ vectiLAToData dvec

  -- Close the dialog when "Close" is pressed
  btnClose `on` buttonActivated $ do
    widgetDestroy khovanovDlg

  -- Run dialog
  dialogRun khovanovDlg
  return ()
