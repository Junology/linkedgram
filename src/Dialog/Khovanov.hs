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
import Control.Monad.ST

import Control.Parallel
import Control.Parallel.Strategies

import Data.Functor ((<&>))
import Data.Maybe
import qualified Data.Text as T
import Data.List as L
import qualified Data.Map.Strict as Map
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import Data.IORef
import Data.STRef

import System.Directory

import Graphics.UI.Gtk
import qualified Graphics.Rendering.Cairo as Cairo

import Config
import ArcGraph
import ArcGraph.Common
import ArcGraph.EnhancedState
import ArcGraph.Cairo
import ArcGraph.TikZ

import qualified Numeric.LinearAlgebra as LA

import Numeric.Algebra.FreeModule
import Numeric.Algebra.Frobenius
import Numeric.Algebra.IntMatrix

import Dialog

{-- for debug
import System.IO
import Debug.Trace
--}

pictSize :: Num a => a
pictSize = 80

-- | Generate Pixbuf of smoothing diagrams
genPixbufArcGraph :: ArcGraph -> Int -> IO ([DiagramState], Map.Map DiagramState Pixbuf)
genPixbufArcGraph ag deg = do
  let states = listStates ag deg
      normAg = normalize (pictSize / 2.0) ag
  pixMapRef <- newIORef Map.empty
  surface <- Cairo.createImageSurface Cairo.FormatRGB24 pictSize pictSize
  forM_ states $ \cs -> do
    Cairo.renderWith surface $ do
      Cairo.setSourceRGB 1.0 1.0 1.0
      Cairo.paint
      Cairo.translate (pictSize / 2.0) (pictSize / 2.0)
      Cairo.scale 1.0 (-1.0)
      Cairo.setSourceRGB 1.0 0.0 0.0
      drawArcGraph (smoothing cs normAg) Nothing
    pixbuf <- pixbufNewFromSurface surface 0 0 pictSize pictSize
    modifyIORef' pixMapRef $ Map.insert cs pixbuf
  Cairo.surfaceFinish surface
  pixMap <- readIORef pixMapRef
  return (states, pixMap)

-- | Create IconView containing smoothings of a designated degree
createArcGraphView :: ArcGraph -> Int -> IO (ListStore DiagramState,IconView)
createArcGraphView ag deg = do
  iconview <- iconViewNew
  (diagSt,pixMap) <- genPixbufArcGraph ag deg
  smoothList <- listStoreNew diagSt
  let colAGPixbuf = makeColumnIdPixbuf 1
      colLabel = makeColumnIdString 2
  treeModelSetColumn smoothList colAGPixbuf (pixMap Map.!)
  treeModelSetColumn smoothList colLabel $ \st ->
    let (AGraph _ cs') = ag
    in map (\c -> if fst c `elem` st then '1' else '0') $ zip [0..] cs'
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
  btnClose <- dialogAddButton khovanovDlg "Close" ResponseOk

  -- Add HBox to contain smoothing diagrams
  hboxSm <- hBoxNew True 0
  listMVec <- MV.unsafeNew (countCross slimAG+1)
  forM_ [0..countCross slimAG] $ \i -> do
    (smthList,arcGraphView) <- createArcGraphView slimAG i
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
  -- Spin button to determine a quantum degree to compute
  btnCompute <- buttonNewWithLabel "Compute"
  btnSummary <- buttonNewWithLabel "Export Summary"
  let maxQDeg = fromIntegral $ let (AGraph ps _) = slimAG in L.length ps
  spinQDeg <-spinButtonNewWithRange (-maxQDeg) maxQDeg 1
  hbtnbox <- hButtonBoxNew
  boxPackStart hbtnbox spinQDeg PackNatural 0
  boxPackStart hbtnbox btnCompute PackNatural 0
  boxPackStart hbtnbox btnSummary PackNatural 0
  boxPackStart vbox hbtnbox PackNatural 0
  widgetShowAll khovanovDlg

  -- Compute (unnormalized) Khovanov homology
  btnCompute `on` buttonActivated $ do
    -- Get q-degree
    qdeg <- spinButtonGetValueAsInt spinQDeg
    -- Get selected states
    agViewListV <- V.unsafeFreeze listMVec
    states <- fmap (V.foldl' (++) []) $ V.forM agViewListV $ \agvlv -> do
      let (smthList,agView) = agvlv
      mapM (listStoreGetValue smthList)
        =<< (map head <$> iconViewGetSelectedItems agView :: IO [Int])

    -- Execute the computation
    let khResult = computeKhovanov slimAG 0 (MV.length listMVec) qdeg states
    -- Print the result
    forM_ (Map.toList khResult) $ \khi -> do
      let ((i,_),KHData freeRk torsion cycleL bndryL) = khi
      let prettyFree = "Z^{" ++ show freeRk ++ "}"
          prettyTor = L.intercalate " (+) " $ map (\x -> "Z/" ++ show x) torsion
      let cohPretty = case (freeRk,torsion) of
                        (0,[]) -> "0"
                        (0,rs) -> prettyTor
                        (n,[]) -> prettyFree
                        (n,rs) -> prettyTor ++ " (+) " ++ prettyFree
      unless (cohPretty == "0") $
        putStrLn $ "Kh^{" ++ show i ++ "," ++ show qdeg ++ "} = " ++ cohPretty

  btnSummary `on` buttonActivated $ do
    -- Write to the file
    mayfname <- showSaveDialog mayparent ""
    when (isJust mayfname) $ do
      -- Get selected states
      agViewListV <- V.unsafeFreeze listMVec
      states <- fmap (V.foldl' (++) []) $ V.forM agViewListV $ \agvlv -> do
        let (smthList,agView) = agvlv
        mapM (listStoreGetValue smthList)
          =<< (map head <$> iconViewGetSelectedItems agView :: IO [Int])
      -- Execute computation on all quantum-degrees
      --{-- parallel version
      let !khMap = V.foldl' Map.union Map.empty $ {- withStrategy (parTraversable  rdeepseq) $ -} V.fromList [-(round maxQDeg)..(round maxQDeg)] <&> \j -> computeKhovanov slimAG 0 (MV.length listMVec) j states
      {-- non-parallel version
      khMapRef <- newIORef Map.empty
      forM_ [-(round maxQDeg)..(round maxQDeg)] $ \j ->
        let khMapj = computeKhovanov slimAG [0..(MV.length listMVec)] j states
        in modifyIORef' khMapRef $ Map.union khMapj
      khMap <- readIORef khMapRef
      --}
      writeFile (fromJust mayfname) (docKhovanovTikz ag "pdftex,a4paper" "scrartcl" khMap)

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
