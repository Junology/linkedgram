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

import Control.DeepSeq (force)
import Control.Parallel
import Control.Parallel.Strategies

import Data.Functor ((<&>))
import Data.Maybe
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Data.List as L
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IMap
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as MV
import Data.IORef

import System.Directory
import System.IO

import Graphics.UI.Gtk
import qualified Graphics.Rendering.Cairo as Cairo

import Data.ChunkedBits
import Config
import Text.TeXout
import ArcGraph
import ArcGraph.Common
import ArcGraph.Component
import ArcGraph.State
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

-----------------------
-- * Orphan instances
-----------------------
instance TeXMathShow (KHData ds e) where
  texMathShow kh = finAbGroup (rank kh) (tors kh)

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
  -- Spin button to input the number of positive crossings
  hlabelCrs <- labelNew (Just "Number of positive crossings")
  miscSetAlignment hlabelCrs 1.0 0.5
  let nCrs = countCross slimAG
  spinPCrs <-spinButtonNewWithRange 0 (fromIntegral nCrs) 1
  -- CheckButton to detemine wheter boundaries are exported or not.
  checkBndry <- checkButtonNewWithLabel "Export boundaries"
  -- CheckButton to determine degree convention for quantum-degree.
  checkSlim <- checkButtonNewWithLabel "Slim Table of Cohomology"
  hboxConf <- hBoxNew False 3
  boxPackStart hboxConf hlabelCrs PackNatural 0
  boxPackStart hboxConf spinPCrs PackNatural 0
  boxPackStart hboxConf checkBndry PackNatural 0
  boxPackStart hboxConf checkSlim PackNatural 0
  boxPackStart vbox hboxConf PackNatural 0
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
      -- Read options
      hasBndry <- toggleButtonGetActive checkBndry
      nPCrs <- spinButtonGetValueAsInt spinPCrs
      modif <- toggleButtonGetActive checkSlim <&> \ifslim -> \i j ->
        if ifslim
        then (i-(nCrs-nPCrs), j-2*i+nPCrs)
        else (i-(nCrs-nPCrs), j-2*(nCrs-nPCrs)+nPCrs)

      -- Compute Khovanov homology
      khMapRef <- newIORef (Map.empty :: Map (Int,Int) (KHData IListState (MapEState (ArcBits ChunkedBits))))
      forM_ [-maxQDeg..maxQDeg] $ \j -> do
        -- Pass ag instead of slimAG.
        let !jkhMap = force $ computeKhovanov ag j states hasBndry
        forM_ (IMap.toList jkhMap) $ \ikh -> do
          let (i,kh) = ikh
          modifyIORef' khMapRef (Map.insert (modif i j) kh)
      khMap <- readIORef khMapRef

      -- Export the TeX source
      handle <- openFile (fromJust mayfname) WriteMode
      T.hPutStrLn handle $ documentClass "pdftex,a4paper" "scrartcl"
      T.hPutStrLn handle $ T.empty
      T.hPutStrLn handle $ usePackage "" "amsmath,amssymb"
      T.hPutStrLn handle $ usePackage "" "tikz"
      T.hPutStrLn handle $ usePackage "" "breqn"
      T.hPutStrLn handle $ macro "allowdisplaybreaks" [OptArg "2"]
      T.hPutStrLn handle $ T.empty
      T.hPutStrLn handle $ beginEnv "document" []
      T.hPutStrLn handle $ beginEnv "center" []
      T.hPutStrLn handle $ beginEnv "tikzpicture" []
      T.hPutStrLn handle $ showArcGraphTikzWithComp 0.15 (normalize 2.0 slimAG)
      T.hPutStrLn handle $ endEnv "tikzpicture"
      T.hPutStrLn handle $ endEnv "center"
      T.hPutStrLn handle $ T.empty
      T.hPutStrLn handle $ macro "section*" [FixArg "Table of Homology Groups"]
      T.hPutStrLn handle $ beginEnv "equation*" []
      T.hPutStrLn handle $ texMathShow khMap
      T.hPutStrLn handle $ endEnv "equation*"
      T.hPutStrLn handle $ T.empty
      T.hPutStrLn handle $ macro "tableofcontents" []
      T.hPutStrLn handle $ T.empty
      forM_ (Map.toList khMap) $ \ijkh -> do
        let ((i,j),kh) = ijkh
        T.hPutStrLn handle $ macro "section*" [FixArg $ "The group $Kh^{" ++ show i ++ "," ++ show j ++ "}$"]
        T.hPutStrLn handle $ beginEnv "equation*" []
        T.hPutStrLn handle $ texMathShow kh
        T.hPutStrLn handle $ endEnv "equation*"
        T.hPutStrLn handle $ T.empty
        T.hPutStrLn handle $ macro "subsection*" [FixArg "Generating Cycle"]
        forM_ (cycleV kh) $ \cyc -> do
          T.hPutStrLn handle $ beginEnv "dmath*" []
          T.hPutStrLn handle $ showStateSumTikzText slimAG cyc
          T.hPutStrLn handle $ endEnv "dmath*"
        T.hPutStrLn handle $ T.empty
        when (isJust (bndryV kh)) $ do
          T.hPutStrLn handle $ macro "subsection*" [FixArg "Boundaries"]
          forM_ (fromMaybe [] (bndryV kh)) $ \bnd -> do
            T.hPutStrLn handle $ beginEnv "dmath*" []
            T.hPutStrLn handle $ showStateSumTikzText slimAG bnd
            T.hPutStrLn handle $ endEnv "dmath*"
          T.hPutStrLn handle $ T.empty
      T.hPutStrLn handle $ endEnv "document"
      hClose handle
      -- writeFile (fromJust mayfname) (docKhovanovTikz ag "pdftex,a4paper" "scrartcl" khMap)

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
