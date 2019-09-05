------------------------------------------------
-- |
-- Module    :  Dialog.Export
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Export dialog
--
------------------------------------------------

module Dialog.Export where

import Data.Text (Text)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Data.List
import Data.Maybe
import Control.Monad

import System.IO
import System.Directory

import Graphics.UI.Gtk

import Text.TeXout
import Config
import Dialog

import ArcGraph
import ArcGraph.Common
import ArcGraph.TikZ

-- | Output format type
data OutFormat = TikzOut -- | XyPicOut
  deriving (Show,Read)

-- | What data to be exported
data OutDataSort = OutLink | OutSmooth
  deriving Show

-- | Export setting
data ExportConfig = ExportConfig {
  format :: OutFormat,
  outSort :: OutDataSort,
  filePath :: String}
  deriving Show

showExportDialog :: ArcGraph -> Maybe Window -> IO ()
showExportDialog ag mayparent = do
  -- Create Dialog and set the parent
  expDlg <- dialogNew
  case mayparent of
    Just win -> set expDlg [windowTransientFor := win]
    Nothing  -> return()

  -- Create child widgets
  -- ComboBox to select the output format
  fmtCombo <- comboBoxNewText
  comboBoxInsertText fmtCombo 0 (T.pack $ show TikzOut)
  --comboBoxInsertText fmtCombo 1 (T.pack $ show XyPicOut)
  comboBoxSetActive fmtCombo 0
  -- CheckButtons to determine whether components are colored or not.
  checkColor <- checkButtonNewWithLabel "Color components"
  -- SpinButton to input the size of the output picture.
  hlabelSize <- labelNew (Just "Output size:")
  miscSetAlignment hlabelSize 1.0 0.5
  spinSize <-spinButtonNewWithRange 0 10  0.1
  spinButtonSetDigits spinSize 2
  spinButtonSetValue spinSize 2.0
  hboxSize <- hBoxNew False 3
  boxPackStart hboxSize hlabelSize PackNatural 0
  boxPackStart hboxSize spinSize PackNatural 0

  -- Add children
  vbox <- castToBox <$> dialogGetContentArea expDlg
  containerAdd vbox fmtCombo
  containerAdd vbox checkColor
  containerAdd vbox hboxSize
  dialogAddButton expDlg "Export" ResponseOk
  dialogAddButton expDlg "Cancel" ResponseCancel
  widgetShowAll expDlg

  -- Run dialog
  resp <- dialogRun expDlg
  when (resp == ResponseOk) $ do
    maypath <- showSaveDialog mayparent ""
    when (isJust maypath) $ do
      -- Get file path name
      let filepath = fromJust maypath
      -- Get configuraions
      mayfmt <- fmap (fmap (read . T.unpack)) $ comboBoxGetActiveText fmtCombo
      colored <- toggleButtonGetActive checkColor
      normAG <- normalize <$> spinButtonGetValue spinSize <*> pure (slimCross ag)
      -- Open the file exported
      handle <- openFile filepath WriteMode
      T.hPutStrLn handle $ documentClass "" "standalone"
      T.hPutStrLn handle $ T.empty
      case mayfmt of
        Just TikzOut -> do
          T.hPutStrLn handle $ usePackage "" "tikz"
          T.hPutStrLn handle $ T.empty
          T.hPutStrLn handle $ beginEnv "document" []
          T.hPutStrLn handle $ beginEnv "tikzpicture" []
          if colored
            then T.hPutStrLn handle $ showArcGraphTikzWithComp 0.15 normAG
            else T.hPutStrLn handle $ showArcGraphTikz 0.15 normAG
          T.hPutStrLn handle $ endEnv "tikzpicture"
        otherwise -> return ()
      T.hPutStrLn handle $ endEnv "document"
      hClose handle
  widgetDestroy expDlg
