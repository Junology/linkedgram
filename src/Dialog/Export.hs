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

import qualified Data.Text as T
import Data.List
import Control.Monad

import System.Directory

import Graphics.UI.Gtk

import Config
import Dialog

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

showExportDialog :: Maybe Window -> IO (Maybe ExportConfig)
showExportDialog mayparent = do
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
  -- Radio buttons to choose what to export
  radLinkOnly <- radioButtonNewWithLabel "Link Diagram only"
  radSmooth <- radioButtonNewWithLabelFromWidget radLinkOnly "All Smoothings"

  -- Add children
  vbox <- castToBox <$> dialogGetContentArea expDlg
  containerAdd vbox fmtCombo
  containerAdd vbox radLinkOnly
  containerAdd vbox radSmooth
  dialogAddButton expDlg "Export" ResponseOk
  dialogAddButton expDlg "Cancel" ResponseCancel
  widgetShowAll expDlg

  -- Run dialog
  resp <- dialogRun expDlg
  mayfmt <- fmap (fmap (read . T.unpack)) $ comboBoxGetActiveText fmtCombo
  isLink <- toggleButtonGetActive radLinkOnly
  -- isSmooth <- toggleButtonGetActive radSmooth
  let sort = if isLink
             then OutLink
             else OutSmooth
  widgetDestroy expDlg
  case (resp, mayfmt) of
    (ResponseOk, Just fmt) -> do
      -- Run SaveDialog to determine an exported file
      maypath <- showSaveDialog mayparent ""
      case maypath of
        Just path
          -> return $ Just $ ExportConfig fmt sort path
        _ {-else-}
          -> return Nothing
    _ {- else -}
      -> return Nothing
