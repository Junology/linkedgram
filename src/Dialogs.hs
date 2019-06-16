------------------------------------------------
-- |
-- Module    :  Dialog
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Dialogs used in the application
--
------------------------------------------------

module Dialogs where

import qualified Data.Text as T
import Data.List
import Control.Monad

import System.Directory

import Graphics.UI.Gtk

import Config

showSaveDialog :: Maybe Window -> String -> IO (Maybe String)
showSaveDialog mwindow fname = do
  saveDlg <- fileChooserDialogNew Nothing mwindow FileChooserActionSave [("Save",ResponseOk),("Cancel",ResponseCancel)]
  if null fname
    then getCurrentDirectory >>= fileChooserSetCurrentFolder saveDlg
    else fileChooserSetFilename saveDlg fname
  resp <- dialogRun saveDlg
  maybeFn <- fileChooserGetFilename saveDlg
  widgetDestroy saveDlg
  case resp of
    ResponseOk -> return maybeFn
    _ {-else-} -> return Nothing

showOpenDialog :: Maybe Window -> String -> IO (Maybe String)
showOpenDialog mwindow fname = do
  openDlg <- fileChooserDialogNew Nothing mwindow FileChooserActionOpen [("Open",ResponseOk),("Cancel",ResponseCancel)]
  if null fname
    then getCurrentDirectory >>= fileChooserSetCurrentFolder openDlg
    else fileChooserSetFilename openDlg fname
  resp <- dialogRun openDlg
  maybeFn <- fileChooserGetFilename openDlg
  widgetDestroy openDlg
  case resp of
    ResponseOk -> return maybeFn
    _ {-else-} -> return Nothing

showAboutDialog :: Maybe Window -> IO ()
showAboutDialog mayparent = do
  aboutDlg <- aboutDialogNew
  case mayparent of
    Just win -> set (toWindow aboutDlg) [ windowTransientFor := win ]
    Nothing -> return ()
  set aboutDlg [
    aboutDialogProgramName := appName,
    aboutDialogVersion := appVersion,
    aboutDialogAuthors := [appAuthor],
    aboutDialogCopyright := appCopyright,
    aboutDialogLicense := Just ("The project is released under " ++ appLicense ++ " license. See LICENSE."),
    aboutDialogComments := "Drawing Link Projection Diagrams",
    aboutDialogWebsite := appURL,
    aboutDialogLogo := Nothing]
  dialogRun aboutDlg
  widgetDestroy aboutDlg

-------------------
-- Export Dialog --
-------------------
data OutFormat = TikzOut -- | XyPicOut
  deriving (Show,Read)

data OutDataSort = OutLink | OutSmooth
  deriving Show

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
  vbox <- dialogGetUpper expDlg
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
