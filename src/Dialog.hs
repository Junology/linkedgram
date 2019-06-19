------------------------------------------------
-- |
-- Module    :  Dialog
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Basic dialogs
--
------------------------------------------------

module Dialog (
  showSaveDialog,
  showOpenDialog,
  showAboutDialog
  ) where

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
