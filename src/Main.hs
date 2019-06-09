module Main where

import Data.Text as T
import Data.List as L
import Data.IORef
import Control.Monad

import Graphics.UI.Gtk as Gtk hiding (changed)
import Graphics.Rendering.Cairo as Cairo

import Debug.Trace

import ArcGraph
import ArcGraphCairo
import ArcGraphTikz
import AppData
import Dialogs
import Config

-- For Debug
testPath1 = APath ClosedPath [(0.0,0.0),(-100.0,100.0),(-200.0,0.0),(-100.0,-100.0)]
testPath2 = APath ClosedPath [(-100.0,0.0),(0.0,100.0),(100.0,0.0),(0.0,-100.0)]

testGraph = mkArcGraph Crossing [testPath1, testPath2]

makeTitle :: AppData -> String
makeTitle appData =
  let fname = fileName appData
      changestar = if changed appData then "*" else ""
  in appName ++ " - " ++ fname ++ changestar

updateTitle :: Gtk.Window -> AppData -> IO ()
updateTitle win appData = do
  set win [windowTitle := makeTitle appData]

mainDraw :: Double -> Double -> AppData -> Cairo.Render ()
mainDraw wid hei appData = do
  translate 0.0 hei
  scale 1.0 (-1.0)

  translate (origX appData) (origY appData)
  setSourceRGB 1 1 1 -- Clear
  paint

  -- A cross on the origin
  newPath
  moveTo 0 (-10)
  lineTo 0 10
  moveTo (-10) 0
  lineTo 10 0
  setSourceRGB 0.8 0.8 0.8
  stroke

  setSourceRGB 1.0 0.0 0.0
  drawArcGraph (arcGraph appData) (Just (0.0, 0.0, 1.0))

  case appMode appData of
    ModeView vs _
      -> forM_ vs (drawVertex 1.0 0.0 0.5 5.0)
    ModeMove (MoveVertex vs _) _
      -> forM_ vs (drawVertex 1.0 0.0 0.5 5.0)
    ModeMove MoveOrigin _
      -> return ()
    ModeDraw ap@(APath _ vs) -> do
      setSourceRGB 0.0 1.0 0.0
      drawArcPath (APath OpenPath vs)
      drawArcVertex 0.0 0.0 1.0 ap

main :: IO ()
main = do
  -- Fix the size of the main window
  let winWidth = 640
  let winHeight = 480

  -- Initialize GTK+
  initGUI

  -- Create the main window
  window <- windowNew
  vbox <- vBoxNew True 10
  set window [windowTitle := appName,
              windowDefaultWidth := winWidth,
              windowDefaultHeight := winHeight,
              containerBorderWidth := 1,
              containerChild := vbox]

  -- Packing in VBox
  set vbox [boxHomogeneous := False]
  tlbar <- toolbarNew
  set tlbar [
    widgetCanFocus := False ]
  canvas <- drawingAreaNew
  boxPackStart vbox tlbar PackNatural 0
  boxPackStart vbox canvas PackGrow 0

  -- Make Toolbar
  tliNew <- toolButtonNewFromStock stockNew
  set tliNew [
    toolButtonLabel := Just "New",
    toolButtonIconName := "label"]
  tliOpen <- toolButtonNewFromStock stockOpen
  set tliOpen [
    toolButtonLabel := Just "Open",
    toolButtonIconName := "label"]
  tliSave <- toolButtonNewFromStock stockSave
  set tliSave [
    toolButtonLabel := Just "Save",
    toolButtonIconName := "label"]
  tliSaveAs <- toolButtonNewFromStock stockSaveAs
  set tliSaveAs [
    toolButtonLabel := Just "Save As",
    toolButtonIconName := "label"]
  tliExport <- toolButtonNewFromStock stockConvert
  set tliExport [
    toolButtonLabel := Just "Export",
    toolButtonIconName := "label"]
  -------------------------------------------
  tliSep1 <- separatorToolItemNew
  -------------------------------------------
  tliDraw <- toggleToolButtonNewFromStock stockEdit
  set tliDraw [
    toolButtonLabel := Just "Draw Mode",
    toolButtonIconName := "label",
    toggleToolButtonActive := False]
  tliDrawCon <- radioToolButtonNewFromStock stockConnect
  set tliDrawCon [
    toolButtonLabel := Just "Closed Path",
    toolButtonIconName := "label"]
  tliDrawDis <- radioToolButtonNewWithStockFromWidget tliDrawCon stockDisconnect
  set tliDrawDis [
    toolButtonLabel := Just "Open Path",
    toolButtonIconName := "label"]
  ---------------------------------------------
  tliSep2 <- separatorToolItemNew
  ---------------------------------------------
  tliAbout <- toolButtonNewFromStock stockAbout
  set tliAbout [
    toolButtonLabel := Just "About",
    toolButtonIconName := "label"]

  toolbarInsert tlbar tliNew (-1)
  toolbarInsert tlbar tliOpen (-1)
  toolbarInsert tlbar tliSave (-1)
  toolbarInsert tlbar tliSaveAs (-1)
  toolbarInsert tlbar tliExport (-1)
  toolbarInsert tlbar tliSep1 (-1)
  toolbarInsert tlbar tliDraw (-1)
  toolbarInsert tlbar tliDrawCon (-1)
  toolbarInsert tlbar tliDrawDis (-1)
  toolbarInsert tlbar tliSep2 (-1)
  toolbarInsert tlbar tliAbout (-1)

  -- Create Application Data
  appRef <- newIORef $ initialAppData
  readIORef appRef >>= updateTitle window
  btnRef <- newIORef False

  -- Signal Handlers for the main Window
  widgetAddEvents canvas []
  window `on` deleteEvent $ tryEvent $ liftIO mainQuit

  -- Signal Handlers for tool buttons
  -- Toolbar Buttons Events
  onToolButtonClicked tliNew $ do
    (wid,hei) <- widgetGetSize canvas
    let (ox,oy) = (fromIntegral wid/2.0, fromIntegral hei/2.0)
    writeIORef appRef $ initialAppData
    modifyIORef' appRef $ moveOrigin ox oy
    readIORef appRef >>= updateTitle window
    widgetQueueDraw canvas
  onToolButtonClicked tliOpen $ do
    maypath <- fmap fileName (readIORef appRef) >>= showOpenDialog (Just window)
    case maypath of
      Just path -> do
        content <- readFile path
        (wid, hei) <- liftIO $ widgetGetSize canvas
        let ox = fromIntegral wid / 2.0
            oy = fromIntegral hei / 2.0
        writeIORef appRef $ openAppData path ox oy (read content)
      Nothing ->
        return ()
    readIORef appRef >>= updateTitle window
    widgetQueueDraw canvas
  onToolButtonClicked tliSave $ do
    appData <- readIORef appRef
    when (changed appData) $ do
      path <- if L.null (fileName appData)
              then do
                mayfname <- showSaveDialog (Just window) ""
                case mayfname of
                  Just fn -> return fn
                  Nothing -> return ""
              else return (fileName appData)
      when (not (L.null path)) $ do
        writeFile path (show $ arcGraph appData)
        modifyIORef' appRef $ setFileName path
        modifyIORef' appRef $ confirmChange
    readIORef appRef >>= updateTitle window
  onToolButtonClicked tliSaveAs $ do
    appData <- readIORef appRef
    mayfname <- showSaveDialog (Just window) (fileName appData)
    case mayfname of
      Just fn -> do
        writeFile fn (show $ arcGraph appData)
        modifyIORef' appRef $ setFileName fn
        modifyIORef' appRef $ confirmChange
        readIORef appRef >>= updateTitle window
      Nothing -> return ()
  onToolButtonClicked tliExport $ do
    graph <- fmap arcGraph $ readIORef appRef
    mayconf <- showExportDialog (Just window)
    case mayconf of
      Just (ExportConfig fmt sort path)
        -> case (fmt,sort) of
             (TikzOut, OutLink)
               -> writeFile path (typesetArcGraphTikz "" "standalone" [graph])
             (TikzOut, OutSmooth)
               -> writeFile path (docArcGraphTikz "" "scrartcl" [graph])
      Nothing
        -> return ()
  onToolButtonToggled tliDraw $ do
    isAct <- toggleToolButtonGetActive tliDraw
    case isAct of
      True -> do
        iscon <- toggleToolButtonGetActive tliDrawCon
        let ty = if iscon then ClosedPath else OpenPath
        modifyIORef' appRef $ changeMode (ModeDraw (APath ty []))
        readIORef appRef >>= updateTitle window
      False -> do
        iscon <- toggleToolButtonGetActive tliDrawCon
        let ty = if iscon then ClosedPath else OpenPath
        modifyIORef' appRef $ \ap ->
          case appMode ap of
            ModeDraw (APath _ vs) ->
              ap {appMode = ModeView [] [],
                  arcGraph = insertArcs Crossing [APath ty vs] (arcGraph ap)}
            _ -> undefined
    readIORef appRef >>= updateTitle window
    widgetQueueDraw canvas
  onToolButtonClicked tliAbout $ do
    showAboutDialog (Just window)

  -- Signal Handlers for DrawingWindow
  widgetAddEvents canvas [PointerMotionHintMask,Button1MotionMask,ButtonPressMask,ButtonReleaseMask]
  widgetDelEvents canvas [PointerMotionMask]
  canvas `on` buttonPressEvent $ tryEvent $ do
    button <- eventButton
    click <- eventClick
    (mx,my') <- eventCoordinates
    (_, hei) <- liftIO $ widgetGetSize canvas
    let my = fromIntegral hei - my'
    modifs <- eventModifier
    case (button,click) of
      (Gtk.LeftButton, Gtk.SingleClick) -> do
        liftIO $ writeIORef btnRef True
        liftIO $ modifyIORef' appRef $ mouseDown mx my (not $ L.null modifs)
      (Gtk.LeftButton, Gtk.DoubleClick) -> do
        liftIO $ modifyIORef' appRef $ mouseDouble mx my
      (Gtk.RightButton, Gtk.SingleClick) -> do
        liftIO $ toggleToolButtonSetActive tliDraw False
        liftIO $ modifyIORef' appRef $ mouseRight mx my Crossing
    liftIO $ readIORef appRef >>= updateTitle window
    liftIO $ widgetQueueDraw canvas
  canvas `on` buttonReleaseEvent $ tryEvent $ do
    Gtk.LeftButton <- eventButton
    liftIO $ writeIORef btnRef False
    (mx,my') <- eventCoordinates
    (_, hei) <- liftIO $ widgetGetSize canvas
    let my = fromIntegral hei - my'
    liftIO $ modifyIORef' appRef $ mouseUp mx my
    liftIO $ readIORef appRef >>= updateTitle window
    liftIO $ widgetQueueDraw canvas
  canvas `on` motionNotifyEvent $ tryEvent $ do
    True <- liftIO $ readIORef btnRef
    (mx,my') <- eventCoordinates
    (_, hei) <- liftIO $ widgetGetSize canvas
    let my = fromIntegral hei - my'
    liftIO $ modifyIORef' appRef $ mouseDrag mx my

  -- Signal Handlers for DrawingWindow
  widgetAddEvents canvas [PointerMotionMask,KeyPressMask]
  canvas `on` exposeEvent $ tryEvent $ do
    drawing <- eventWindow
    appData <- liftIO $ readIORef appRef
    (wid,hei) <- liftIO $ widgetGetSize canvas
    let (w,h) = (fromIntegral wid, fromIntegral hei)
    liftIO $ renderWithDrawable drawing $ mainDraw w h appData

  -- Show Main Window
  widgetShowAll window
  (wid,hei) <- widgetGetSize canvas
  let (ox,oy) = (fromIntegral wid/2.0, fromIntegral hei/2.0)
  modifyIORef' appRef $ moveOrigin ox oy

  -- Start Main Loop
  mainGUI
