{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
{-# LANGUAGE FlexibleInstances, TypeSynonymInstances #-}

------------------------------------------------
-- |
-- Module    :  Dialog.Khovanov.Common
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Defining several data types used in exporting Khovanov homology.
--
------------------------------------------------

module Dialog.Khovanov.Common where

import GHC.Generics (Generic)
import Control.DeepSeq (NFData)

import System.IO (Handle,openFile,IOMode(..))

import Control.Monad (when)

import Data.Text (Text)
import qualified Data.Text as T (empty,singleton)
import qualified Data.Text.IO as T (hPutStrLn)

import Text.TeXout

import ArcGraph
import ArcGraph.Common
import ArcGraph.Component
import ArcGraph.State (IListState)
import ArcGraph.EnhancedState (KHData)
import ArcGraph.TikZ

-- | Export Configuration
data ExportConfig = ExpConfig {
  filePath :: String,
  targetArcGraph :: ArcGraph,
  statesUsed :: [IListState],
  qdegBnd :: Int,
  numPCrs :: Int,
  hasBndry :: Bool,
  degModifier :: Int -> Int -> (Int,Int),
  hasTOC :: Bool,
  isColored :: Bool
  } deriving (Generic,NFData)


-------------------------------------------
-- * Common part in exporting TeX sources
-------------------------------------------
openWithCfg :: ExportConfig -> IO Handle
openWithCfg cfg = openFile (filePath cfg) WriteMode

writeTeXBegining :: Handle -> IO ()
writeTeXBegining handle = do
  T.hPutStrLn handle $ documentClass "pdftex,a4paper" "scrartcl"
  T.hPutStrLn handle $ T.empty
  T.hPutStrLn handle $ usePackage "" "amsmath,amssymb"
  T.hPutStrLn handle $ usePackage "" "tikz"
  T.hPutStrLn handle $ macro "allowdisplaybreaks" [OptArg "2"]
  T.hPutStrLn handle $ T.empty
  T.hPutStrLn handle $ beginEnv "document" []

writeTeXLink :: ExportConfig -> Handle -> IO ()
writeTeXLink cfg handle = do
  let normAG = normalize 2.0 (slimCross (targetArcGraph cfg))
  T.hPutStrLn handle $ beginEnv "center" []
  T.hPutStrLn handle $ beginEnv "tikzpicture" []
  if isColored cfg
    then T.hPutStrLn handle $ showArcGraphTikzWithComp 0.15 normAG
    else T.hPutStrLn handle $ showArcGraphTikz 0.15 normAG
  T.hPutStrLn handle $ endEnv "tikzpicture"
  T.hPutStrLn handle $ endEnv "center"
  T.hPutStrLn handle $ T.empty

writeTeXCohomologyTable :: (TeXMathShow a) => a -> Handle -> IO ()
writeTeXCohomologyTable kh handle = do
  T.hPutStrLn handle $ macro "section*" [FixArg "Table of Homology Groups"]
  T.hPutStrLn handle $ beginEnv "equation*" []
  T.hPutStrLn handle $ texMathShow kh
  T.hPutStrLn handle $ endEnv "equation*"
  T.hPutStrLn handle $ T.empty

writeTeXMaybeTOC :: ExportConfig -> Handle -> IO ()
writeTeXMaybeTOC cfg handle = do
  when (hasTOC cfg) $ do
    T.hPutStrLn handle $ macro "tableofcontents" []
    T.hPutStrLn handle $ T.empty

writeTeXEnding :: Handle -> IO ()
writeTeXEnding handle = do
  T.hPutStrLn handle $ endEnv "document"
