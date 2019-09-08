{-# LANGUAGE Strict, StrictData #-}
{-# LANGUAGE FlexibleInstances, TypeSynonymInstances #-}

------------------------------------------------
-- |
-- Module    :  Dialog.Khovanov.F2
-- Copyright :  (c) Jun Yoshida 2019
-- License   :  BSD3
--
-- Dialog to indicate the progress of the computation of Khvanov homology with coefficients in F2.
--
------------------------------------------------
module Dialog.Khovanov.F2 where


import Control.DeepSeq (force)
import Control.Monad

import System.IO (hClose)
import Data.IORef

import Data.Maybe
import Data.Map.Strict (Map)
import qualified Data.Map.Strict as Map
import Data.IntMap.Strict (IntMap)
import qualified Data.IntMap.Strict as IMap

import qualified Data.Text as T
import qualified Data.Text.IO as T

import Numeric.Natural (Natural)

import Text.TeXout
import ArcGraph
import ArcGraph.Common
import ArcGraph.Component
import ArcGraph.State
import ArcGraph.EnhancedState
import ArcGraph.Cairo
import ArcGraph.TikZ

import Dialog.Khovanov.Common

import Numeric.F2

-----------------------
-- * Orphan instances
-----------------------
instance TeXMathShow (KHData F2 ds e) where
  texMathShow kh = (mathBB "F" <> T.pack "_2") `ofRank` rank kh

-- | Export Khovanov homology with coefficients in F2.
exportKhovanovF2 :: ExportConfig -> IO ()
exportKhovanovF2 cfg@(ExpConfig _ ag states maxQDeg nPCrs hasBndry modif _ _) = do
  let slimAG = slimCross $ normalize 1.0 ag
      numCrs = countCross ag
  -- Compute Khovanov homology
  khMapRef <- newIORef (Map.empty :: Map (Int,Int) (KHData F2 IListState (MapEState (ArcBits Natural))))
  forM_ [-maxQDeg..maxQDeg+numCrs] $ \j -> do
    -- Pass ag instead of slimAG.
    let !jkhMap = force $ computeKhovanov ag j states hasBndry
    forM_ (IMap.toList jkhMap) $ \ikh -> do
      let (i,kh) = ikh
      modifyIORef' khMapRef (Map.insert (modif i j) kh)
  khMap <- readIORef khMapRef

  -- Export TeX source
  handle <- openWithCfg cfg
  writeTeXBegining handle
  writeTeXLink cfg handle
  writeTeXCohomologyTable khMap handle
  writeTeXMaybeTOC cfg handle
  forM_ (Map.toList khMap) $ \ijkh -> do
    let ((i,j),kh) = ijkh
        sectionName = "The group $Kh^{" ++ show i ++ "," ++ show j ++ "}$"
    T.hPutStrLn handle $ macro "section*" [FixArg sectionName ]
    T.hPutStrLn handle $ macro "addcontentsline" [
      FixArg "toc", FixArg "section", FixArg sectionName ]
    T.hPutStrLn handle $ beginEnv "equation*" []
    T.hPutStrLn handle $ texMathShow kh
    T.hPutStrLn handle $ endEnv "equation*"
    T.hPutStrLn handle $ T.empty
    T.hPutStrLn handle $ macro "subsection*" [FixArg "Generating Cycle"]
    forM_ (cycleV kh) $ \cyc -> do
      let (isMult,tex) = showStateSumTikzMultlined slimAG cyc 5
          envName = if isMult
                    then "multline*"
                    else "equation*"
      T.hPutStrLn handle $ beginEnv envName []
      T.hPutStrLn handle $ tex
      T.hPutStrLn handle $ endEnv envName
    T.hPutStrLn handle $ T.empty
    when (isJust (bndryV kh)) $ do
      T.hPutStrLn handle $ macro "subsection*" [FixArg "Boundaries"]
      forM_ (fromMaybe [] (bndryV kh)) $ \bnd -> do
        let (isMult,tex) = showStateSumTikzMultlined slimAG bnd 5
            envName = if isMult
                      then "multline*"
                      else "equation*"
        T.hPutStrLn handle $ beginEnv envName []
        T.hPutStrLn handle $ tex
        T.hPutStrLn handle $ endEnv envName
      T.hPutStrLn handle $ T.empty
  writeTeXEnding handle
  hClose handle
