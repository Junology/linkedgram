{-# LANGUAGE CPP #-}

module Config where

appName :: String
appName = PROJECT_NAME -- defined in package.yaml

appVersion :: String
appVersion = PROJECT_VERSION -- defined in package.yaml

appAuthor :: String
appAuthor = PROJECT_AUTHOR_GIVEN ++ " " ++ PROJECT_AUTHOR_FAMILY -- defined in package.yaml

appCopyright :: String
appCopyright = PROJECT_COPYRIGHT_YEAR ++ " " ++ appAuthor -- defined in package.yaml

appLicense :: String
appLicense = PROJECT_LICENSE -- defined in package.yaml

appURL :: String
appURL = PROJECT_HOMEPAGE -- defined in package.yaml
