import System.Environment
import System.Console.GetOpt
import TSMParser
import Control.Monad
import Data.List

buildGeneList mf =
    do
      mats <- loadTRANSFACMatrices mf
      let tfs = nub $ map fst mats
      forM_ tfs putStrLn

data ProgramOptions = ProgramOptions { poshowHelp :: Bool, poMatrixFile :: Maybe String }
emptyProgramOptions = ProgramOptions { poshowHelp = False, poMatrixFile = Nothing }

options :: [OptDescr (ProgramOptions -> ProgramOptions)]
options =
    [
     Option ['h'] ["help"] (NoArg (\opts -> opts {poshowHelp = True } )) "Displays this help message",
     Option ['m'] ["matrix file"] (ReqArg (\val opts -> opts { poMatrixFile = Just val} ) "FILENAME") "The matrix file to process"
    ]

usageString = "Usage: BuildGeneList [opts...]"
commandLineInfo x = do
  putStrLn x
  putStrLn (usageInfo usageString options)

main = do
  args <- getArgs
  case getOpt RequireOrder options args of
    (o,_,[]) ->
        let
            opts = foldl (\opt f -> f opt) emptyProgramOptions o
        in
          case opts of
            ProgramOptions { poshowHelp = True } ->
                commandLineInfo "Help information"
            ProgramOptions { poMatrixFile = Nothing } ->
                commandLineInfo "Matrix file not specified"
            ProgramOptions { poMatrixFile = Just mf } ->
                do
                  buildGeneList mf
                  return ()
    (_, _, e) -> commandLineInfo $ "Command line parse error:\n" ++ (concat e)
