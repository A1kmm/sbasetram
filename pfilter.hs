{-# LANGUAGE RankNTypes #-}
import Text.ParserCombinators.Parsec
import Text.ParserCombinators.Parsec.Token
import Text.ParserCombinators.Parsec.Language
import Data.Char
import System.Environment
import System.Console.GetOpt
import Control.Monad

data ProgramOptions = ProgramOptions { poshowHelp :: Bool, poPosteriorCutoff :: Maybe Double, poSiteFile :: Maybe String }
emptyProgramOptions = ProgramOptions { poshowHelp = False, poPosteriorCutoff = Nothing, poSiteFile = Nothing }

options :: [OptDescr (ProgramOptions -> ProgramOptions)]
options =
    [
     Option ['h'] ["help"] (NoArg (\opts -> opts {poshowHelp = True } )) "Displays this help message",
     Option ['p'] ["posterior-cutoff"] (ReqArg (\val opts -> opts { poPosteriorCutoff = Just (read val)} ) "DOUBLE") "Minimum value required for log posterior probability in order for data to be included in output",
     Option ['s'] ["site-file"] (ReqArg (\val opts -> opts {poSiteFile = Just val}) "FILE") "File containing the sites and probabilities identfied by sbasetram"
    ]

filterOneProbeByCutoff cutoff = 
    filter ((>=cutoff) . snd)
filterByCutoff cutoff vals =
    let
        x = map (\(f, s) -> (f, filterOneProbeByCutoff cutoff s)) vals
    in
      filter (not . null . snd) x
    
valuesParser l =
    (
     do
       probe <- oneProbeParser
       valuesParser (probe:l)
    ) <|> 
    (
     do
       eof
       return l
    )
oneProbeParser =
    do
      char '>'
      name <- manyTill (noneOf "\n") (char '\n')
      l <- tfListParser []
      return (name, l)

tfListParser l =
    (
     do
       tf <- manyTill (noneOf ">\n") (try $ string " -")
       prob <- float haskell
       tfListParser $ (tf, -prob):l
    ) <|> (return l)

readValues sitefile =
    do
      result <- parseFromFile (valuesParser []) sitefile
      case result of
        Left err -> error $ (showString "Cannot parse site file: " . shows err) ""
        Right val -> return val

writeFASTA probes =
    forM probes $ \(name, vals) ->
        do
          putStrLn $ '>':name
          forM vals $ \(tfname, _) ->
              putStrLn tfname

pfilter postco sitefile =
    do
      vals <- readValues sitefile
      writeFASTA (filterByCutoff postco vals)

usageString = "Usage: pfilter [opts...]"
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
            ProgramOptions { poPosteriorCutoff = Nothing } ->
                commandLineInfo "Posterior cutoff not specified"
            ProgramOptions { poSiteFile = Nothing } ->
                commandLineInfo "Site file not specified"
            ProgramOptions { poPosteriorCutoff = Just postco, poSiteFile = Just sitefile } ->
                do
                  pfilter postco sitefile
                  return ()
    (_, _, e) -> commandLineInfo $ "Command line parse error:\n" ++ (concat e)
