import qualified Data.Array.Unboxed as A
import Data.Array.IArray ((!))
import Data.Maybe
import TSMParser
import FSAParser
import Data.List
import Control.Monad
import qualified Data.ByteString as B
import Control.Exception
import System.Console.GetOpt
import System.Environment

data TSMMatrix = TSMMatrix {
      matrixFreqs :: A.UArray (Int, Int) Double,
      matrixInformationVector :: A.UArray Int Double,
      matrixMinFreq :: A.UArray Int Double,
      matrixMaxFreq :: A.UArray Int Double,
      matrixCore :: [Int],
      coreMinScore :: Double,
      coreMaxScore :: Double,
      matrixMinScore :: Double,
      matrixMaxScore :: Double
    }

informationContribByFrequency 0 = 0.0
informationContribByFrequency f = f * log(4 * f)

combineFreqs :: ([Double] -> Double) -> A.UArray (Int, Int) Double -> A.UArray Int Double
combineFreqs f m =
    let
        (_, (_, n)) = A.bounds m
    in
      A.listArray (0,n) $ map (\j -> f (map (\i -> m ! (i, j)) [0..3])) [0..n]

coreSize = 5
findCoreSequences f =
    let
        cutoff = (last . take coreSize . sort . A.elems) f
    in
      (take 5 . findIndices (>=cutoff) . A.elems) f

computeScoreOverIndices idxs information freqs =
    sum $ map (\i -> (information ! (assert (i <= ((snd . A.bounds) information)) i)) * (freqs ! i) ) idxs

summariseMatrix freqs =
    let
        infoVector = combineFreqs sum (A.amap informationContribByFrequency freqs)
        minFreqs = combineFreqs minimum freqs
        maxFreqs = combineFreqs maximum freqs
        core = findCoreSequences maxFreqs
        (_, (_, n)) = A.bounds freqs
    in
      TSMMatrix { matrixFreqs = freqs,
                  matrixInformationVector = infoVector,
                  matrixMinFreq = minFreqs,
                  matrixMaxFreq = maxFreqs,
                  matrixCore = core,
                  coreMinScore = computeScoreOverIndices core infoVector minFreqs,
                  coreMaxScore = computeScoreOverIndices core infoVector maxFreqs,
                  matrixMinScore = computeScoreOverIndices [0..n] infoVector minFreqs,
                  matrixMaxScore = computeScoreOverIndices [0..n] infoVector maxFreqs
                }
summariseMatrices = map (\(n, m) -> (n, summariseMatrix m))

searchForMatrices s bs mats =
    let
        l = B.length bs
    in
      concatMap (checkForMatchingMatricesAt s bs mats) [0..(l-1)]

checkForMatchingMatricesAt s bs mats start =
   mapMaybe (\(i,m) -> (liftM $ const (i, start)) (checkMatrixMatch s bs start m)) mats

checkMatrixMatch :: Cutoffs -> B.ByteString -> Int -> TSMMatrix -> Maybe Double
checkMatrixMatch s bs start mat =
    let
        lrem = (B.length bs) - start
        (_, lmat) = A.bounds (matrixInformationVector mat)
    in
      if lmat < lrem
      then
          unsafeCheckMatrixMatch s bs start mat lmat
      else
          Nothing

data Cutoffs = Cutoffs {
      mssCutoff :: Double,
      cssCutoff :: Double
};

normaliseScore minScore maxScore currentScore =
    (currentScore - minScore) / (maxScore - minScore)

computeScoreOver :: [Int] -> B.ByteString -> Int -> A.UArray (Int, Int) Double -> A.UArray Int Double -> Double
computeScoreOver idxs seqbs offset freqs iv =
    sum $ map (\idx -> (freqs!((fromIntegral $ B.index seqbs (idx + offset)), idx)) * (iv!idx)) idxs

unsafeCheckMatrixMatch :: Cutoffs -> B.ByteString -> Int -> TSMMatrix -> Int -> Maybe Double
unsafeCheckMatrixMatch s bs start mat n =
    let
        css = normaliseScore (coreMinScore mat) (coreMaxScore mat)
                (computeScoreOver (matrixCore mat) bs start (matrixFreqs mat) (matrixInformationVector mat))
        mss = normaliseScore (matrixMinScore mat) (matrixMaxScore mat)
              (computeScoreOver [0..n] bs start (matrixFreqs mat) (matrixInformationVector mat))
    in
      if {- trace ((showString "CSS = " . shows css)"") $ -} css > (cssCutoff s) && mss > (mssCutoff s)
      then
          Just mss
      else
          Nothing

data ProgramOptions = ProgramOptions {
      poshowHelp :: Bool,
      pocssCutoff :: Maybe Double,
      pomssCutoff :: Maybe Double,
      pofastaFile :: Maybe String,
      pofsaFile :: Maybe String
    }
defaultOptions = ProgramOptions { poshowHelp = False, pocssCutoff = Nothing, pomssCutoff = Nothing,
                                  pofastaFile = Nothing, pofsaFile = Nothing }

options :: [OptDescr (ProgramOptions -> ProgramOptions)]
options =
    [ 
      Option ['h'] ["help"] (NoArg (\opts -> opts {poshowHelp = True } )) "Displays this help message"
    , Option ['c'] ["css-cutoff"] (ReqArg (\val -> \opts -> opts { pocssCutoff = Just (read val) }) "DOUBLE") "Cut-off for core similarity score"
    , Option ['m'] ["mss-cutoff"] (ReqArg (\val -> \opts -> opts { pomssCutoff = Just (read val) }) "DOUBLE") "Cut-off for matrix similarity score"
    , Option ['M'] ["matrix-file"] (ReqArg (\val -> \opts -> opts { pofastaFile = Just val }) "FILE") "File contain TF matrices"
    , Option ['f'] ["probe-file"] (ReqArg (\val -> \opts -> opts { pofsaFile = Just val }) "FILE") "File containing all probe sequences"
    ]
usageString = "Usage: gmatim [opts...]"

commandLineInfo x = do
  putStrLn x
  putStrLn (usageInfo usageString options)

main = do
  args <- getArgs
  case getOpt RequireOrder options args of
    (o,_,[]) ->
        let
            opts = foldl (\opt -> \f -> f opt) defaultOptions o
        in
          case opts of
            ProgramOptions { poshowHelp = True } ->
                commandLineInfo "Help information"
            ProgramOptions { pocssCutoff = Nothing } ->
                commandLineInfo "CSS Cutoff not specificed"
            ProgramOptions { pomssCutoff = Nothing } ->
                commandLineInfo "MSS Cutoff not specificed"
            ProgramOptions { pofastaFile = Nothing } ->
                commandLineInfo "FASTA (matrix) file not specificed"
            ProgramOptions { pofsaFile = Nothing } ->
                commandLineInfo "FSA (sequence) file not specificed"
            ProgramOptions { pofastaFile = Just fasf, pofsaFile = Just fsaf, pomssCutoff = Just mssco, pocssCutoff = Just cssco } ->
                gmatimMain fasf fsaf (Cutoffs mssco cssco)
    (_, _, e) -> commandLineInfo $ "Command line parse error:\n" ++ (concat e)
-- "/home/andrew/Documents/TSM/matrices.fasta"
-- "/home/andrew/tgz/yeast_Young_6k.fsa"
-- defaultCutoffs = Cutoffs { mssCutoff = 0.75, cssCutoff = 0.7 }


gmatimMain fastaFile fsaFile cutoffs =
    do
      matricesAsWeights <- loadTRANSFACMatrices fastaFile
      fsaData <- loadFSAData fsaFile
      let
          matrixSummary = summariseMatrices matricesAsWeights
        in
          forM_ fsaData $ \(probe, seqbs) ->
          let
              matches = nub (map fst $ searchForMatrices cutoffs seqbs matrixSummary)
            in
              if null matches
              then
                  return []
              else
                  do
                    putStrLn $ showString ">" probe
                    forM matches $ \matrix ->
                        putStrLn matrix
