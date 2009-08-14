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
import Numeric.GSL.Special.Gamma
import qualified Data.Array.ST as S
import qualified Control.Monad.ST as S
import qualified Data.STRef as S
import Data.Word

data TSMMatrix = TSMMatrix {
      matrixFreqs :: A.UArray (Int, Int) Double,
      matrixProbs :: A.UArray (Int, Int) Double
    }

prependProbsOneEntry freqs n i l j =
    let
        f = freqs!(j, i)
        p = (beta (f + 2.0) (n - f + 1.0)) / (beta (f + 1.0) (n - f + 1.0))
    in
      ((j, i), p):l

prependProbsOnePosition freqs l i =
    let
        n = sum (\j -> freqs!(j, i)) [0..3]
    in
      foldl' (prependProbsOneEntry freqs n i) l [0..3]

computeProbs freqs = 
    let
        (_, n) = A.bounds freqs
    in
      A.array $ foldl' (prependProbsOnePosition freqs) [] [0..(n-1)]

summariseMatrix freqs = TSMMatrix {
                          matrixFreqs = freqs,
                          matrixProbs = computeProbs freqs
                        }

summariseMatrices = map (\(n, m) -> (n, summariseMatrix m))

byteStringForM :: Monad m => B.ByteString -> (Word8 -> m ()) -> m ()
byteStringForM s f = B.foldl (\m -> \w -> m >> f w) (return ()) s
byteStringPairsForM :: Monad m => B.ByteString -> (Word8 -> Word8 -> m ()) -> m ()
byteStringPairsForM s f =
    case B.uncons s of
      Nothing -> return ()
      Just (w1, s1) ->
          snd $ B.foldl (\(wp, m) -> \w -> m >> (w, f wp w)) (w1, return ()) s

modifyMArrayAtIdx :: (Monad m, S.MArray a e m, S.Ix i) => a i e -> i -> (e -> e) -> m ()
modifyMArrayAtIdx a idx f =
    do
      e <- S.readArray a idx
      S.writeArray a idx (f e)

countBackground probes =
    do
      basecount <- S.newSTRef 0
      transcount <- S.newSTRef 0
      bases <- (S.newArray 4 0)::m S.STUArray Int Int
      transitions <- (S.newArray 4 0)::m S.STUArray (Int, Int) Int

      forM probes $ \probe ->
          do
            byteStringForM probe $ \b ->
                do
                  S.modifySTRef basecount (+1)
                  modifyMArrayAtIdx bases b (+1)
            byteStringPairsForM probe $ \bp -> \b ->
                do
                  S.modifySTRef transcount (+1)
                  modifyMArrayAtIdx transitions (bp, b) (+1)

      fbases <- S.unsafeFreeze bases
      ftransitions <- S.unsafeFreeze transitions
      fbasecount <- S.readSTRef basecount
      ftranscount <- S.readSTRef transcount
      return (fbases, fbasecount, ftransitions, ftranscount)

computeBackgroundProbabilities probes =
    let
        (bases, nbases, transitions, ntransitions) = S.runST $ countBackground probes
        baseprobs = (A.amap ((/ (fromIntegral nbases)) . fromIntegral) bases)
        condprobs = A.array $ concatMap (\i -> (map (\j -> ((i,j), (fromIntegral (transitions!(i, j))) / (baseprobs!i))) [0..3])) [0..3]
    in
      (baseprobs, condprobs)

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

data Params = Params {
      priorProb :: Double,
      posteriorCutoff :: Double
}

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
      poPosteriorCutoff :: Maybe Double,
      poPrior :: Maybe Double,
      pofastaFile :: Maybe String,
      pofsaFile :: Maybe String
    }
defaultOptions = ProgramOptions { poshowHelp = False, poPosteriorCutoff = Nothing, poPrior = Nothing,
                                  pofastaFile = Nothing, pofsaFile = Nothing }

options :: [OptDescr (ProgramOptions -> ProgramOptions)]
options =
    [
      Option ['h'] ["help"] (NoArg (\opts -> opts {poshowHelp = True } )) "Displays this help message"
    , Option ['p'] ["posterior-cutoff"] (ReqArg (\val -> \opts -> opts { poPoteriorCutoff = Just (read val) }) "DOUBLE") "Minimum value required for posterior cut-off in order for data to be retained"
    , Option ['0'] ["prior-prob"] (ReqArg (\val -> \opts -> opts { poPrior = Just (read val) }) "DOUBLE") "The prior probability that a given transcription factor is at any particular site"
    , Option ['M'] ["matrix-file"] (ReqArg (\val -> \opts -> opts { pofastaFile = Just val }) "FILE") "File contain TF matrices"
    , Option ['f'] ["probe-file"] (ReqArg (\val -> \opts -> opts { pofsaFile = Just val }) "FILE") "File containing all probe sequences"
    ]
usageString = "Usage: sbasetram [opts...]"

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
            ProgramOptions { poPosteriorCutoff = Nothing } ->
                commandLineInfo "Posterior cutoff not specified"
            ProgramOptions { poPrior = Nothing } ->
                commandLineInfo "Prior probability not specified"
            ProgramOptions { pofastaFile = Nothing } ->
                commandLineInfo "FASTA (matrix) file not specified"
            ProgramOptions { pofsaFile = Nothing } ->
                commandLineInfo "FSA (sequence) file not specified"
            ProgramOptions { pofastaFile = Just fasf, pofsaFile = Just fsaf, poPosteriorCutoff = Just postco, poPrior = Just prior } ->
                gmatimMain fasf fsaf (Params prior postco)
    (_, _, e) -> commandLineInfo $ "Command line parse error:\n" ++ (concat e)
-- "/home/andrew/Documents/TSM/matrices.fasta"
-- "/home/andrew/tgz/yeast_Young_6k.fsa"

gmatimMain fastaFile fsaFile params =
    do
      matricesAsWeights <- loadTRANSFACMatrices fastaFile
      fsaData <- loadFSAData fsaFile
      let
          matrixSummary = summariseMatrices matricesAsWeights
          backgroundModel = buildBackgroundModel fsaData
        in
          forM_ fsaData $ \(probe, seqbs) ->
          let
              allHits = searchForMatrices params seqbs matrixSummary
              -- If we get multiple hits per probe, we simply take the one with the highest probability.
              groupedHits = groupBy (\(x,_)(y,_) -> x==y) $ sortBy (comparing fst) allHits
              bestHits = map (maximumBy comparing fst) groupedHits
            in
              if null matches
              then
                  return []
              else
                  do
                    putStrLn $ showString ">" probe
                    forM bestHits $ \(prob, matrix) ->
                        putStrLn $ (showString matrix . showString " " . shows prob) ""
