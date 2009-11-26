{-# LANGUAGE RankNTypes,BangPatterns #-}
import Prelude hiding (sum, map, (++), concatMap, concatMap, concat, null, head, tail)
import qualified Data.Array.Unboxed as A
import Data.Array.IArray ((!))
import Data.Maybe
import TSMParser
import FSAParser
import Data.List.Stream
import Control.Monad
import qualified Data.ByteString as B
import System.Console.GetOpt
import System.Environment
import Numeric.GSL.Special.Gamma
import qualified Data.Array.ST as S
import qualified Control.Monad.ST.Strict as S
import qualified Data.STRef.Strict as S
import Data.Word
import Data.Ord

type TSMMatrix = A.UArray (Word8, Int) Double

prependProbsOneEntry freqs n i l j =
    let
        f = freqs!(fromIntegral j, i)
        logp = log ((beta (f + 2.0) (n - f + 1.0)) / (beta (f + 1.0) (n - f + 1.0)))
    in
      ((j, i), logp):l

prependProbsOnePosition freqs l i =
    let
        n = sum $ map (\j -> freqs!(j, i)) [0..3]
    in
      foldl' (prependProbsOneEntry freqs n i) l [0..3]

computeProbs freqs = 
    let
        (_, (_, n)) = A.bounds freqs
    in
      A.array ((0, 0), (3, n)) $ (foldl' (prependProbsOnePosition freqs) [] [0..n])

summariseMatrix freqs = computeProbs freqs

summariseMatrices = map (\(n, m) -> (n, summariseMatrix m))

byteStringForM :: Monad m => B.ByteString -> (Word8 -> m ()) -> m ()
byteStringForM s f = B.foldl (\m -> \w -> m >> f w) (return ()) s
byteStringPairsForM :: Monad m => B.ByteString -> (Word8 -> Word8 -> m ()) -> m ()
byteStringPairsForM s f =
    case B.uncons s of
      Nothing -> return ()
      Just (w1, s1) ->
          snd $ B.foldl (\(wp, m) w -> (w, m >> f wp w)) (w1, return ()) s1

modifyMArrayAtIdx :: (Monad m, S.MArray a e m, S.Ix i) => a i e -> i -> (e -> e) -> m ()
modifyMArrayAtIdx a idx f =
    do
      e <- S.readArray a idx
      S.writeArray a idx (f e)

countBackground :: [(String, B.ByteString)] -> forall s . S.ST s (A.UArray Word8 Int, Int, A.UArray (Word8, Word8) Int, Int)

strictIncrement !x = x + 1

modifySTRef' r f =
    do
      lastv <- S.readSTRef r
      let !nextv = f lastv
      S.writeSTRef r nextv

countBackground probes =
    do
      basecount <- S.newSTRef 0
      transcount <- S.newSTRef 0
      bases <- (S.newArray (0, 3) 0) ::S.ST s (S.STUArray s Word8 Int)
      transitions <- (S.newArray ((0, 0), (3, 3)) 0)::S.ST s (S.STUArray s (Word8, Word8) Int)

      forM probes $ \(_, probe) ->
          do
            byteStringForM probe $ \b ->
                do
                  modifySTRef' basecount (+1)
                  modifyMArrayAtIdx bases (fromIntegral b) (+1)
            byteStringPairsForM probe $ \bp -> \b ->
                do
                  modifySTRef' transcount strictIncrement
                  modifyMArrayAtIdx transitions ((fromIntegral bp), (fromIntegral b)) (+1)

      fbases <- S.unsafeFreeze bases
      ftransitions <- S.unsafeFreeze transitions
      fbasecount <- S.readSTRef basecount
      ftranscount <- S.readSTRef transcount
      return (fbases, fbasecount, ftransitions, ftranscount)

computeBackgroundProbabilities :: [(String, B.ByteString)] -> (A.UArray Word8 Double, A.UArray (Word8, Word8) Double)
computeBackgroundProbabilities probes =
    let
        (bases, nbases, transitions, ntransitions) = S.runST $ countBackground probes
        baseprobs = (A.amap (log . (/ (fromIntegral nbases)) . fromIntegral) bases)
        condprobs = A.array ((0, 0), (3,3)) $
                      foldl' (\l i ->
                                  (foldl' (\l' j ->
                                            ((i,j),
                                             (
                                              log ((fromIntegral (transitions!(i, j))) /
                                                   (fromIntegral ntransitions)) -
                                              (baseprobs!i)
                                             )
                                            ):l'
                                          ) l [0..3]
                                  )
                             ) [] [0..3]
    in
      (baseprobs, condprobs)

buildPartialSumsForSeq :: A.UArray (Word8, Word8) Double -> B.ByteString -> A.UArray Int Double
buildPartialSumsForSeq bgcond seqn =
    let
        !n = B.length seqn
    in
      A.array (0, n-1) $ snd (
                              foldl' (\(!s, !l) !i ->
                                          let !s' = s + (bgcond ! ((B.index seqn (i-1)), (B.index seqn i)))
                                          in (s', (i, s'):l)
                                     ) (0, [(0, 0)]) [1..(n-1)]
                             )

 
buildPartialSums0OrderForSeq :: A.UArray Word8 Double -> B.ByteString -> A.UArray Int Double
buildPartialSums0OrderForSeq bgprob seqn =
    let
        !n = B.length seqn
    in
      A.array (0, n) $ snd (
                            foldl' (\(!s, !l) !i ->
                                        let !s' = s + bgprob!(B.index seqn i)
                                        in (s', (i + 1, s'):l)
                                   ) (0, [(0,0)]) [0..(n-1)]
                           )


searchForMatrices :: Params -> (A.UArray Word8 Double, A.UArray (Word8, Word8) Double) -> B.ByteString -> [(String, TSMMatrix)] -> [(String, Double)]
searchForMatrices params (baseprobs, condprobs) bs mats =
    let
        l = B.length bs
        bgsums = {- buildPartialSums0OrderForSeq baseprobs bs -} buildPartialSumsForSeq condprobs bs
    in
      concatMap (checkForMatchingMatricesAt params baseprobs bgsums bs mats) [0..(l-1)]

checkForMatchingMatricesAt :: Params -> A.UArray Word8 Double -> A.UArray Int Double -> B.ByteString -> [(String, TSMMatrix)] -> Int -> [(String, Double)]
checkForMatchingMatricesAt !params !baseprobs !bgsums !bs !mats !start =
   mapMaybe (\(i,m) -> liftM ((,)i) (checkMatrixMatch params bs baseprobs bgsums start m)) mats

checkMatrixMatch :: Params -> B.ByteString -> A.UArray Word8 Double -> A.UArray Int Double -> Int -> TSMMatrix -> Maybe Double
checkMatrixMatch !params !bs !baseprobs !bgsums !start !mat =
    let
        lrem = (B.length bs) - start
        (_, (_, lmat)) = A.bounds mat
    in
      if lmat < lrem
      then
          unsafeCheckMatrixMatch params bs baseprobs bgsums start mat lmat
      else
          Nothing

logplus a b =
    let
        c = max a b
    in
      c + (log ((exp (a - c)) + (exp (b - c))))

-- Computes the probability that at least one of two independent events occurs.
logProbAtLeastOne logp1 logp2 =
    let
        c = max logp1 logp2
    in
      c + (log (exp (logp1 - c)) + (exp (logp2 - c)) - (exp (logp1 + logp2 - c)))

unsafeCheckMatrixMatch !params !bs !baseprobs !bgsums !start !mat !lmat =
    let
        !pdgivenh1 = sum (map (\ !i -> mat ! (B.index bs (i + start), i)) [0..lmat])
        !pdgivenh0 = bgsums!(start + lmat) - bgsums!start + baseprobs!(B.index bs start) -- bgsums!(start + lmat + 1) - bgsums!start -- (log 0.25) * (fromIntegral $ lmat + 1)
        !pdandh1 = pdgivenh1 + (logPriorProb params)
        !pdandh0 = pdgivenh0 + (logPriorNotProb params)
        !posterior = pdandh1 - (pdandh1 `logplus` pdandh0)
    in
      if posterior >= (logPosteriorCutoff params)
      then
         Just posterior
      else
          Nothing

data Params = Params {
      logPriorProb :: !Double,
      logPosteriorCutoff :: !Double,
      logPriorNotProb :: !Double
}

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
    , Option ['p'] ["posterior-cutoff"] (ReqArg (\val -> \opts -> opts { poPosteriorCutoff = Just (read val) }) "DOUBLE") "Minimum value required for posterior probability in order for data to be retained"
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
            opts = foldl' (\opt -> \f -> f opt) defaultOptions o
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
                sbasetramMain fasf fsaf (Params (log prior) (log postco) (log (1.0 - prior)))
    (_, _, e) -> commandLineInfo $ "Command line parse error:\n" ++ (concat e)
-- "/home/andrew/Documents/TSM/matrices.fasta"
-- "/home/andrew/tgz/yeast_Young_6k.fsa"

showsNegative v =
    if v < 0 then
        shows v
    else
        showString "-0"

sbasetramMain fastaFile fsaFile params =
    do
      matricesAsWeights <- loadTRANSFACMatrices fastaFile
      fsaData <- loadFSAData fsaFile
      let
          matrixSummary = summariseMatrices matricesAsWeights
          backgroundModel = computeBackgroundProbabilities fsaData
        in
          forM_ fsaData $ \(probe, seqbs) ->
          let
              allHits = (searchForMatrices params backgroundModel seqbs matrixSummary) ++
                        (searchForMatrices params backgroundModel (reverseComplement seqbs) matrixSummary)
              -- If we get multiple hits per probe, we simply take the one with the highest probability.
              groupedHits = groupBy (\(x,_)(y,_) -> x==y) $ sortBy (comparing fst) allHits
              bestHits = map (maximumBy (comparing snd)) groupedHits
              -- bestHits = map (\h -> (
              --                        (fst . head) h,
              --                         foldl' (\p1 (_,p2) -> logProbAtLeastOne p1 p2) ((snd . head) h) (tail h)
              --                       )
              --                ) groupedHits
            in
              if null bestHits
              then
                  return ()
              else
                  do
                    putStrLn $ showString ">" probe
                    forM_ bestHits $ \(matrix, prob) ->
                        putStrLn $ (showString matrix . showString " " . showsNegative prob) ""
