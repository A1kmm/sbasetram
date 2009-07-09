import qualified Data.Array.Unboxed as A
import Data.Array.IArray ((!))
import Data.Maybe
import TSMParser
import FSAParser
import Data.List
import Control.Monad
import qualified Data.ByteString as B
import Control.Exception
import Debug.Trace

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

searchForMatrices bs mats =
    let
        l = B.length bs
    in
      concatMap (checkForMatchingMatricesAt bs mats) [0..(l-1)]

checkForMatchingMatricesAt bs mats start =
   mapMaybe (\(i,m) -> (liftM $ const (i, start)) (checkMatrixMatch bs start m)) mats

checkMatrixMatch :: B.ByteString -> Int -> TSMMatrix -> Maybe Double
checkMatrixMatch bs start mat =
    let
        lrem = (B.length bs) - start
        (_, lmat) = A.bounds (matrixInformationVector mat)
    in
      if lmat < lrem
      then
          unsafeCheckMatrixMatch bs start mat lmat
      else
          Nothing

mssCutoff = 0.75
cssCutoff = 0.7

normaliseScore minScore maxScore currentScore =
    (currentScore - minScore) / (maxScore - minScore)

computeScoreOver :: [Int] -> B.ByteString -> Int -> A.UArray (Int, Int) Double -> A.UArray Int Double -> Double
computeScoreOver idxs seqbs offset freqs iv =
    sum $ map (\idx -> (freqs!((fromIntegral $ B.index seqbs (idx + offset)), idx)) * (iv!idx)) idxs

unsafeCheckMatrixMatch :: B.ByteString -> Int -> TSMMatrix -> Int -> Maybe Double
unsafeCheckMatrixMatch bs start mat n =
    let
        css = normaliseScore (coreMinScore mat) (coreMaxScore mat)
                (computeScoreOver (matrixCore mat) bs start (matrixFreqs mat) (matrixInformationVector mat))
        mss = normaliseScore (matrixMinScore mat) (matrixMaxScore mat)
              (computeScoreOver [0..n] bs start (matrixFreqs mat) (matrixInformationVector mat))
    in
      if {- trace ((showString "CSS = " . shows css)"") $ -} css > cssCutoff && mss > mssCutoff
      then
          Just mssCutoff
      else
          Nothing

main =
    do
      matricesAsWeights <- loadTRANSFACMatrices "/home/andrew/Documents/TSM/matrices.fasta"
      fsaData <- loadFSAData "/home/andrew/tgz/yeast_Young_6k.fsa"
      let
          matrixSummary = summariseMatrices matricesAsWeights
        in
          forM_ fsaData $ \(probe, seqbs) ->
          let
              matches = nub (map fst $ searchForMatrices seqbs matrixSummary)
            in
              if null matches
              then
                  return []
              else
                  do
                    putStrLn $ showString ">" probe
                    forM matches $ \matrix ->
                        putStrLn matrix
