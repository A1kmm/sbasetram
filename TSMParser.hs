module TSMParser (transfacMatrixParser, loadTRANSFACMatrices)
where
  import Text.ParserCombinators.Parsec
  import Data.Char
  import qualified Data.Array.Unboxed as A
  import Control.Monad

  loadTRANSFACMatrices :: String -> IO [(String, A.UArray (Int, Int) Double)]
  loadTRANSFACMatrices filename = do
    result <- parseFromFile transfacMatrixParser filename
    case result of
      Left err -> error $ (showString "Cannot parse TRANSFAC file: " . shows err) ""
      Right val -> return val

  buildMatrix matrixRaw =
      A.array ((0,0), (3, (maximum . (map $ snd.fst)) matrixRaw)) matrixRaw

  transfacMatrixParser = processMatrices []
  processMatrices matrices =
      (
       do
         (maybeName, maybeMatrixRaw) <- processMatrix (Nothing, [])
         string "//\n"
         processMatrices $ case (maybeName, maybeMatrixRaw) of
                             (Just name, Just matrixRaw) -> (name, buildMatrix matrixRaw):matrices
                             _ -> matrices
      ) <|> (eof >> return matrices)

  processMatrix :: (Maybe String, [((Int, Int), Double)]) -> Parser (Maybe String, Maybe [((Int,Int),Double)])
  processMatrix x@(name, values) =
      (
       do
         try $ string "BF T"
         manyTill digit (char ';')
         char ' '
         name' <- manyTill (noneOf ";\n") (char ';')
         manyTill (noneOf "\n") (char '\n')
         processMatrix (Just name', values)
      ) <|>
      (
       do
         entry <- try $ do
                    char 'M'
                    processMatrixEntry (Just values)
         case entry of
           Just matrixEntry ->
               processMatrix $ (name, matrixEntry)
           Nothing -> consumeMatrix
      ) <|>
      (
       do
         alphaNum
         alphaNum
         manyTill (noneOf "\n") (char '\n')
         processMatrix x
      ) <|> (return (name, Just values))

  consumeMatrix =
      (
       do
         alphaNum
         alphaNum
         manyTill (noneOf "\n") (char '\n')
         consumeMatrix
       ) <|> (return $ (Nothing, Just []))

  nucleotideCode =
      (char 'A' >> return 0) <|>
      (char 'G' >> return 1) <|>
      (char 'C' >> return 2) <|>
      (char 'T' >> return 3)
  
  processMatrixEntry values =
      do
        base <- nucleotideCode
        char ' '
        processMatrixValues base 0 values
  
  parseIntegerOrIgnoreFloat l =
      do
        (sign, ls) <- parseOptionalSign l
        (whole, lw) <- parseNumber ls 0.0
        case lw of
          0 ->  return $ Just (sign * whole)
          _ ->
              (
               do
                 char '.'
                 parseFracNumber (lw - 1) 0.0 0.1
                 return $ Nothing {- sign * (whole + frac) -}
              ) <|> (return $ Just (sign * whole))
  
  parseOptionalSign l =
      (char '+' >> return (1.0, l - 1)) <|> (char '-' >> (return $ (-1.0, l - 1))) <|> (return (1.0, l))
  
  parseNumber :: Int -> Double -> Parser (Double, Int)
  parseNumber l n =
      case l of
        0 -> return (n, l)
        _ ->
            (
             do
               d <- digit
               parseNumber (l - 1) (n * 10.0 + ((fromIntegral . digitToInt) d))
            ) <|> (return (n, l))
  
  parseFracNumber :: Int -> Double -> Double -> Parser (Double, Int)
  parseFracNumber l s m =
      case l of
        0 -> return (s, 0)
        _ ->
            (
             do
               d <- digit
               parseFracNumber (l - 1) (((fromIntegral . digitToInt) d) * m + s) (m*0.1)
            ) <|> return (s, l)
  
  processMatrixValues base pos values =
      (char '\n' >> return values) <|>
      (
       do
         v <- many $ char ' '
         maybeval <- parseIntegerOrIgnoreFloat (5 - (length v))
         case (values, maybeval) of
           (Just vs, Just val) ->
               processMatrixValues base (pos + 1) (Just (((base, pos), val):vs))
           _ ->
               processMatrixValues base (pos + 1) Nothing
      )
