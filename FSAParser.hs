module FSAParser (loadFSAData, reverseComplement)
where
  import Text.ParserCombinators.Parsec
  import qualified Data.ByteString as B

  loadFSAData :: String -> IO [(String, B.ByteString)]
  loadFSAData filename = do
    result <- parseFromFile fsaParser filename
    case result of
      Left err -> error $ (showString "Cannot parse FSA file: " . shows err) ""
      Right val -> return val

  fsaParser = processSeqs []
  processSeqs seqs =
      (
       do
         char '>'
         name <- manyTill (noneOf "\n") (char '\n')
         seqv <- processSequence []
         processSeqs $ (name, B.concat (reverse seqv)):seqs
      ) <|> (return seqs)

  nucleotideCode =
      (char 'A' >> return 0) <|>
      (char 'G' >> return 1) <|>
      (char 'C' >> return 2) <|>
      (char 'T' >> return 3)

  processSequence seqs =
      (
       do
         seqData <- manyTill nucleotideCode (char '\n')
         processSequence $ (B.pack seqData):seqs
      ) <|> (return seqs)

  reverseComplement bs =
      B.pack $ map (\ i ->
                        case (B.index bs i)
                        of
                          0 -> 3 -- A => T
                          1 -> 2 -- G => C
                          2 -> 1 -- C => G
                          3 -> 0 -- T => A
                          _ -> error "Invalid base passed to reverseComplement"
                   ) [l,(l-1)..0]
          where
            l = (B.length bs) - 1

