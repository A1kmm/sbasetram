module FSAParser (loadFSAData)
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
