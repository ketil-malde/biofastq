{- |
   Module: Bio.Sequence.FastQ

   Support the FastQ format that combines sequence and quality. See:

   * <http://www.bioperl.org/wiki/FASTQ_sequence_format>

   Of course, this is yet another vaguely defined pseudo-standard with
   conflicting definitions.  Of course Solexa had to go and invent not one, but two 
   different, and indistinguishably so, ways to do it:

   * <http://www.bcgsc.ca/pipermail/ssrformat/2007-March/000137.html>

   * <http://maq.sourceforge.net/fastq.shtml>

   * <http://en.wikipedia.org/wiki/FASTQ_format>

   Sanger-style FastQ-format is supported with the (h)read/writeSangerQ functions,
   and the new Illumina/Solexa-style with (h)read/writeIllumina.

   As far as I know, FastQ is only used for nucleotide sequences, never amino acid.
-}


module Bio.Sequence.FastQ
    (
     -- * Reading FastQ
    readFastQ, hReadFastQ, parse
     -- * Writing FastQ
    , writeFastQ, hWriteFastQ

    -- * use Sanger-style quality information
    , readSangerQ, hReadSangerQ
    , writeSangerQ, hWriteSangerQ
    
    -- * use Illumina (>v1.3)-style quality information
    , readIllumina, hReadIllumina
    , writeIllumina, hWriteIllumina
    
    -- * Sequence data structure
    , Sequence(..)
    ) where

import System.IO
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.ByteString.Lazy as BB
import Data.List (unfoldr)
import Data.Maybe (fromJust)
import Bio.Core.Sequence

data Sequence = Seq SeqLabel SeqData QualData deriving Show

instance BioSeq Sequence where
  seqid (Seq sl _ _) = sl
  seqheader (Seq sl _ _) = sl
  seqdata (Seq _ sd _) = sd  
  seqlength = Offset . B.length . unSD . seqdata -- should be default
  
instance BioSeqQual Sequence where
  seqqual (Seq _ _ q) = q
  
{-# DEPRECATED readFastQ, hReadFastQ, writeFastQ, hWriteFastQ  "FastQ assumes Sanger-style quality info use {read,write}SangerQ or -Illumina instead" #-}

readSangerQ, readIllumina :: FilePath -> IO [Sequence]
readSangerQ = readFastQ
readIllumina f = addQual (negate 31) `fmap` readFastQ f

hReadSangerQ, hReadIllumina :: Handle -> IO [Sequence]
hReadSangerQ = hReadFastQ
hReadIllumina h = addQual (negate 31) `fmap` hReadFastQ h

writeSangerQ, writeIllumina :: FilePath -> [Sequence] -> IO ()
writeSangerQ = writeFastQ
writeIllumina f = writeFastQ f . addQual 31 

hWriteSangerQ, hWriteIllumina :: Handle -> [Sequence] -> IO ()
hWriteSangerQ = hWriteFastQ
hWriteIllumina h = hWriteFastQ h . addQual 31

addQual :: Qual -> [Sequence] -> [Sequence]
addQual (Qual q) = map (\(Seq h d mq) -> (Seq h d $ qmap (+q) mq))
  where qmap f (QualData qd) = QualData (BB.map f qd)

readFastQ :: FilePath -> IO [Sequence]
readFastQ f = (go . B.lines) `fmap` B.readFile f

hReadFastQ :: Handle -> IO [Sequence]
hReadFastQ h = (go . B.lines) `fmap` B.hGetContents h

go :: [B.ByteString] -> [Sequence]
go = map (either error id) . unfoldr parse

-- | Parse one FastQ entry, suitable for using in 'unfoldr' over
--   'B.lines' from a file
parse :: [B.ByteString] -> Maybe (Either String (Sequence), [B.ByteString])
parse (h1:sd:h2:sq:rest) =
    case (B.uncons h1,B.uncons h2) of
      -- The fast path: four-line format
      (Just ('@',h1name), Just ('+',h2name))
          | h1name == h2name || B.null h2name
            -> Just (Right $ Seq (SeqLabel h1name) (SeqData sd) (QualData (BB.map (subtract 33) sq)), rest)
          | otherwise
            -> Just (Left $ "Bio.Sequence.FastQ: name mismatch:" ++ showStanza, rest)
      (Just ('@',h1name), Just (_,_)) -- not the '+' quality start header
        -> let ls = (sd:h2:sq:rest)
               ss = takeWhile (\l -> fst (fromJust (B.uncons l)) /= '+') ls
               qs = take (length ss) $ drop (length ss+1)     ls
               rs = drop (length ss*2+1)                      ls
               s' = B.concat ss
               q' = B.concat qs
               in if B.length s' /= B.length q' 
                  then Just (Left $ "Bio.Sequence.FastQ: lenght of sequence data doesn't match qualdata for "++B.unpack h1name,rs)
                  else Just (Right $ Seq (SeqLabel h1name) (SeqData s') (QualData (BB.map (subtract 33) q')), rs)
      _     -> Just (Left $ "Bio.Sequence.FastQ: illegal FastQ format:" ++ showStanza, rest)
    where showStanza = unlines $ map B.unpack [ h1, sd, h2, sq ]
parse [] = Nothing
parse fs = let showStanza = unlines (map B.unpack fs)
               err = Left $ "Bio.Sequence.FastQ: illegal number of lines in FastQ format: " ++ showStanza
           in Just (err, [])

writeFastQ :: FilePath -> [Sequence] -> IO ()
writeFastQ f = B.writeFile f . B.concat . map toFastQ

hWriteFastQ :: Handle -> [Sequence] -> IO ()
hWriteFastQ h = B.hPut h . B.concat . map toFastQ
