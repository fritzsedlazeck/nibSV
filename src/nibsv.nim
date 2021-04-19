from nibpkg/compose import nil
from nibpkg/classify import nil
from nibpkg/captain import nil

when isMainModule:
  import cligen
  dispatchMulti(
        [compose.compose_variants, cmdName = "compose"],
        [classify.buildSvIndex, cmdName = "lookup"],
        [classify.main_classify, cmdName = "classify"],
        [captain.main_runner, cmdName = "main",
        help={
        "variants-fn": "long read VCF SV calls",
        "refSeq-fn": "reference genome FASTA, compressed OK",
        "reads-fn": "input short-reads in BAM/SAM/CRAM/FASTQ",
        "prefix"  : "output prefix",
        "kmer-size" : "kmer size, for spaced seeds use <=16 otherwise <=32",
        "spaced-seeds" : "turn on spaced seeds",
        "space" : "width between spaced kmers",
        "flank" : "number of bases on either side of ALT/REF in VCF records",
        "max-ref-kmer-count" : "max number of reference kmers allowed in SV event"
        }
        ],
  )
