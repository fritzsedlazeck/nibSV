# vim: sw=4 ts=4 sts=4 tw=0 et:
import refmers
import svidx
import strformat
import classify
import reporter
#from ./read import `$`
from os import nil
from tables import len

proc main_runner*(variants_fn, refSeq_fn, reads_fn: string, prefix = "test", kmer_size: int = 29, space: int = 0, maxRefKmerCount : uint32 = 0 ) =
    ## Generate a SV kmer database, and genotype new samples.
    ## If a file called "{prefix}.sv_kmers.msgpack" exists, use it.
    ## Otherwise, generate it.
    var idx: SvIndex

    # flank is just 1 less than kmer size.
    # for spaced, it's kmer_size + kmersize + spaces - 1
    let flank = kmer_size + int(space > 0) * (space + kmer_size) - 1
    var index_fn = "{prefix}.{kmer_size}.{space}.sv_kmers.msgpck".fmt

    if not os.fileExists(index_fn):
        echo "building an SV kmer DB."
        idx = buildSvIndex(refSeq_fn, variants_fn, flank, kmer_size, space)
        echo "updating reference kmer counts."
        updateSvIndex(refSeq_fn, idx, kmer_size, 10_000_000, space)
        echo "dumpIndexToFile:'", index_fn, "'"
        dumpIndexToFile(idx, index_fn)
    else:
        echo "loadIndexFromFile:'", index_fn, "'"
        idx = loadIndexFromFile(index_fn, kmer_size)

    echo "final idx contains: {idx.len} forward and reverse SV kmers.".fmt


    filterRefKmers(idx, maxRefKmerCount)


    #echo dumpIndexToJson(idx)
    let classifyCount = classify_file(reads_fn, idx, kmer_size, space, refseq_fn)

    #echo "classifyCount:"
    #echo classifyCount


    echo "reporting variants."

    report(variants_fn, classifyCount, idx, prefix)

    echo "nibbleSV finished without problems, goodbye!"


when isMainModule:
  import cligen
  dispatch(main_runner)
