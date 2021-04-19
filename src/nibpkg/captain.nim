# vim: sw=4 ts=4 sts=4 tw=0 et:
import refmers
import svidx
import strformat
import classify
import reporter
#from ./read import `$`
from os import nil
from tables import len

proc main_runner*(variants_fn, refSeq_fn, reads_fn: string, prefix = "test", kmer_size: int = 25, spaced_seeds : bool = false, space: int = 0, flank: int = 100, maxRefKmerCount : uint32 = 1 ) =
    ## Generate a SV kmer database, and genotype new samples.
    ## If a file called "{prefix}.sv_kmers.msgpack" exists, use it.
    ## Otherwise, generate it.
    var index_fn = "{prefix}.sv_kmers.msgpck".fmt
    var idx: SvIndex

    if not os.fileExists(index_fn):
        echo "building an SV kmer DB."
        let sp = if spaced_seeds:
          space
        else:
          0
        idx = buildSvIndex(refSeq_fn, variants_fn, flank, kmer_size, sp)
        echo "updating reference kmer counts."
        updateSvIndex(refSeq_fn, idx, kmer_size, 1000000, sp)
        echo "dumpIndexToFile:'", index_fn, "'"
        dumpIndexToFile(idx, index_fn)
    else:
        echo "loadIndexFromFile:'", index_fn, "'"
        idx = loadIndexFromFile(index_fn, kmer_size)

    echo "final idx contains: {idx.len} forward and reverse SV kmers.".fmt


    filterRefKmers(idx, maxRefKmerCount)


    #echo dumpIndexToJson(idx)


    let classifyCount = classify_file(reads_fn, idx, kmer_size, spaced_seeds, space)

    #echo "classifyCount:"
    #echo classifyCount


    echo "reporting variants."

    report(variants_fn, classifyCount, idx, prefix)

    echo "nibbleSV finished without problems, goodbye!"


when isMainModule:
  import cligen
  dispatch(main_runner)
