# vim: sw=4 ts=4 sts=4 tw=0 et:
import hts
import kmers
import tables
import svidx

type
    Chunk = object
        chrom_name: string
        chrom_start: int
        chrom_end: int


iterator createdChunks(fai: Fai, chunk_size: int): Chunk =
    for i in 0..<fai.len:
        let chrom_name = fai[i]
        let chrom_len = fai.chrom_len(chrom_name)
        let step = if chunk_size <= 0:
            chrom_len # typically small in this case
        else:
            chunk_size
        for j in countup(0, chrom_len, step):
            yield Chunk(chrom_name: chrom_name, chrom_start: j, chrom_end: j + step)


proc addRefCount(svKmers: var SvIndex, full_sequence: string, kmer_size: int = 25, space: int = 0) =
    ## Use spaced-seeds if space > 0. (Try 50.)
    var convertedKmers: pot_t = dna_to_kmers(full_sequence, kmer_size)
    if space > 0:
        convertedKmers = spacing_kmer(convertedKmers, space)
    #for seed in convertedKmers.seeds:
    #    echo "btd:", bin_to_dna(seed.kmer, convertedKmers.word_size, seed.strand), ' ', seed.kmer

    for km in convertedKmers.seeds:
        if km.kmer in svKmers.counts:
            svKmers.counts[km.kmer].refCount.inc

proc updateChunk(svKmers: var SvIndex, fai: Fai, chunk: Chunk, kmer_size: int, space: int = 0) =
    var sub_seq = fai.get(chunk.chrom_name, chunk.chrom_start, chunk.chrom_end)
    addRefCount(svKmers, sub_seq, kmer_size, space)

proc updateSvIndex*(input_ref_fn: string, svKmers: var SvIndex, kmer_size: int = 25, chunk_size: int = 1_000_000, space: int = 0) =
    ## Walk over reference sequences and count kmers.
    ## Update any existing svIdx entries with these counts.
    ## Use spaced-seeds if space > 0. (Try 50.)
    var fai: Fai
    if not fai.open(input_ref_fn):
        quit "couldn't open fasta"

    for i in createdChunks(fai, chunk_size):
        echo " chunk i=", i
        updateChunk(svKmers, fai, i, kmer_size, space)

when isMainModule:
    import hts
    var fai: Fai
    import times

    if not fai.open("/data/human/g1k_v37_decoy.fa"):
        quit "bad"

    var s = fai.get("22")
    var svkmers: svIdx
    new(svkmers)
    echo "starting"
    for i in countup(0, 100_000_000, 10):
        svkmers[i.uint64] = (0'u32, 0'u32, newSeq[uint32]())

    var t0 = cpuTime()
    svKmers.addRefCount(s)
    echo "time:", cpuTime() - t0
