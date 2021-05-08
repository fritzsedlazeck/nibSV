import ./kmers
import tables
export tables
import ./svidx

type Read* = object
    ## key of svid, count of supporting kmers
    compatible_SVs*: CountTable[uint32]

proc process_read*(s: string, idx: SvIndex, k: int, space: int): Read =
    # find SVs with kmers intersecting with those from this read.
    var kmers = (cast[Dna](s)).dna_to_kmers(k)
    if(space > 0):
        kmers = spacing_kmer(kmers, space)
    for kmer in kmers.seeds:
        var matching_svs = idx.lookupKmer(kmer.kmer)
        for svId in matching_svs:
            result.compatible_SVs.inc(svId)


proc filter_read_matches*(read: var Read, min_matches: int, winner_takes_all: bool) =
    ## track sv with most kmer matches
    var removables: seq[uint32]
    var max_sv = int.high
    var max_kcnt = 0
    for sv, kcnt in read.compatible_SVs:
        if kcnt < min_matches:
            removables.add(sv)
        if kcnt > max_kcnt:
            max_sv = sv.int
            max_kcnt = kcnt

    if winner_takes_all:
        clear(read.compatible_SVs)
        read.compatible_SVs.inc(max_sv.uint32, max_kcnt)
    else:
        for r in removables:
            read.compatible_SVs.del(r)
