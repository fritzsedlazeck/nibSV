# vim: sw=4 ts=4 sts=4 tw=0 et:
import tables
from strutils import nil
from strformat import fmt
import msgpack4nim, streams, json
import ./kmers

type
    #SvValue* = tuple[refCount: uint32, altCount: uint32, svs: seq[uint32]]
    SvValue* = object
        refCount*: uint32
        altCount*: uint32
        svs*: seq[uint32]

    ## A map from KMER ID -> (number of time kmer appears in a ref seq, number of times kmer appears in an alt seq, list(SVs) that kmer is contained in )
    #svIdx* = TableRef[uint64, SvValue]
    SvIndex* = object
        counts*: Table[uint64, SvValue]
        kmerSize*: uint8

proc len*(idx: SvIndex): int =
    return idx.counts.len


#Cost savings on allocations?
var empty: seq[uint32]

proc lookupKmer*(idx: SvIndex, kmer: seed_t): seq[uint32] {.noInit.} =
    if kmer.kmer in idx.counts:
        return idx.counts[kmer.kmer].svs
    return empty

proc dumpIndexToFile*(idx: SvIndex, fn: string) =
    let strm = openFileStream(fn, fmWrite)
    strm.pack(idx)
    strm.close()

proc loadIndexFromFile*(fn: string, kmerSize: int): SvIndex =
    let strm = openFileStream(fn, fmRead)
    strm.unpack(result)
    strm.close()
    if kmerSize != result.kmerSize.int:
        echo "ERROR: Inconsistent SvIndex file '{fn}'\nkmerSize={kmerSize} != SvIndex.kmerSize={result.kmerSize}".fmt
        doAssert(kmerSize == result.kmerSize.int)

proc `%`(idx: SvIndex): JsonNode =
    result = json.newJObject()
    result["kmerSize"] = %idx.kmerSize
    result["counts"] = json.newJObject()
    for k, v in idx.counts.pairs():
        let val = SvValue(refCount: v.refCount, altCount: v.altCount, svs: v.svs)
        result["counts"][$k] = %val

proc dumpIndexToJson*(idx: SvIndex): string =
    return json.pretty(%idx)

proc loadIndexFromJson*(js: string): SvIndex =
    ## This painful method might become simple if SvIndex values
    ## switched from tuple to object.
    let j = json.parseJson(js)
    result.kmerSize = j["kmerSize"].getInt().uint8
    for key, val in j["counts"]:
        let k: uint64 = strutils.parseBiggestUint(key)
        let v = json.to(val, SvValue)
        result.counts[k] = v

proc insert*(idx: var SvIndex, sequence: string, k: int, sv_idx: int = -1, space: int = 0) =
    ## when inserting reference sequences leave sv_idx as -1
    #doAssert(k == idx.kmerSize.int);
    var l = Dna(sequence).dna_to_kmers(k.int)
    if space > 0:
       l = spacing_kmer(l, space)

    # inserting alternates
    if sv_idx >= 0:
        for kmer in l.seeds:
            var kc = idx.counts.getOrDefault(kmer.kmer)
            kc.altCount.inc
            kc.svs.add(sv_idx.uint32)
            idx.counts[kmer.kmer] = kc

        return

    # inserting reference counts iff the kmer was already found as alternate.
    for kmer in l.seeds:
        # note: sometimes doing double lookup.
        if kmer.kmer notin idx.counts: continue
        idx.counts[kmer.kmer].refCount.inc

proc filterRefKmers*(svKmers: var SvIndex, maxRefCount: uint32) =
    ## Remove entries in the SV index that have a ref count higher than specified
    echo "before:", svKmers.len, " maxRefCount:", maxRefCount
    var toRemove: seq[uint64]
    for k, v in pairs(svKmers.counts):
        if v.refCount > maxRefCount:
            toRemove.add(k)
    for k in toRemove:
        svKmers.counts.del(k)
    echo "after:", svKmers.len, " maxRefCount:", maxRefCount
