# vim: sw=4 ts=4 sts=4 tw=0 et:
import deques
import tables
#from sets import nil
from algorithm import sort
from hashes import nil
from strutils import format
from ./util import raiseEx, PbError

export PbError

type
    Dna* = string # someday, this might be an array
    Bin* = uint64 # compact bitvector of DNA
    ##  In bitvector, A is 0, C is 1, G is two, and T is 3.

    Min* = uint64 # minimizer
    Strand* = enum
        forward, reverse

    ##  kmer - a uint64 supporting a maximum of 32 DNA bases.
    ##  pos  - position along the sequence
    seed_t* = object
        kmer*: Bin
        pos*: uint32
        strand*: Strand

    minimizer_t* = object
        minimizer*: Min
        pos*: uint32
        strand*: Strand

    ##  a & b are two seed_t's designed for matching in the hash lookup
    seed_pair_t* = object
        a*: seed_t
        b*: seed_t

    Hash* = int

    ##  seeds - a pointer to the kmers
    pot_t* = ref object of RootObj
        word_size*: uint8 # <=32
        seeds*: seq[seed_t]

    ##  searchable seed-pot
    spot_t* = ref object of pot_t
        ht*: tables.TableRef[Bin, int]

var seq_nt4_table: array[256, int] = [
        0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]


##  @return uninitialized
#
proc newDna(size: int): Dna =
    return newString(size)

# hashes for sets and tables

proc hash*(s: kmers.seed_t): hashes.Hash =
    #hashes.hash(s.pos) + hashes.hash(s.kmer shl 8) + hashes.hash(s.strand)
    hashes.hash([s.pos.int64, s.kmer.int64, s.strand.int64])

proc hash*(p: kmers.seed_pair_t): hashes.Hash =
    hashes.hash([hash(p.a), hash(p.b)])

# convenience for C coders

template `<<`(a, b: uint64): uint64 =
    a shl b

template `>>`(a, b: uint64): uint64 =
    a shr b

## Return binary version of kmer exactly matching Dna.
## Mostly for testing, as this is not efficient for sliding windows.
## NOT FINISHED.
#
proc encode(sq: Dna) =
    assert sq.len() <= 32
    let k = sq.len()
    let
        shift1: uint64 = 2'u64 * (k - 1).uint64
        mask: uint64 = (1'u64 << (2 * k).uint64) - 1
    var
        forward_bin: Bin = 0
        reverse_bin: Bin = 0
    for i in 0 ..< k:
        let ch = cast[uint8](sq[i])
        let c = seq_nt4_table[ch].uint64
        assert c < 4
        forward_bin = (forward_bin << 2 or c) and mask
        reverse_bin = (reverse_bin >> 2) or (
                3'u64 xor c) << shift1

##  Converts a char * into a set of seed_t objects.
##  @param  sq  - sequence
##  @param  k - kmer size (<=32)
##  @return pot
#
proc dna_to_kmers*(sq: Dna; k: int): pot_t =
    if k > 32:
        raiseEx("k > 32")

    let
        shift1: uint64 = 2'u64 * (k - 1).uint64
        mask: uint64 = (1'u64 << (2 * k).uint64) - 1
    #echo format("shift1=$# mask=$#", shift1, mask)

    var forward_kmer: seed_t
    var reverse_kmer: seed_t

    forward_kmer.kmer = 0
    forward_kmer.pos = 0
    reverse_kmer.kmer = 0
    reverse_kmer.pos = 0
    forward_kmer.strand = forward
    reverse_kmer.strand = reverse

    var kmers: pot_t
    new(kmers)
    kmers.seeds = newSeqOfCap[seed_t](max(sq.len - int(k) + 1,0))
    kmers.word_size = k.uint8

    ##  lk is the length of the kmers being built on the fly. The variable n is the total number of
    var
        i: int
        lk: int
        n: int
    i = 0
    lk = 0
    n = 0

    while i < sq.len():
        let ch = cast[uint8](sq[i])
        let c = seq_nt4_table[ch].uint64
        if c < 4:
            forward_kmer.kmer = (forward_kmer.kmer << 2 or c) and mask
            reverse_kmer.kmer = (reverse_kmer.kmer >> 2) or (
                    3'u64 xor c) << shift1
            #echo format("[$#]=$# $#==$#($# $#) f:$# r:$#",
            #    i, sq[i], ch, c, (3'u8 xor c), (3'u8 xor c).uint64 shl shift1, forward_kmer.kmer, reverse_kmer.kmer)
            inc(lk)
        else:
            ##  advance the window beyond the unknown character
            lk = 0
            inc(i, k)
            inc(forward_kmer.pos, k)
            forward_kmer.kmer = 0
            inc(reverse_kmer.pos, k)
            reverse_kmer.kmer = 0

        if lk >= k:
            inc(n, 2)
            kmers.seeds.add(forward_kmer)
            kmers.seeds.add(reverse_kmer)
            inc(forward_kmer.pos, 1)
            inc(reverse_kmer.pos, 1)
        inc(i)




    return kmers

##  A function to convert the binary DNA back into character
##  @param kmer   up to 32 2-bit bases
##  @param k      kmer length
##  @param strand If reverse, start at kth bit and go backwards.
##
##  Zero is A, one is C, G is two, and T is 3
#
proc bin_to_dna*(kmer: Bin; k: uint8; strand: Strand = forward): Dna =
    var lookup: array[4, char] = ['A', 'C', 'G', 'T']
    var mask: uint64 = 3
    var i: uint8 = 0
    var tmp: uint64 = 0
    var offset: uint64 = 0

    var dna = newDna(k.int)
    i = 0
    while i < k:
        if strand == forward:
            offset = (k - i - 1) * 2
            tmp = kmer >> offset
            dna[i.int] = lookup[mask and tmp]
        else:
            offset = i * 2
            tmp = kmer >> offset
            dna[i.int] = lookup[mask and not tmp]
        inc(i)

    return dna

proc nkmers*(pot: pot_t): int =
    return len(pot.seeds)

##  Prints the pot structure to STDOUT
##  @param pot a ref to the pot
#
proc print_pot*(pot: pot_t) =
    var i: int = 0

    while i < pot.seeds.len():
        let dna = bin_to_dna(pot.seeds[i].kmer, pot.word_size,
                         pot.seeds[i].strand)
        echo format("pos:$# strand:$# seq:$# bin:$#",
            pot.seeds[i].pos, pot.seeds[i].strand, dna, pot.seeds[i].kmer)
        inc(i, 1)

proc get_dnas*(pot: pot_t): seq[Dna] =
    for i in 0 ..< pot.seeds.len():
        let dna = bin_to_dna(pot.seeds[i].kmer, pot.word_size,
                         pot.seeds[i].strand)
        result.add(dna)

proc cmp_seeds(a, b: seed_t): int =
    let c = a.kmer
    let d = b.kmer

    if c < d:
        return -1

    if c == d:
        if a.pos < b.pos:
            return -1
        else:
            return 0

    return 1

# Actual implementation, private.
#
proc make_searchable(seeds: var seq[seed_t]; ht: var tables.TableRef[Bin, int]) =
    seeds.sort(cmp_seeds)
    ht = newTable[Bin, int]()
    #let dups = sets.initHashSet[Bin]()
    var ndups = 0

    var i: int = 0
    while i < seeds.len():
        let key = seeds[i].kmer
        if ht.hasKeyOrPut(key, i):
            ndups += 1
            #echo format("WARNING: Duplicate seed $# @$#, not re-adding @$#",
            #        key, i, ht[key])
        inc(i)
#[
    if ndups > 0:
        echo format("WARNING: $# duplicates in kmer table", ndups)

]#

##  Construct searchable-pot from pot.
##  Move construct seeds (i.e. original is emptied).
##
##  Sort the seeds and load the kmers into a hash table.
##  For any dups, the table refers to the first seed with that kmer.
#
proc initSpot*(kms: var pot_t): spot_t =
    new(result)
    result.word_size = kms.word_size
    shallowCopy(result.seeds, kms.seeds)
    #kms.seeds = @[]
    kms = nil # simpler, obvious move-construction
    make_searchable(result.seeds, result.ht)

##  Check for the presence or absence of a kmer in a
##  pot regardless of the position.
##  @param  pot_t * - a pointer to a pot_t
##  @return false if kmer doesn't exist
#
proc haskmer*(target: spot_t; query: Bin): bool =
    if target.ht.hasKey(query):
        return true
    return false

## Counts the number of shared unique kmers
## @param pot_t * - a pointer to a pot_t
## @param pot_t * - a pointer to a pot_t
## @return int - number of shared kmers
#
proc uniqueShared*(a, b: spot_t): int =
    result = 0

    for k in a.ht.keys():
        if(haskmer(b, k)):
            inc(result)

## Find (target - remove), without altering target.
#
proc difference*(target: pot_t; remove: spot_t): pot_t =
    new(result)
    result.word_size = target.word_size

    var kmer_stack = newSeq[seed_t]()

    for i in 0 ..< target.seeds.len():
        if(not haskmer(remove, target.seeds[i].kmer)):
            kmer_stack.add(target.seeds[i])

    result.seeds = kmer_stack

## Return the seeds in the intersection of target and query.
#
proc search*(target: spot_t; query: pot_t): deques.Deque[seed_pair_t] =
    echo format("Searching through $# kmers", query.seeds.len())
    var hit_stack = deques.initDeque[seed_pair_t](128)
    var hit: seed_pair_t
    var hit_index: int

    var i: int = 0
    #echo format("target.ht=$#", target.ht)
    #echo format("query.ht=$#", query.ht)
    while i < query.seeds.len():
        let key = query.seeds[i].kmer
        if key in target.ht:
            hit_index = target.ht[key]
            #echo format("For $# ($#), ql=$# tl=$#, hit_index=$#", i, key, query.seeds.len(), target.seeds.len(), hit_index)
            while (hit_index < target.seeds.len() and key == target.seeds[
                    hit_index].kmer):
                #echo format("--For $# ($#), ql=$# tl=$#, hit_index=$#", i, key, query.seeds.len(), target.seeds.len(), hit_index)
                hit.a = query.seeds[i]
                hit.b = target.seeds[hit_index]
                deques.addLast(hit_stack, hit)
                inc(hit_index, 1)
        inc(i)

    return hit_stack

## This function counts the number of uniq kmers in the pot if searchable if not
## the function calls make searchable.
## @param pot_t - a ref to a pot_t
## TODO: add test coverage
#
proc nuniq*(pot: spot_t): int =
    return len(pot.ht)

proc spacing_kmer*(pot: pot_t; space: int): pot_t =
    #doAssert(space > int(pot.word_size)) # typical, but not necessary
    doAssert(pot.word_size <= 16)

    new(result) #default return knwos the type from function header
    result.word_size = pot.word_size*2

    for i in (0 ..< pot.seeds.len - 2*(space + pot.word_size.int)):
        let j = i + 2*(space + pot.word_size.int)
        assert(j < pot.seeds.len)

        let k1 = pot.seeds[i]
        let k2 = pot.seeds[j]
        assert(k1.strand == k2.strand)

        # new kmer
        var k: seed_t
        let (left, right) = if k1.strand == forward:
            (k1, k2)
        else:
            (k2, k1)
        k.kmer = left.kmer
        k.kmer = k.kmer << 2*pot.word_size
        k.kmer = k.kmer or right.kmer
        k.strand = left.strand
        k.pos = left.pos
        result.seeds.add(k)
