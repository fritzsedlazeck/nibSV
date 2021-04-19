# vim: sw=4 ts=4 sts=4 tw=0 et:
from nibpkg/kmers import hash # avoiding "*" imports
import unittest
import deques
import sequtils
import sets
#from strformat import fmt

suite "kmers":
    test "bin_to_dna":
        check kmers.bin_to_dna(0, 1, kmers.forward) == "A"
        check kmers.bin_to_dna(1, 1, kmers.forward) == "C"
        check kmers.bin_to_dna(2, 1, kmers.forward) == "G"
        check kmers.bin_to_dna(3, 1, kmers.forward) == "T"

        check kmers.bin_to_dna(0b10001111000100, 7, kmers.forward) == "GATTACA"
        check kmers.bin_to_dna(0b10001111000100, 7, kmers.reverse) == "TGTAATC"

    test "dna_to_kmers":
        check kmers.dna_to_kmers("AAAA", 2).seeds.len() == 6

    test "sorted_kmers":
        let
            sq = "ATCGGCTACTATT"
            expected = [
                "TCGGCTACTATT",
                "ATCGGCTACTAT",
                "ATCGGCTACTAT",
                "TCGGCTACTATT",
            ]
            k = 12
        var
            kms = kmers.dna_to_kmers(sq, k)
        check kms != nil
        let spot = kmers.initSpot(kms) # sort
        check kms == nil
        let got = sequtils.mapIt(spot.seeds, kmers.bin_to_dna(it.kmer,
                k.uint8,
                it.strand))
        check got == expected
        #check kmers.haskmer("AGCCGATGATAA", kms)

    test "search":
        let
            sq = "ATCGGCTACTATT"
            k = 12
            qms = kmers.dna_to_kmers(sq, k)
        var
            kms = kmers.dna_to_kmers(sq, k)
        let spot = kmers.initSpot(kms)
        let hits = kmers.search(spot, qms)
        check hits.len() == 4
        #check sets.toSet(seqUtils.toSeq(hits)).len() == 4 # 4 unique items
        check sets.len(sets.toHashSet(seqUtils.toSeq(deques.items(hits)))) ==
                4 # same as above

suite "kmers difference":
    let
        sq = "ATCGGCTACTATT"
        k = 12

    test "difference_of_self_is_nothing":
        let kms = kmers.dna_to_kmers(sq, k)
        #let qms = deepCopy(kms)
        var qms: kmers.pot_t
        deepCopy(qms, kms)
        check qms[] == kms[]
        check kmers.nkmers(qms) == 4
        check kmers.nkmers(kms) == 4
        let qspot = kmers.initSpot(qms)
        check qms == nil

        let kms0 = kmers.difference(kms, qspot)

        check kmers.nkmers(kms0) == 0
        check kmers.nkmers(kms) == 4
        check kmers.nkmers(qspot) == 4

        let
            expected: array[0, string] = []
            got = kmers.get_dnas(kms0)
        check got == expected

    test "difference_of_nothing_is_self":
        let kms = kmers.dna_to_kmers(sq, k)
        var qms = kmers.dna_to_kmers("", k)
        #let orig = deepCopy(kms)
        var orig: kmers.pot_t
        deepCopy(orig, kms)
        check kmers.nkmers(qms) == 0
        check kmers.nkmers(kms) == 4
        let qspot = kmers.initSpot(qms)
        check qms == nil

        let kms4 = kmers.difference(kms, qspot)

        check kmers.nkmers(kms4) == 4
        check kmers.nkmers(kms) == 4
        let got = kmers.get_dnas(kms4)
        let expected = kmers.get_dnas(orig)
        check got == expected

suite "kmer order":
    let dna = "ATGCGGACAGAAATATATACATAGAGACATACTCCCNAAAAAAAACTCAGAAGACACACATGCGCCC"
    let kms = kmers.dna_to_kmers(dna, 11)

    test "paired position of neg/pos strand":
        for i in countup(0, kms.seeds.len - 2, 2):
            check kms.seeds[i + 0].pos == kms.seeds[i + 1].pos

    test "increasing position":
        var lv: uint32 = 0
        for i in countup(2, kms.seeds.len - 2, 2):
            check (kms.seeds[i].pos > lv)
            lv = kms.seeds[i].pos

#proc test_FS*(pot: pot_t) =
#    for i in pot.seeds:
#        echo bin_to_dna(i.kmer,pot.word_size,i.strand)

suite "sparse_seeds":

    let dna = "CCCGAAAGTTT"
    let kms = kmers.dna_to_kmers(dna, 4)

    #test_FS(kms)

    test "spaced-seeds":
        let test_seeds = kmers.spacing_kmer(kms, 3)

        #echo dna.len, " ", test_seeds.seeds.len

        # 11 (dna len) - 2*4 (2 kmer lens) == 3
        # So with a space of 3, we have exactly 1 forward and 1 reverse seed.
        check(test_seeds.seeds.len == 2)

        block: # forward seed
            let x = 0
            let i = test_seeds.seeds[x]
            var myDNA = kmers.bin_to_dna(i.kmer, test_seeds.word_size, i.strand)
            #echo fmt"{myDNA} {dna} {test_seeds.word_size} {i.strand}"
            check((myDNA.len mod 2) == 0)
            check(myDNA == "CCCGGTTT")
            check(i.strand == kmers.forward)
        block: # reverse seed
            let x = 1
            let i = test_seeds.seeds[x]
            var myDNA = kmers.bin_to_dna(i.kmer, test_seeds.word_size, i.strand)
            #echo fmt"{myDNA} {dna} {test_seeds.word_size} {i.strand}"
            check((myDNA.len mod 2) == 0)
            check(myDNA == "CCCGGTTT")
            check(i.strand == kmers.reverse)
