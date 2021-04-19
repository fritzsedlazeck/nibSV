import unittest
import nibpkg/compose
import hts

suite "compose suite":
  # TODO this test only partially covers the function. We need to check the full
  #  sequences. For now some tests are better than none.

  test "check that the variant haplotypes are correctly constructed":
    var variants: VCF
    doAssert(open(variants, "../test-data/GIAB_PBSV_TRIO_CALLS_TEST2.vcf"))


    for v in variants:
      var variant_type: string
      doAssert v.info.get("SVTYPE", variant_type) == Status.OK
      var variant_seq = composePositioned(v, "AAAAA", "TTTTT")
      if variant_type == "DEL":
        check(variant_seq.sequences.alt_seq == "AAAAATTTTT")
      if variant_type == "INS":
        check(variant_seq.sequences.alt_seq != "AAAAATTTTT")

  test "TODO check positions are correct!":
    echo "Please fill me in."
  test "TODO check full alleles are ok!":
    echo "Please fill me in."
