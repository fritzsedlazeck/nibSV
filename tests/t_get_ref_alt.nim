# vim: sts=2:ts=2:sw=2:et:tw=0
import nibsv
import unittest

suite "nibsv":
  test "foo":
    let
      s = "T[chr1_KI270709v1_random:461["
    check parse_sv_allele(s) == 1
