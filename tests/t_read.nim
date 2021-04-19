import unittest
import nibpkg/read

suite "read suite":
  test "test that filter read works":

    var r = Read(compatible_SVs: initCountTable[uint32]())
    r.compatible_SVs.inc(23, 5)
    r.compatible_SVs.inc(22, 1)

    r.filter_read_matches()

    check r.compatible_SVs.len == 1
    check 23 in r.compatible_SVs

