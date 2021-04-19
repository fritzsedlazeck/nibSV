# vim: sw=4 ts=4 sts=4 tw=0 et:
import nibpkg/util
import unittest
from strformat import fmt

from os import nil
from sequtils import nil
from strutils import nil

suite "util":
    test "thousands":
        check thousands(0) == "0"
        check thousands(1) == "1"
        check thousands(10) == "10"
        check thousands(100) == "100"
        check thousands(1_000) == "1,000"
        check thousands(10_000) == "10,000"
        check thousands(100_000) == "100,000"
        check thousands(1_000_000) == "1,000,000"
        check thousands(-1_000_000) == "-1,000,000"
        check thousands(-10_000) == "-10,000"
        check thousands(-1_000) == "-1,000"
        check thousands(-1) == "-1"
        check thousands(-0) == "0"
    test "splitWeighted":
        check splitWeighted(0, @[]) == []
        check splitWeighted(0, @[42]) == []
        check splitWeighted(1, @[42]) == [1]
        check splitWeighted(2, @[42, 2]) == [1, 1]
        check splitWeighted(3, @[1, 1, 1]) == [1, 1, 1]
        check splitWeighted(3, @[1, 1, 1, 1]) == [2, 1, 1]
        check splitWeighted(3, @[1, 1]) == [1, 1]
        check splitWeighted(1, @[1, 2, 3, 4]) == [4]
        check splitWeighted(2, @[1, 2, 3, 4]) == [3, 1]
        check splitWeighted(3, @[1, 2, 3, 4]) == [3, 1] # greedy
        check splitWeighted(4, @[1, 2, 3, 4]) == [2, 1, 1] # greedy, so order matters
        check splitWeighted(4, @[4, 3, 2, 1]) == [1, 1, 1, 1] # see?
        check splitWeighted(3, @[4, 3, 2, 1]) == [1, 1, 2]
        check splitWeighted(2, @[4, 3, 2, 1]) == [2, 2]
        check splitWeighted(1, @[4, 3, 2, 1]) == [4]
        check splitWeighted(4, @[4191650, 4009608, 4154778, 4096102]) == [1, 2, 1] # not very good
    test "partitionWeighted":
        check partitionWeighted(4, @[4191650, 4009608, 4154778, 4096102]) == @[@[1], @[3], @[2], @[0]]
        check partitionWeighted(4, @[4191650, 4009608, 4154778, 4096102, 99]) == @[@[1, 4], @[3], @[2], @[0]]
        check partitionWeighted(4, @[4_191_650, 4_009_608, 4_154_778, 4_096_102, 500_000]) == @[@[3], @[2], @[0], @[1, 4]]
        check partitionWeighted(2, @[1, 3, 5, 2, 4, 6]) == @[@[2, 4, 0], @[5, 1, 3]]
        check partitionWeighted(3, @[1, 2, 2]) == @[@[0], @[2], @[1]]
    test "combineToTarget":
        proc icombineToTarget(t: int, weights: seq[int]): seq[seq[int]] =
            return combineToTarget(t, sequtils.mapIt(weights, int64(it)))
        check icombineToTarget(3, @[2, 2, 2, 2]) == @[@[0, 1], @[2, 3]]
        check icombineToTarget(3, @[2, 2, 2]) == @[@[0, 1], @[2]]
        check icombineToTarget(3, @[2, 2]) == @[@[0, 1]]
        check icombineToTarget(2, @[2, 2, 2, 2]) == @[@[0], @[1], @[2], @[3]]
        check icombineToTarget(1, @[2, 2, 2, 2]) == @[@[0], @[1], @[2], @[3]]
        check icombineToTarget(4, @[2, 2, 2, 2]) == @[@[0, 1], @[2, 3]]
        check icombineToTarget(3, @[1, 2, 3, 4]) == @[@[0, 1], @[2], @[3]]
        check icombineToTarget(3, @[1, 1, 2, 1]) == @[@[0, 1, 2], @[3]]
        check icombineToTarget(3, @[1, 2, 1, 1]) == @[@[0, 1], @[2, 3]]

    test "sscanf":
        let s_frmt = strutils.format("%ld %$#[^\n]",
            (util.MAX_HEADROOM - 1))
        var
            bufAname: util.Headroom
            name: string
            val: int32
            line: string
        line = "123 abc def"
        let scanned = util.sscanf(line.cstring, s_frmt.cstring,
            addr val, addr bufAname)
        check val == 123
        check scanned == 2
        util.toString(bufAname, name, line)
        check "abc def" == name

    test "isEmptyFile":
        let fn = "empty.txt"
        check os.execShellCmd("rm -f {fn}".fmt) == 0
        check os.execShellCmd("touch {fn}".fmt) == 0
        check isEmptyFile(fn)
        check os.execShellCmd("echo fuller >> {fn}".fmt) == 0
        check not isEmptyFile(fn)
        check os.execShellCmd("rm -f {fn}".fmt) == 0

    test "isOlderFile":
        let afn = "a.txt"
        let bfn = "b.txt"
        check os.execShellCmd("touch {afn}".fmt) == 0
        check os.execShellCmd("touch {bfn}".fmt) == 0
        # We cannot reliably test for strictly older because of fsys probs.
        check not isOlderFile(bfn, afn)
        check os.execShellCmd("touch {afn}".fmt) == 0
        check not isOlderFile(afn, bfn)
        check os.execShellCmd("rm -f {afn} {bfn}".fmt) == 0

    test "getNthWord":
        check getNthWord("a", 0, ' ') == "a"
        check getNthWord("a b", 0, ' ') == "a"
        check getNthWord("a b", 1, ' ') == "b"
        check getNthWord("a b ", 1, ' ') == "b"
        check getNthWord("ax bx cx", 0, ' ') == "ax"
        check getNthWord("ax bx cx", 1, ' ') == "bx"
        check getNthWord("ax bx cx", 2, ' ') == "cx"
        expect PbError:
            discard getNthWord("ax bx cx", 3, ' ')
