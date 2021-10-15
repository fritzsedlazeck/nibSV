import nibsv
import unittest

suite "nibsv":
    test "foo":
        check parse_sv_allele("hi") == 1
