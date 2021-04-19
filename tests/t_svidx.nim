# vim: sw=1 ts=1 sts=1 tw=0 et:
import unittest
import nibpkg/svidx
import tables

suite "SvIndex suite":
 test "that sv insertion works":

  var idx: SvIndex
  idx.insert("ATCGGCTACTATT", 11, 2)

  for kmer, t in idx.counts:
   check t.svs == @[2'u32]

 test "that no ref insertion occurs unless kmer matches":

  var idx: SvIndex
  idx.insert("ATCGGCTACTATT", 11, -1)
  check idx.len == 0

  idx.insert("ATCGGCTACTATT", 11, 2)
  idx.insert("ATCGGCTACTATT", 11, -1)

  for kmer, t in idx.counts:
   check t.svs == @[2'u32]
   check t.refCount == 1'u32

 test "that filter removes SV entries with refcount gt zero":
  var idx: SvIndex

  idx.insert("ATCGGCTACTATT", 11, 2)
  idx.insert("ATCGGCTACTATT", 11, -1)

  filterRefKmers(idx, 0)
  check idx.len == 0;

test "that filter does not SV entries with refcount <= two":
 var idx: SvIndex

 idx.insert("ATCGGCTACTATT", 11, 2)
 idx.insert("ATCGGCTACTATT", 11, -1)
 idx.insert("ATCGGCTACTATT", 11, -1)

 filterRefKmers(idx, 2)
 check idx.len == 6;

test "that filter does not SV entries with refcount <= 0":
 var idx: SvIndex

 idx.insert("ATCGGCTACTATT", 11, 2)


 filterRefKmers(idx, 2)
 check idx.len == 6;
