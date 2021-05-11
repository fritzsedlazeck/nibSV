import argparse
import sets
import tables
import strformat
import sequtils
import strutils
import algorithm
import hts
import kmer

#[
steps:
1. generate ref and alt kmers for given SV VCF
2. scan entire reference and remove any:
   a. alt kmers present
   b. ref kmers present more than 1x.
3. iterate over bam and count kmers.
]#


type Sv* = object
  chrom*: string
  pos*: int # 0-based position
  k: uint8
  ref_allele*: string
  alt_allele*: string
  ref_kmers*:seq[uint64]
  alt_kmers*:seq[uint64]

  # these correspond to the counts of kmers seen for each kmer in ref and alt
  # kmers respectively.
  ref_counts*: seq[uint32]
  alt_counts*: seq[uint32]

proc stop*(sv:Sv): int {.inline.} =
  result = sv.pos + sv.ref_allele.len

proc add_set(s:var seq[uint64], vals:HashSet[uint64]) =
  ## efficient add of a set to a seq
  doAssert s.len == 0
  s.setLen(vals.len)
  var i = 0
  for v in vals:
    s[i] = v
    i.inc

proc generate_kmers*(sv:var Sv, fai:Fai) =
  var ref_sequence = fai.get(sv.chrom, max(0, sv.pos - sv.k.int + 1), sv.stop + sv.k.int - 1)
  when defined(debug):
    doAssert sv.ref_allele in ref_sequence, $(sv.ref_allele, ref_sequence, sv.ref_allele.len, ref_sequence.len)
  # TODO AA: this might be .. < sv.k - 1]
  var alt_sequence = ref_sequence[0 ..< (sv.k - 1)]
  if sv.ref_allele.len == 1: # INS
    alt_sequence &= sv.alt_allele
  else:
    alt_sequence &= ref_sequence[(sv.k - 1)]

  # TODO AA: this might be .. ^sv.k with above change
  alt_sequence &= ref_sequence[^(sv.k.int) ..< ^0]
  when defined(debug):
    if sv.ref_allele.len == 1 or sv.alt_allele.len == 1:
      doAssert sv.alt_allele in alt_sequence, $(sv.alt_allele, alt_sequence, sv.alt_allele.len, alt_sequence.len, sv)

  # find unique ref and alt kmers.
  var refs = initHashSet[uint64]()
  var alts = initHashSet[uint64]()
  var refexcl = initHashSet[uint64]()
  for km in ref_sequence.slide(sv.k.int):
    refs.incl(km.enc)
  for km in alt_sequence.slide(sv.k.int):
    if km.enc in refs:
      # it's not unique to refs so we exclude from there.
      refexcl.incl(km.enc)
      continue
    alts.incl(km.enc)

  for km in refexcl:
    refs.excl(km)

  sv.ref_kmers.add_set(refs)
  sv.alt_kmers.add_set(alts)

proc get_all_kmers(svs:seq[Sv]): tuple[ref_kmers: HashSet[uint64], alt_kmers: HashSet[uint64], exclude: HashSet[uint64]] =
  # kmers are stored in each sv. we sometimes need them all together.
  # here we merge across svs into a single call-set wide Hash (for ref and alt)
  result.ref_kmers.init()
  result.alt_kmers.init()
  result.exclude.init()
  for sv in svs:
    for k in sv.ref_kmers:
      if k in result.alt_kmers or k in result.ref_kmers:
        result.exclude.incl(k)
      else:
        result.ref_kmers.incl(k)

    for k in sv.alt_kmers:
      if k in result.alt_kmers or k in result.ref_kmers:
        result.exclude.incl(k)
      else:
        result.alt_kmers.incl(k)

proc get_exclude(fai:Fai, all_ref_kmers: var HashSet[uint64], all_alt_kmers: var HashSet[uint64], exclude: var HashSet[uint64], k:int): tuple[refs_to_exclude:HashSet[uint64], alts_to_exclude: HashSet[uint64]] =
  # we want to exclude any kmers that are shared between ref and alt
  result.alts_to_exclude = exclude.union(all_ref_kmers.intersection(all_alt_kmers))
  result.refs_to_exclude = result.alts_to_exclude # value semantics so this makes a copy
  # now we want to iterate over the entire reference genome and exclude any
  # alt kmers that are in the reference and any ref kmers that are in more than
  # once.
  var ref_counts = newTable[uint64, int]()
  for k in all_ref_kmers: ref_counts[k] = 0

  stderr.write "[nibsv] "
  let chunk_size = 20_000_000
  for i in 0..<fai.len:
    var chrom = fai[i]
    var chrom_len = fai.chrom_len(chrom)
    if chrom_len > 10_000_000:
      stderr.write chrom, " "
      flushFile(stderr)
    for start in countup(0, chrom_len, 20_000_000):
      var start = max(0, start - k + 1) # redo to account for edge effects
      var sequence = fai.get(chrom, start, start + chunk_size)
      for kmer in sequence.slide(k):
        if kmer.enc in all_alt_kmers:
          result.alts_to_exclude.incl(kmer.enc)
        # we count every kmer that was one of our possible reference kmers and
        # we check that it's only seen once below
        if kmer.enc in ref_counts:
          ref_counts[kmer.enc].inc
  # now exclude any ref with a count > 1 as not unique.
  for k, cnt in ref_counts:
    if cnt > 1: result.refs_to_exclude.incl(k)
  stderr.write_line ""

proc remove(kmers:var seq[uint64], excludes:var HashSet[uint64]) =
  var excluded = 0
  for i, k in kmers:
    if k in excludes:
      kmers[i] = 0'u64
      excluded += 1
  if excluded != 0:
    var keep = newSeqOfCap[uint64](kmers.len - excluded)
    for v in kmers:
      if v == 0'u64: continue
      keep.add(v)
    kmers = keep

  # NOTE! we randomly remove about half of the kmers here
  # if we have more than 15. mostly these will be redundant anyway.
  #[
  if kmers.len >= 15:
    kmers.sort()
    var keep = newSeqOfCap[uint64](int(kmers.len/2) + 1)
    for i, k in kmers:
      if i mod 2 == 0: continue
      keep.add(k)
    kmers = keep
   ]#

proc remove_reference_kmers(svs:var seq[Sv], refs_to_exclude:var HashSet[uint64], alts_to_exclude:var HashSet[uint64]) =
  # now go back through svs and update each to remove any kmers that were 
  # present in reference or presnt > 1 time for refs
  for sv in svs.mitems:
    sv.ref_kmers.remove(refs_to_exclude)
    sv.alt_kmers.remove(alts_to_exclude)

proc to_kmer_cnt_table(svs: seq[Sv]): TableRef[uint64, int] =
  result = newTable[uint64, int]()
  for sv in svs:
    for k in sv.ref_kmers:
      when defined(check_unique_kmers):
        doAssert k notin result
      result[k] = 0
    for k in sv.alt_kmers:
      when defined(check_unique_kmers):
        doAssert k notin result
      result[k] = 0

proc count(svs:var seq[Sv], bam:Bam) =

  # kmer => count
  var kmer_cnts = svs.to_kmer_cnt_table()
  var sequence: string
  var k = svs[0].k.int
  stderr.write "[nibsv] "
  for tgt in bam.hdr.targets:
    stderr.write tgt.name, " "
    flushFile(stderr)
    for aln in bam.query(tgt.name):
      if aln.flag.supplementary or aln.flag.secondary: continue
      aln.sequence(sequence)

      for km in sequence.slide(k):
        # we only care about kmers that are in our sv set.
        if km.enc in kmer_cnts:
          kmer_cnts[km.enc] += 1
  stderr.write_line ""

  # transfer counts from single, global table back to per-variant counts.
  for sv in svs.mitems:
    for k in sv.ref_kmers:
      sv.ref_counts.add(kmer_cnts[k].uint32)
    for k in sv.alt_kmers:
      sv.alt_counts.add(kmer_cnts[k].uint32)

proc check_unique_kmers(svs: seq[Sv]) =
  # this is just an internal debugging function that checks assumptions
  stderr.write_line "[nibsv] checking that kmers are unique"
  var seen: HashSet[uint64]
  seen.init()
  for s in svs:
    for k in s.ref_kmers:
      #doAssert k notin seen
      seen.incl(k)
    for k in s.alt_kmers:
      doAssert k notin seen
      seen.incl(k)

type Stat = object
  unique_ref_only: int
  unique_alt_only: int
  unique_ref_and_alt: int
  no_unique: int

import hts/private/hts_concat

proc write(svs: seq[Sv], ivcf:VCF, output_path:string, sample_name:string) =
  ## write an output vcf with the counts.
  var ovcf:VCF
  if not ovcf.open(output_path, mode="w"):
    quit &"couldn't open output vcf: {output_path}"
  ovcf.copy_header(ivcf.header)
  discard ovcf.header.hdr.bcf_hdr_set_samples(nil, 0)
  discard ovcf.header.add_format("NIR", "1", "Integer", "nibsv: max reference counts. a value of -1 means that there were no suitable ref kmers for this sv")
  discard ovcf.header.add_format("NIA", "1", "Integer", "nibsv: max alternate counts. a value of -1 means that there were no suitable alt kmers for this sv")
  ovcf.add_sample(sample_name)
  doAssert ovcf.write_header()

  var i = 0
  for variant in ivcf:
    variant.vcf = ovcf
    let sv = svs[i]
    var max_ref = if sv.ref_counts.len > 0: @[sv.ref_counts.max.int32] else: @[-1'i32]
    doAssert variant.format.set("NIR", max_ref) == Status.OK
    var max_alt = if sv.alt_counts.len > 0: @[sv.alt_counts.max.int32] else: @[-1'i32]
    doAssert variant.format.set("NIA", max_alt) == Status.OK
    doAssert ovcf.write_variant(variant)
    i.inc
  doAssert i == svs.len
  ovcf.close()

proc stats(svs:seq[Sv]): Stat =
  for sv in svs:
    if sv.ref_kmers.len > 0:
      if sv.alt_kmers.len > 0:
        result.unique_ref_and_alt += 1
      else:
        result.unique_ref_only += 1
    else:
      if sv.alt_kmers.len > 0:
        result.unique_alt_only += 1
      else:
        result.no_unique += 1

proc sample_name(ibam:Bam): string =
  for l in ($ibam.hdr).split('\n'):
      if not l.startswith("@RG"): continue
      for t in l.split('\t'):
        if t.startswith("SM:"):
          return t.split(':')[1]

proc main() =

  var p = newParser("nibsv"):
    option("-k", default="31", help="kmer-size")
    option("-o", default="nibsv.vcf.gz", help="output vcf")
    arg("vcf", help="SV vcf with sites to genotype")
    arg("bam", help="bam or cram file for sample")
    arg("ref", help="reference fasta file")

  try:
    var a = p.parse()
  except UsageError as e:
    stderr.write_line(p.help)
    stderr.write_line(getCurrentExceptionMsg())
    quit(1)
  var a = p.parse()

  if a.help:
    quit 0
  var output_vcf = a.o
  let k = parseInt(a.k)
  if k > 31:
    quit "-k must be < 32"
  var ibam:Bam
  if not ibam.open(a.bam, threads=2, fai=a.ref, index=true):
    quit &"[nibsv] couldn't open cram file:{a.bam}"
  # options to do less work on cram decoding. we only need sequence and flag.
  var opts = SamField.SAM_FLAG.int or SamField.SAM_RNAME.int or SamField.SAM_POS.int or SamField.SAM_SEQ.int
  discard ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, opts)
  discard ibam.set_option(FormatOption.CRAM_OPT_DECODE_MD, 0)

  var ivcf:VCF
  if not ivcf.open(a.vcf, threads=2):
    quit &"[nibsv] couldn't open vcf file:{a.vcf}"

  var fai:Fai
  if not fai.open(a.ref):
    quit &"[nibsv] couldn't open fai file:{a.ref}"

  stderr.write_line "[nibsv] generating per-sv kmers"
  var svs: seq[Sv]
  for v in ivcf:
    #if v.REF.len != 1: continue
    var sv = Sv(ref_allele:v.REF, alt_allele:v.ALT[0], pos: v.start.int, chrom: $v.CHROM, k:k.uint8)

    # we have to skip bad variants. but leave them in so we keep the order.
    if v.ALT[0][0] == v.REF[0]:
      sv.generate_kmers(fai)
    svs.add(sv)

  ivcf.close()

  # get all kmers we've identified across break-points
  stderr.write_line "[nibsv] merging per-sv kmers to single set"
  var (all_ref_kmers, all_alt_kmers, exclude_kmers) = svs.get_all_kmers()

  stderr.write_line "[nibsv] finding kmers in reference in initial Sv set"
  var (refs_to_exclude, alts_to_exclude) = fai.get_exclude(all_ref_kmers, all_alt_kmers, exclude_kmers, k)

  all_ref_kmers.clear()
  all_alt_kmers.clear()
  exclude_kmers.clear()
  GC_fullCollect()

  stderr.write_line "[nibsv] removing reference kmers from sv sets"
  svs.remove_reference_kmers(refs_to_exclude, alts_to_exclude)

  refs_to_exclude.clear()
  alts_to_exclude.clear()
  GC_fullCollect()


  # now, we should only see each kmer once.
  when defined(check_unique_kmers):
    svs.check_unique_kmers()

  var stats = svs.stats()
  stderr.write_line &"[nibsv] stats: {stats} ({stats.unique_ref_and_alt} are genotype-able)"

  stderr.write_line "[nibsv] counting kmers in bam"
  svs.count(ibam)

  #for sv in svs:
  #  echo sv

  # re-open so we can write a new file
  if not ivcf.open(a.vcf, threads=2):
    quit &"[nibsv] couldn't open vcf file:{a.vcf}"

  svs.write(ivcf, output_vcf, ibam.sample_name)
  stderr.write_line &"[nibsv] wrote: {svs.len} variants to {output_vcf}"
  ivcf.close()

when isMainModule:
  main()

