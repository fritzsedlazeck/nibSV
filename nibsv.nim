# vim: sts=2:ts=2:sw=2:et:tw=0
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


const nibsvVersion* = "0.0.2"
const nibsvGitCommit* = staticExec("git rev-parse --verify HEAD")

type Sv* = object
  chrom*: string
  pos*: int # 0-based position
  k: uint8
  space: seq[uint8]
  ref_allele*: string
  alt_allele*: string
  ref_kmers*:seq[uint64]
  alt_kmers*:seq[uint64]

  # these correspond to the counts of kmers seen for each kmer in ref and alt
  # kmers respectively.
  ref_counts*: seq[uint32]
  alt_counts*: seq[uint32]


proc add_set(s:var seq[uint64], vals:HashSet[uint64]) =
  ## efficient add of a set to a seq
  doAssert s.len == 0
  var i = s.len
  s.setLen(s.len + vals.len)
  for v in vals:
    s[i] = v
    i.inc

proc kmer_size(sv:Sv): seq[int] =
  ## size of sequence for kmer, including space
  for s in sv.space:
    result.add(sv.k.int + s.int + int(s > 0) * sv.k.int)
  if result.len == 0: result.add(sv.k.int)

proc update_kmers(sv:var Sv, ref_sequences:seq[string], alt_sequences:seq[string], step:int) =
  # find unique ref and alt kmers.
  var refs = initHashSet[uint64]()
  var alts = initHashSet[uint64]()
  var refexcl = initHashSet[uint64]()

  for ref_sequence in ref_sequences:
    let skip_start = if ref_sequence.len > 100: max(sv.kmer_size).int + 10 else: ref_sequence.len
    let skip_stop = if ref_sequence.len > 100: ref_sequence.len - max(sv.kmer_size).int - 10 else: ref_sequence.len
    var ki = -1
    for space in sv.space:
      let space = space.uint64
      for km in ref_sequence.slide_space(sv.k.int, space):
        ki.inc
        if ki >= skip_start and ki <= skip_stop: continue
        refs.incl(km.enc)

  for alt_sequence in alt_sequences:
    # for large insertion sequences we don't need all of the internal kmers so here, we skip the internal kmers.
    let skip_start = if alt_sequence.len > 100: max(sv.kmer_size).int + 10 else: alt_sequence.len
    let skip_stop = if alt_sequence.len > 100: alt_sequence.len - max(sv.kmer_size).int - 10 else: alt_sequence.len
    for space in sv.space:
      var ki = -1
      for km in alt_sequence.slide_space(sv.k.int, space):
        ki.inc
        if ki >= skip_start and ki <= skip_stop: continue
        if km.enc in refs:
          # it's not unique to refs so we exclude from there.
          refexcl.incl(km.enc)
          continue
        alts.incl(km.enc)

  for km in refexcl:
    refs.excl(km)

  sv.ref_kmers.add_set(refs)
  sv.alt_kmers.add_set(alts)

proc stop*(sv:Sv): int {.inline.} =
  result = sv.pos + sv.ref_allele.len

proc parse_sv_allele*(sv_allele: string): int =
  var first_parens_index: int = 0
  var second_parens_index: int = 0
  var pre_bases: string
  var post_bases: string
  var chrom: string
  var pos: int

  while sv_allele[first_parens_index] != '[' and sv_allele[first_parens_index] != ']':
    first_parens_index.inc

  while sv_allele[second_parens_index] != '[' and sv_allele[second_parens_index] != ']':
    second_parens_index.inc
  
  if first_parens_index > 0:
    pre_bases = sv_allele[0 ..< first_parens_index]
  if second_parens_index > 0:
    post_bases = sv_allele[second_parens_index ..< sv_allele.len]
  
  var chrom_pos_splits: seq[string] = (sv_allele[first_parens_index ..< second_parens_index]).split(sep=':')
  chrom = chrom_pos_splits[0]
  pos = strUtils.parseInt(chrom_pos_splits[1])
  
  

proc generate_ref_alt*(sv:var Sv, fai:Fai, overlap:uint8=6): tuple[ref_sequence:seq[string], alt_sequence:seq[string]] =
  let overlap = overlap.int

  for i, kmer_size in sv.kmer_size:
    result.ref_sequence.add( fai.get(sv.chrom, max(0, sv.pos - kmer_size + overlap), sv.stop + kmer_size - overlap))
    # reference goes back from start:
    #     k + space + k
    when defined(debug):
      doAssert sv.ref_allele in result.ref_sequence[i], $(sv.ref_allele, result.ref_sequence, sv.ref_allele.len, result.ref_sequence.len)

    result.alt_sequence.add(result.ref_sequence[i][0 ..< (kmer_size - overlap)])
    doAssert sv.ref_allele[0] == result.ref_sequence[i][(kmer_size - overlap)]
    if sv.ref_allele.len == 1: # INS
      #if sv.alt_allele.len < 100:
      result.alt_sequence[i] &= sv.alt_allele
      #else:
      #  result.alt_sequence[i] &= sv.alt_allele[0..(kmer_size + 10)] & repeat('N', kmer_size) & sv.alt_allele[^(sv.alt_allele.len - kmer_size - 10)..<sv.alt_allele.len]
    else:
    # NOTE: this doesn't work for REF and ALT lengths > 1
      result.alt_sequence[i] &= result.ref_sequence[i][kmer_size - overlap]

    ## Example: kmer size = 5 and overlap is 3, ref_seq is of len 10:
    ## 10 - (5 - 3 + 1) -> 10 - (3 - 1)
    ## TODO: possible double substraction at ^(overlap - 1)
    ## Number of elements - (value)
    ## e.g. ^1 == n_elements -1 == last_value_in_seq
    result.alt_sequence[i] &= result.ref_sequence[i][^(kmer_size.int - overlap + 1) ..< ^(overlap - 1)]
    when defined(debug):
      if sv.ref_allele.len == 1 or sv.alt_allele.len == 1:
        doAssert sv.alt_allele in result.alt_sequence[i], $(sv.alt_allele, result.alt_sequence, sv.alt_allele.len, result.alt_sequence.len, sv)


proc generate_kmers*(sv:var Sv, fai:Fai, step:int) =
  let (ref_sequences, alt_sequences) = sv.generate_ref_alt(fai)
  sv.update_kmers(ref_sequences, alt_sequences, step)

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

type Excluder = ref object
  refs_to_exclude:HashSet[uint64]
  alts_to_exclude:HashSet[uint64]

proc get_exclude(fai:Fai, all_ref_kmers: var HashSet[uint64], all_alt_kmers: var HashSet[uint64], exclude: var HashSet[uint64], k:int, spaces:seq[uint8]): Excluder =
  result = new(Excluder)
  # we want to exclude any kmers that are shared between ref and alt
  result.alts_to_exclude = exclude.union(all_ref_kmers.intersection(all_alt_kmers))
  exclude.clear()
  result.refs_to_exclude = result.alts_to_exclude # value semantics so this makes a copy
  # now we want to iterate over the entire reference genome and exclude any
  # alt kmers that are in the reference and any ref kmers that are in more than
  # once.
  var ref_counts = newTable[uint64, uint8](all_ref_kmers.len)
  for km in all_ref_kmers: ref_counts[km] = 0'u8

  stderr.write "[nibsv] "
  let chunk_size = 20_000_000
  for i in 0..<fai.len:
    var chrom = fai[i]
    var chrom_len = fai.chrom_len(chrom)
    if chrom_len > 10_000_000:
      stderr.write chrom, " "
      flushFile(stderr)
    for ostart in countup(0, chrom_len, 20_000_000):
      var start = max(0, ostart - k + 1) # redo to account for edge effects
      var sequence = fai.get(chrom, start, ostart + chunk_size)
      for space in spaces:
        for kmer in sequence.slide_space(k, space.uint64):
          if kmer.enc in all_alt_kmers:
            result.alts_to_exclude.incl(kmer.enc)
          # we count every kmer that was one of our possible reference kmers and
          # we check that it's only seen once below
          # if it's 2, we already know we can exclude and if it's not in the
          # table (default of 2, we can also skip)
          if ref_counts.getOrDefault(kmer.enc, 2'u8) != 2'u8:
            ref_counts[kmer.enc].inc
  # now exclude any ref with a count > 1 as not unique.
  for k, cnt in ref_counts:
    if cnt == 2'u8: result.refs_to_exclude.incl(k)
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
  #if kmers.len >= 100:
  #  let imod = int(kmers.len / 50 + 1)
  #  kmers.sort()
  #  var keep = newSeqOfCap[uint64](int(kmers.len/2) + 1)
  #  for i, k in kmers:
  #    if i mod imod == 0: continue
  #    keep.add(k)
  #  kmers = keep

proc remove_reference_kmers(svs:var seq[Sv], ex:Excluder) =
  # now go back through svs and update each to remove any kmers that were 
  # present in reference or presnt > 1 time for refs
  for sv in svs.mitems:
    sv.ref_kmers.remove(ex.refs_to_exclude)
    sv.alt_kmers.remove(ex.alts_to_exclude)

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
  var spaces = newSeq[uint64](svs[0].space.len)
  for s in svs[0].space:
    spaces.add(s.uint64)

  stderr.write "[nibsv] "
  for tgt in bam.hdr.targets:
    stderr.write tgt.name, " "
    flushFile(stderr)
    for aln in bam.query(tgt.name):
      if aln.flag.supplementary or aln.flag.secondary: continue
      aln.sequence(sequence)

      for space in spaces:
        for km in sequence.slide_space(k, space):
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

proc argmax(a:seq[uint32], maxval:uint32): int =
  if len(a) == 0: return -1
  result = -1
  var m = -1
  if a[0] <= maxval:
      m = a[0].int
      result = 0
  for i, v in a:
    if v.int > m and v <= maxval:
      result = i
      m = v.int

type vali = object
  val: uint32
  i: uint32

proc vali_cmp(a, b:vali): int =
  if a.val < b.val: return -1
  elif a.val == b.val: return 0
  return 1

proc argmed(a:seq[uint32], maxval:uint32): int =
  if len(a) == 0: return -1

  var ai = newSeqOfCap[vali](a.len)
  for i, val in a:
    if val > maxval: continue
    ai.add(vali(i:i.uint32, val: val))
  if ai.len == 0:
    return -1
  ai.sort(vali_cmp)
  return ai[len(ai) /% 2].i.int

proc get_ith_kmer(kmers:seq[uint64], k:int, space:int, idx:int): string =
  if kmers.len == 0 or idx < 0: return ""
  var km = kmers[idx]
  result = newString(if space == 0: k else: 2 * k)
  km.decode(result)

proc write(svs: seq[Sv], ivcf:VCF, output_path:string, sample_name:string, maxval:uint32, use_med:bool) =
  ## write an output vcf with the counts.
  var ovcf:VCF
  if not ovcf.open(output_path, mode="w"):
    quit &"couldn't open output vcf: {output_path}"

  ovcf.copy_header(ivcf.header)
  discard ovcf.header.hdr.bcf_hdr_set_samples(nil, 0)
  discard ovcf.header.add_format("NIRK", "1", "String", "nibsv: reference kmer (this and reverse-complement are used)")
  discard ovcf.header.add_format("NIAK", "1", "String", "nibsv: alternate kmer (this and reverse-complement are used)")
  discard ovcf.header.add_format("NIR", "1", "Integer", "nibsv: max reference counts. a value of -1 means that there were no suitable ref kmers for this sv")
  discard ovcf.header.add_format("NIA", "1", "Integer", "nibsv: max alternate counts. a value of -1 means that there were no suitable alt kmers for this sv")
  ovcf.add_sample(sample_name)
  doAssert ovcf.write_header()

  var i = 0
  for variant in ivcf:
    variant.vcf = ovcf
    let sv = svs[i]
    let ref_max = if use_med: sv.ref_counts.argmed(maxval=maxval) else: sv.ref_counts.argmax(maxval=maxval)
    var kms = @[get_ith_kmer(sv.ref_kmers, sv.k.int, sv.space[0].int, ref_max)]
    var max_ref = if sv.ref_counts.len > 0 and ref_max >= 0: @[sv.ref_counts[ref_max].int32] else: @[-1'i32]
    doAssert variant.format.set("NIR", max_ref) == Status.OK
    if kms[0] != "":
      doAssert variant.format.set("NIRK", kms) == Status.OK

    let alt_max = if use_med: sv.alt_counts.argmed(maxval=maxval) else: sv.alt_counts.argmax(maxval=maxval)
    kms = @[get_ith_kmer(sv.alt_kmers, sv.k.int, sv.space[0].int, alt_max)]
    var max_alt = if sv.alt_counts.len > 0 and alt_max >= 0: @[sv.alt_counts[alt_max].int32] else: @[-1'i32]
    doAssert variant.format.set("NIA", max_alt) == Status.OK

    if kms[0] != "":
      doAssert variant.format.set("NIAK", kms) == Status.OK

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
  stderr.write_line "[nibsv] no read-groups found, using 'SAMPLE' as sample id"

proc main() =

  when not defined(danger):
    stderr.write_line "[nibsv] WARNING: nibsv compiled without optimizations; will be slow"

  stderr.write_line &"[nibsv] version: {nibsvVersion} commit: {nibsvGitCommit}"

  var p = newParser("nibsv"):
    option("-k", default="27", help="kmer-size must be <= 15 if space > 0 else 31")
    option("--space", default="", help="space between kmers", multiple=true)
    flag("--use-med", help="use median instead of max for choosing kmer count")
    option("--step", default="1", help="step between generated reference kmers (larger values save more memory)")
    option("-o", default="nibsv.vcf.gz", help="output vcf")
    option("--cram-ref", help="optional reference fasta file for cram if difference from reference fasta")
    arg("vcf", help="SV vcf with sites to genotype")
    arg("bam", help="bam or cram file for sample")
    arg("ref", help="reference fasta file")

  let maxval = 500'u32 # TODO: make this a parameter
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
  var spaces = newSeqOfCap[uint8](a.space.len)
  for sp in a.space:
    let space = parseInt(sp)
    if space > uint8.high.int:
      quit "--space must be < 256"
    if space > 0:
      if k >= 16:
        quit "-k must be < 16 when space is > 0"
    spaces.add(space.uint8)
  if spaces.len == 0: spaces.add(0'u8)
  var ibam:Bam
  let step = parseInt(a.step)
  if a.cram_ref == "":
    a.cram_ref = a.ref
  if not ibam.open(a.bam, threads=2, fai=a.cram_ref, index=true):
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
    #if v.FILTER notin ["PASS", "LongReadHomRef"]: continue
    #if v.REF.len != 1: continue
    var sv = Sv(ref_allele: $v.REF, alt_allele: $v.ALT[0], pos: v.start.int, chrom: $v.CHROM, k:k.uint8, space:spaces)

    # we have to skip bad variants. but leave them in so we keep the order.
    if v.ALT[0][0] == v.REF[0]: # and (v.ALT[0].len > ml or v.REF.len > ml):
      sv.generate_kmers(fai, step)
    svs.add(sv)

  ivcf.close()

  var bstats = svs.stats()
  stderr.write_line "[nibsv] before stats:", bstats
  # get all kmers we've identified across break-points
  stderr.write_line "[nibsv] merging per-sv kmers to single set"
  var (all_ref_kmers, all_alt_kmers, exclude_kmers) = svs.get_all_kmers()

  stderr.write_line "[nibsv] finding kmers in reference in initial Sv set:"
  var ex = fai.get_exclude(all_ref_kmers, all_alt_kmers, exclude_kmers, k, spaces)

  all_ref_kmers.clear()
  all_alt_kmers.clear()
  exclude_kmers.clear()
  GC_fullCollect()

  stderr.write_line "[nibsv] removing reference kmers from sv sets"
  svs.remove_reference_kmers(ex)

  ex.refs_to_exclude.clear()
  ex.alts_to_exclude.clear()
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

  if ivcf.header.add_string(&"##nibsv-info=\"version:{nibsvVersion} commit:{nibsvGitCommit} k:{a.k} space:{a.space} ref:{a.ref} cram_ref:{a.cram_ref}\"") != Status.OK:
    stderr.write_line "[nibsv] warning! couldn't add nibsv-info to the header of output vcf"

  svs.write(ivcf, output_vcf, ibam.sample_name, maxval, a.use_med)
  stderr.write_line &"[nibsv] wrote: {svs.len} variants to {output_vcf}"
  ivcf.close()

when isMainModule:
  main()

