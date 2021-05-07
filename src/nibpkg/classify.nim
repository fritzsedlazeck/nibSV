import strutils
import tables
import hts
import ./read
import ./svidx
from ./compose import nil

proc buildSvIndex*(reference_path: string, vcf_path: string, flank: int, k: int, space: int): SvIndex =
  ## Open FASTA index
  var fai: Fai
  doAssert fai.open(reference_path), "Failed to open FASTA file: " & reference_path

  var variants: VCF
  doAssert(open(variants, vcf_path))

  result.kmerSize = k.uint8

  var sv_idx = 0
  echo "flank:", flank
  for v in variants:
    let sv_chrom = $v.CHROM

    let flanks = compose.retrieve_flanking_sequences_from_fai(fai, $v.CHROM, v.start.int, v.stop.int, flank)
    var p = compose.composePositioned(v, flanks.left, flanks.right, k, space)

    result.insert(p.sequences.alt_seq, k, sv_idx)
    # The insert function allows us to add to the ref count, but refmer also
    #  adds the same counts, so for now i'm commenting this out to minimize the
    #  lines of code we are debugging. --Zev
    # result.insert(p.sequences.ref_seq, k, -1)

    sv_idx.inc

proc classify_bam(filename: string, idx: SvIndex, k: int, spacedSeeds: bool, space: int, threads: int, fasta:cstring): CountTableRef[uint32] =
    new(result)

    var bamfile: Bam
    open(bamfile, filename, index = false, threads=threads, fai=fasta)
    var sequence: string
    var last_tid = -1

    for record in bamfile:
        if record.tid != last_tid:
          last_tid = record.tid
          stderr.write_line "on chrom:", record.chrom
        # NOTE: we may also want to filter record.flag.dup in the future, but
        # that will make results differ between bam and fastq
        if record.flag.secondary or record.flag.supplementary: continue
        record.sequence(sequence)

        var read_classification = process_read(sequence, idx, k, spacedSeeds, space)

        #if read_classification.compatible_SVs.len != 0:
        #    echo read_classification

        filter_read_matches(read_classification, winner_takes_all=false)
        for svId, count in read_classification.compatible_SVs:
            result.inc(svId)

    #echo result


proc classify_file*(filename: string, idx: SvIndex, k: int, spacedSeeds: bool, space: int, fasta:cstring): CountTableRef[uint32] =
    if endsWith(filename, ".bam") or endsWith(filename, ".cram"):
        return classify_bam(filename, idx, k, spacedSeeds, space, threads=2, fasta=fasta)
    else:
        quit("Error: only BAM input currently supported.")
