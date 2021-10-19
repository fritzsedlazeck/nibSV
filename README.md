<img align="right" width="300" height="300" src="https://github.com/collaborativebioinformatics/nibSV/blob/main/SVNibbler.png">

# NibblerSV

Fading out to https://github.com/brentp/nibsv <-- 
## Intro 

**nibsv** genotypes base-resolved structural variants by counting the occurence of
novel k-mers generated at the break-points.

This means that `nibsv` has very **high specificity**, but because there are
some SVs that do not generate novel k-mers, it can have a **lower
sensitivity**.

# Usage

To run nibsv, download the static binary from the releases and run as:
```
nibsv -o ${genotyped_vcf} \
  -k 13 --space 17 --space 18 \
  ${sites_vcf}  ${bam_or_cram} ${reference_fasta}
```

Full usage:
```
[nibsv] version: 0.0.2 commit: b1b9e06937d519fed27f9032b26a285367106354
nibsv

Usage:
  nibsv [options] vcf bam ref

Arguments:
  vcf              SV vcf with sites to genotype
  bam              bam or cram file for sample
  ref              reference fasta file

Options:
  -k=K                       kmer-size must be <= 15 if space > 0 else 31 (default: 27)
  --space=SPACE              space between kmers (can be specified multiple times)
  -o=O                       output vcf (default: nibsv.vcf.gz)
  --cram-ref=CRAM_REF        optional reference fasta file for cram if difference from reference fasta
  -h, --help                 Show this help

```

## Contributors

Brent Pedersen<sup>1</sup>, Christopher Dunn<sup>2</sup>, Eric Dawson<sup>3</sup>, Fritz Sedlazeck<sup>4</sup>, Peter Xie<sup>5</sup>, and Zev Kronenberg<sup>2</sup>

<sup>1</sup> University of Utrecht; <sup>2</sup>PacBio; <sup>3</sup>Nvidia Corporation; <sup>4</sup>Baylor College of Medicine; <sup>5</sup>JBrowse (UC Berkeley);
