The implementation is does the following operations:

1. Generate reference and alternate k-mers for each break-point in the SV set.
2. Iterate over the fasta and limit k-mers from 1. to those that are unique across the genome and across svs
3. Iterate over the bam file and count occurrences of kmers from 2.

( After 3, each SV has 0 or more kmers that are unique to the alternate allele for the SV and 0 or more kmers
   that are unique to the reference allele that spans the SV breakpoints. Only sites with at least 1 unique
   reference kmer and alternate kmer are genotypable.)

4. Finally, iterate over each sv and report the max count (since there might be multiple kmers per sv).
   Note that this is not ideal since different samples could be genotyped at the same site with different kmers
   (though both will always have the info for the same kmers collected).

Note that *2.* is both critical and non-trivial as we must track which SV k-mers are seen
across the reference so we know they are not unique.
Much of the code involves switching back and forth between the kmers stored with their respective
Svs and the kmers stored globally to avoid `n^2` lookups.

