# seqerror

Estimate **sequence-context-specific sequencing error rates** from a
`samtools mpileup` stream (2013).

## Why

Sequencing errors are not uniform. Some sequence contexts, such as particular
trinucleotides or bases after a homopolymer, produce many more errors than
others, and the error profile differs between the forward and reverse strands.
Variant and de novo mutation callers that assume a single global error rate
call false positives at these error-prone contexts. Measuring the error rate
per context gives a better background model.

## What it does

`seqError` reads mpileup output on standard input and, for every position with
read depth above 10:

1. Builds the reference **triplet** (the base and the two bases before it) for
   forward-strand reads, and the reverse-complement triplet for reverse-strand
   reads.
2. Counts reads that disagree with the reference on each strand.
3. Accumulates mismatch counts and per-site mismatch fractions per triplet.

At the end it prints one row per triplet seen more than 10 times.

## Build

```bash
make        # g++ seqError.cc
```

## Usage

```bash
samtools mpileup -f reference.fa sample.bam | ./a.out > triplet_error_rates.tsv
```

## Output

| Column | Meaning |
| --- | --- |
| `Triplet` | Reference trinucleotide context (strand-aware) |
| `Triplet_cumul_mismatch_count` | Total mismatching reads at this context |
| `Triplet_cumul_mismatch_fraction` | Sum of per-site mismatch fractions |
| `Triplet_count` | Number of sites with this context |
| `Normalized_mismatch_fraction` | Mean mismatch fraction per site, the estimated error rate |

The total number of indels seen is reported on standard error.

## Related

This idea of site- and context-specific error modelling later grew into
[Arnav](https://github.com/gatoravi/arnav), which learns per-site error rates
from a panel of samples to call variants from RNA-seq.

## License

GPL-3.0 (see the header of `seqError.cc`).
