# `vcfdo`: do stuff with a VCF file

## Installation
To use `vcfdo` requires a working **Python 3** installation as well as the following packages:

* [`numpy`](http://www.numpy.org/)
* [`cyvcf2`](https://github.com/brentp/cyvcf2)
* [`pyfaidx`](https://pypi.org/project/pyfaidx/)
* [`scikit-allel`](https://github.com/cggh/scikit-allel)

Assuming the user either can write to the global `site-packages` directory (ie. has root access or a `conda`-type environment), `vcfdo` can be installed straight from Github with `pip`:
```
pip install git+https://github.com/IDEELResearch/vcfdo
```
Dependencies will be automatically installed or updated as required. The main executable, also called `vcfdo`, will be made available on the user's `$PATH`.

## Use patterns
Utilities in `vcfdo` are compartmentalized into sub-commands in the style of `samtools` or `bcftools`. To see a list, run `vcfdo --help` to see a message like:

```
usage:
	vcfdo <command> [options]

	Commands:

	  --- Basic statistics
		  missing           rates of missing calls by site or by sample
		  depthsummary      quantiles of read depth (at called sites) per sample
		  filtersummary     tallies of filter status by site or by sample

	  --- Downsampling
		  thin              emit every nth site, possibly with random starting point
		  prune             greedy LD-pruning in sliding windows

	  --- Ancestral alleles
		  polarize          define ancestral alleles based on outgroup sequence
		  titv              annotate sites as transitions or transversions; flag gBGC candidates
		  derived           count derived alleles per sample, using either read counts or hard calls

	  --- Allele frequencies
		  wsaf              calculate within-sample allele frequency (WSAF) and related quantities from read counts
		  fws               estimate pseudo-inbreeding coefficient F_ws ("within-sample relatedness")
		  private           annotate sites where non-ref allele is only found in certain subset of samples
		  sfs               approximate the n-dimensional unfolded SFS by sampling

	  --- Relatedness/ordination
		  dist              calculate LD-weighted pairwise distances from within-sample allele frequencies
		  ibs               calculate simple identity-by-state (IBS) matrix from hard calls at polymorphic sites only
		  pca               perform PCA using within-sample allele frequencies instead of hard calls

```

Some options are common to all commands; these are
```
-h, --help            show this help message and exit
-i INFILE, --infile INFILE
					  VCF file, possibly gzipped [default: stdin]
-n NSITES, --nsites NSITES
					  stop after processing this many sites [default: no
					  limit]
--chunk-size CHUNK_SIZE
					  when loading genotypes into memory, process this many
					  sites per chunk [default: 2000]
--seed SEED           random number seed; only affects stuff involving
					  sampling [default: 1]
-q, --quiet           suppress logging messages
-v, --verbose		  see more logging messages
```

All of the utilities in `vcfdo` work on streams, and are intended to be used that way. The goal is to store on disk only those attributes of a VCF that are invariant under subsetting by rows (sites) or columns (samples). Everything else can be calculated on the fly as needed. Simple filtering operations such as those available in `bcftools view` are not duplicated in `vcfdo` -- wherever possible, it will be more efficient to do what you can with `bcftools view` and stream the result to `vcfdo`.

On that note, the VCF parser used in `vcfdo` is `cyvcf2`, which itself is a thin Python wrapper around `htslib`, so any VCF that is compatible with `bcftools` should work with `vcfdo`.

For example, suppose we have a (gzipped) VCF `file.vcf.gz` and a list of samples of interest in `samples.txt`. To perform PCA on the within-sample allele frequencies (WSAF), keeping only sites with population-level minor allele frequency (PLMAF) > 0.01, the command would be:

```
bcftools view --samples-file samples.txt file.vcf.gz | \
vcfdo wsaf | \
bcftools view -i 'INFO/PLMAF > 0.01' | \
vcfdo pca -o my_pca
```
This produces a PCA result in a pair of files `my_pca.pca` (sample projections) and `my_pca.pca.ev` (eigenvalues).

## General remarks

* All utilities can be limited to the first handful of sites in the file by passing `-n {number of sites to visit}`. It is recommended to test-drive commands on a small number of sites before unleashing them in a multi-gigabyte file.
* All utilities emit rather verbose log messages, including regular progress updates while streaming through the input VCF file. To suppress this, pass `--quiet`.
* All utilities that produce VCF output leave a "breadcrumb" in the VCF header documenting the command-line call, current working directory and timestamp.
* All utilities that produce VCF output will dump uncompressed VCF, starting with a valid header, to `stdout`.
* Default input for all utilities is `stdin`. An unfortunate side-effect that I haven't been able to fix is that calling `vcfdo <command>` with no arguments and no piped input will cause the command-line to hang until the user interrupts it with `Ctrl-D` (not `Ctrl-C`; why??)
* When a list of samples is required at the command line, `vcfdo` will check it against the VCF header and remove duplicates. Samples requsted but not present in the VCF are ignored with a warning.

## Definitions and conventions
As the target audience for this tool is people working on _Plasmodium_ spp, heavy use is made of some quantities that are bespoke to malariologists. These include:

* Genotype-level quantities (these live in the `FORMAT` field)
	* **WSAF** (within-sample allele frequency) -- proportion of reads at this site supporting _all_ non-reference alleles -- float [0,1]
	* **WSMAF** (within-sample _MAJOR_ allele frequency) -- proportion of reads at this site supporting the most-abundant allele, whichever it is -- float [0,1]
* Site-level quantities (these live in the `INFO` field)
	* **PLAF** (population-level allele frequency) -- average over per-sample WSAFs at this site -- float [0,1]
	* **PLMAF** (population-level _MINOR_ allele frequency) -- min(PLAF, 1 - PLAF) -- float [0,1]
	* **UNW** (unweighted non-reference allele count) -- absolute evidence for non-reference allele -- int [0,inf)
* Sample-level quantities (don't live in VCF):
	* **_F_<sub>ws</sub>** (pseudo-inbreeding coefficient within a sample) -- normalized difference between observed WSAF and expected WSAF if complex infections occurred Hardy-Weinberg-style -- float [0,1]

These were introduced (I think) in [Manske _et al_ (2012) _Nature Genetics_ **487**: 375-379 (PMC3738909)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3738909/).

To avoid spawning not-a-number warnings, `vcfdo` encodes missing values for these quantities as -1.

## Order of operations
The design of `vcfdo` is modular: multiple tools may need to be strung together to accomplish a given task. In most cases this comes down to creating and then using site- or genotype-level annotations in the VCF. The table below summarises the annotations required and/or produced by each tool in `vcfdo`. Annotations in the **recalculates** column are re-computed on the fly _before_ doing any work.

| tool | requires | adds | recalculates
--- | --- | --- | ---
| `thin` | none | none | none
| `prune` | none | none | none
| `wsaf` | `FORMAT/AD` | `FORMAT/WSAF`, `FORMAT/WSMAF`, `INFO/PLAF`, `INFO/PLMAF` | none
| `fws` | `FORMAT/WSAF` | none | `INFO/PLMAF`
| `pca` | `FORMAT/WSAF` | none | `INFO/PLMAF`
| `dist` | `FORMAT/WSAF` | none | none
| `ibs` | none | none | none
| `polarize` | none | `INFO/AA` | none
| `titv` | `INFO/AA` | `INFO/Transversion`, `INFO/StrongWeak`, `INFO/BGC` | none
| `derived` | `INFO/AA` | none | none
| `private` | `FORMAT/WSAF` (with `--reads`) | `INFO/{flag}`, where `{flag}` is specified by user | none
| `sfs` | `INFO/AA` (unless `--ref-is-ancestral`), `FORMAT/WSAF` (with `--reads`) | none | none

Given the constraints implied above, a reasonable compromise between time wasted on re-calculation and space wasted on a bloated VCF would be to run `vcfdo wsaf` and `vcfdo polarize` on a "master" call set and save the result. Then, since `INFO/AA` (ancestral allele) and `FORMAT/WSAF` (within-sample allele frequency) are invariant to subsetting on rows or columns, all other analyses can be done on a stream later.

## Reproducibility
Most operations in `vcfdo` are deterministic up to floating-point weirdness.

The allele-count sampler in `vcfdo sfs` uses the Mersenne twister RNG in Numpy. The seed can be set on the command line with `--seed`, and is fixed to `1` by default -- which means successive calls with the same input will return the same output. That's good for reproducibility but definitely _not_ the desired behavior if one actually wants to sample from the SFS. To allow auto-seeding of the RNG (ie. allow different output each run) pass a seed strictly less than 0.

## Bug reports
If you uncover an error or unexpected behavior, of which there will be plenty, please file a Github issue and contact the author via Slack.
