# SeleScan: Guide

*By Pavel Salazar-Fernandez (epsalazarf@gmail.com), in collaboration with Marc Pybus (marcpybus@gmail.com) and Andres Moreno-Estrada (andres.moreno@cinvestav.mx)*

*Laboratory of Human Population and Evolutionary Genomics (LANGEBIO) & Institut de Biologia Evolutiva (PRBB)*

## About
This document explains how to set up a selection scan using this repository. Phased VCF files are split into chromosomes and summary statistics are calculated with the R packages `PopGenome` and `rEHH`. The results are compiled into a table and then compared to results derived from neutral simulated data to identify outlier regions that may have been affected by selective pressure.

## Overview
1. Split phased VCF file into chromosomes.
2. Obtain summary statistics with `PopGenome`: neutrality, Fst, SFS.
3. Obtain haplotype statistics with `eRHH`: iHS, iHH.
4. Generate neutral simulations with `ms`.
5. Obtain summary and haplotype statistics (steps 2 and 3) for the simulations.
6. Merge all results in a summary table.
7. Scan table for regions with significant stats suggesting selection.

## Requirements
**Software:**
* Unix OS: bash
* [Perl 5+](http://www.perl.org/)
* [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html)
* [vcftools](http://vcftools.github.io/index.html)
* [R 3.0+](http://cran.r-project.org/)
    * Packages: [PopGenome](https://cran.r-project.org/web/packages/PopGenome/index.html), [rEHH](https://cran.r-project.org/web/packages/rehh/index.html), and dependencies.

## Simulations
**Recommended reading:**
> Hoban, S., Bertorelle, G., & Gaggiotti, O. E. (2011). Computer simulations: tools for population and evolutionary genetics. Nat Rev Genet, 13(2), 110–122. http://doi.org/10.1038/nrg3130

### MS
[OmicsTools: ms](http://omictools.com/ms-tool)

A Monte Carlo computer program is available to generate samples drawn from a population evolving according to a Wright-Fisher neutral model. The program assumes an infinite-sites model of mutation, and allows recombination, gene conversion, symmetric migration among subpopulations, and a variety of demographic histories. The samples produced can be used to investigate the sampling properties of any sample statistic under these neutral models.

**Related software:**
* [msHOT](http://omictools.com/mshot-tool) (Hellenthal and Stephens, 2007): Both crossover and gene conversion hotspots have been incorporated to *ms*.
* [MSMS](http://omictools.com/msms-tool) (Ewing and Hermisson, 2010): Includes the functionality of *ms* to model population structure and demography, but adds a model for deme- and time-dependent selection using forward simulations.

**Usage:**
*Basic:*
`ms nsam nreps -t θ > output.ms`

* `nsam`: number of copies of the locus in each sample.
* `nreps`: number of independent samples to generate.
* `-t θ`: mutation parameter (*θ= 4 ∗ N ∗ μ ∗ Sequence Length*, where μ is the mutation rate).

*Advanced:*
`ms nsam nreps -t θ -r ρ -I POPS N1 N2 -ej t i j -en t i x > output.ms`

* `-r ρ`: cross-over parameter (*ρ= 4 ∗ N ∗ r ∗ Sequence Length*, where r is the recombination rate).
* `-I pops n1 n2`: number of populations and number of locus for each (must add exactly to `nsam`).
* `-ej t i j`: moves all lineages from population *i* to population *j* at time *t* with migration rate zero.
* `-en t i x`: sets population i size to *x* ∗ *N* at time *t* and growth rate zero.

**Output:**
A text file containing all the generated haplotypes generated. An example output is:

```
./ms 3 3 -t 3.0
56917 16239 12250

//
segsites: 2
positions: 0.4493 0.7948
01
00
10

//
segsites: 4
positions: 0.2499 0.4466 0.6762 0.8853
0100
0011
1000

//
segsites: 0
```
Where:
* Line 1: *ms* command used.
* Line 2: Random number seed(s) used. Values saved as `seedms`.
* `//`: Set of lines separating each simulation.
* `segsites: #`: Segregating sites generated. Nothing else printed if zero.
* `positions: ### ### ### ...`: Positions for each polymorphic site, on a scale (0,1).
* Haplotypes generated, `0` represents ancestral allele and `1` derived allele.

## Population Genomics Software

### PopGenome
[CRAN: PopGenome](https://cran.r-project.org/web/packages/PopGenome/index.html)

PopGenome is a new package for population genomic analyses and method development. PopGenome includes, e.g., a wide range of polymorphism, neutrality statistics, and FST estimates; these can be applied to sequence data stored in alignment format, as well as to whole genome SNP data, e.g., from the 1000/1001 Genome projects. The full range of methods can be applied to whole alignments, sets of sub-sequences, and sliding windows based on either nucleotide positions or on SNP counts. PopGenome is also able to handle GFF/GTF annotation files and automatically specifies the SNPs located in, e.g., exon or intron regions. Those subsites can be analyzed together (e.g., all introns together) or each region seperately (e.g., one value per intron). The PopGenome frame- work is linked to Hudson’s MS and Ewing’s MSMS programs for significance tests using coalescent simulations.

### rEHH
[CRAN: rEHH](https://cran.r-project.org/web/packages/rehh/index.html)

Functions for the detection of footprints of selection on dense SNP data using Extended Homozygosity Haplotype (EHH) based tests. The package includes computation of EHH, iHS (within population) and Rsb (across pairs of populations) statistics. Various plotting functions are also included to facilitate visualization and interpretation of the results.

## Selection Scan Pipeline
### Directory Tree
* `SelectionScan/`
    * `README.md`: this guide.
    * `file_index.txt`: extended file index.
    * `PopGenome/`: neutrality tests and results of empirical data.
    * `rEHH/`: haplotype analyses and results of empirical data.
    * `simulations/`: creation and analysis of simulated data.
        * `PopGenome/`: neutrality tests and results of simulated data.
        * `rEHH/`: haplotype analyses and results of simulated data.

### Empirical Data Analysis
#### Input
* VCF file requirements:
    * All autosomes must be **phased** in a single concatenated file
    * Only  biallelic SNPs*
    * No indels*
    * No monomorphic positions*
    * No singletons*

* Sample IDs files: a list of sample IDs for each population, one `[pop].ids` per population.

*Note*: requirements with \* are filtered by `split_by_*.sh` scripts.

#### Pipeline
##### PopGenome
1. Move or link the VCF file to `PopGenome/`.
2. Modify `split_by_chr.sh` for your data. If needed, add a sample filter.
3. Run `split_by_chr.sh`. Individual directories for each chromosome will be created, containing a VCF file for that chromosome.
4. Run `RUN_PopGenome.sh`. Results for each chromosome will be written to its corresponding directory.

##### rEHH
1. Move or link the VCF file and sample IDs files to `rEHH/`
2. Modify `split_by_chr.sh` for your data, declaring the populations to be analyzed and their respective sample IDs files.
3. Run `split_by_pop_and_chr.sh`. Individual directories for each chromosome will be created, containing a VCF file for each population and the corresponding chromosome.
4. Run `RUN_PIPELINE.sh`. Results for each population and chromosome will be written to its corresponding directory.

#### Output
* **PopGenome** - Subdirectories per chromosome (chr*/) containing:
    * `neutrality.results.[chr].[pop].txt`
    * `SFS.[chr].ALL.txt`
    * `SFS.[chr].[pop].txt`
    * `fst.[chr].ALL.txt`
    * `fst.[chr].[popA]x[popB].txt`
* **rEHH (Version 2.0+ required)** - Subdirectories per chromosome containing:
    * `[pop].[chr].hap`
    * `[pop].[chr].map`
    * `[pop].[chr].ihs`
    * `[pop].[chr].pdf`

##### Sample Output
* `neutrality.results.[chr].[pop].txt`

|*(region)*|Tajima.D|Segregating Sites|Rozas.R_2|Fu.Li.F|Fu.Li.D|Fu.F_S|Fay.Wu.H|Zeng.E|Strobeck.S|
|---|---|---|---|---|---|---|---|---|---|
|### - ###|###|###|NA|###|###|NA|NA*|NA*|NA|
|...|...|...|...|...|...|...|...|...|...|

NA*: Calculated only if ancestral alleles are provided

* `SFS.[chr].txt`

|*(SNP)*|pop 1|pop 2|pop 3|
|---|---|---|---|
|###|###|###|###|
|...|...|...|...|

* `fst.[chr].txt` / `fst.[chr].[popA]x[popB].txt`

|*(SNP)*|*(Fst)*|
|---|---|
|###|###|
|...|...|

* `[pop].[chr].ihs`

|ID|CHR|POSITION|iHS|Pvalue|
|---|---|---|---|---|
|rs###|##|#####|###|###|
|...|...|...|...|...|

* `[pop].[chr].pdf`
    1. iHS distribution (kernel density)
    2. iHS (dual manhattan plot)
    3. P-value (manhattan plot)


### Simulated Data  Analysis
#### Input
##### ms Parameters for ALL

Parameters:
* Loci= 2n (n diploid samples)
* Mutation Rate (μ) = 1.25e-8
* Population Size (N) = ?
* Sequence Length= 100,000 bp
* Theta (*θ = 4 ∗ N ∗ μ ∗ Sequence Length*) = ?
* Recombination rate range (*r=min:max*)= 1e-9, 1e-8
* Rho_min (*ρ= 4 ∗ N ∗ r min ∗ Sequence Length*)= ?
* Rho_max (*ρ= 4 ∗ N ∗ r max ∗ Sequence Length*)= ?
* Populations (number)= POP1 (?), POP2 (?), POP3 (?), POP4 (?)

Command:
`ms [loci] 1 -t [theta] -r $(( ( RANDOM % [rho_max - rho min + 1])  + rho_min )) [seq length +1] -I 4 [pop1] [pop2] [pop3] [pop4] > output.ms`

#### Pipeline
1. Setup the parameters in these files:
    * `ms_commands_100Kb.sh`
    * `add_ancestral.pl`
    * `pop_genome_analysis_simulations.R`
    * `RUN_PIPELINE_2.sh`
    * `run_conversion_and_parsing.pl`
    * `run_conversion_and_parsing2.pl`
2. Run `bash RUN_PIPELINE_1.sh` (~10 hours)
3. Run `bash RUN_PIPELINE_2.sh` (~1 hour)
4. Run `bash RUN_PIPELINE_3.sh` (>4 days)
    * Runtime estimation: (Populations * Sequence Length * Simulations) / 4.16e8 = Hours

#### Output
* **`RUN_PIPELINE_1.sh`**
    * `sim_*.ms`: simulation files.
    * `sim_*.ms.pos`: positions for segregating sites in the simulation file.
    * `neutrality.sim_*.[pop].txt`
    * `SFS.sim_*.ALL.txt`
    * `sim_*.[pop].hap`
    * `sim_*.[pop].map`
    * `sim_*.[pop].ihh`
* **`RUN_PIPELINE_2.sh`**
    * `sim_*.[pop].ihh`: summary iHH values of all simulations for a single population.
* **`RUN_PIPELINE_3.sh`**
    * `sim_*.[pop].ihs`: summary iHS values of all simulations for a single population.
    * `sim_*.[pop].pdf`: plot.

##### Sample Output
* `neutrality.sim_*.[pop].txt`

|*(region)*|Tajima.D|Segregating Sites|Rozas.R_2|Fu.Li.F|Fu.Li.D|Fu.F_S|Fay.Wu.H|Zeng.E|Strobeck.S|
|---|---|---|---|---|---|---|---|---|---|
|### - ###|###|###|NA|###|###|NA|###|###|NA|
|...|...|...|...|...|...|...|...|...|...|

* `SFS.sim_*.ALL.txt` (transposed)

||snp1|snp2|snp3|...|
|---|---|---|---|---|
|pop1|###|###|###|...|
|pop2|###|###|###|...|
|pop3|###|###|###|...|

* `fst.results.sim_*.[pop].txt` / `fst.results.sim_*.[popA]x[popB].txt`

|*(SNP)*|*(Fst)*|
|---|---|
|###|###|
|...|...|

* `sim_*.[pop].ihs`

|ID|CHR|POSITION|iHS|Pvalue|
|---|---|---|---|---|
|chr1_###|##|#####|###|###|
|...|...|...|...|...|


* `[pop].[chr].pdf`
    1. iHS distribution (kernel density)
    2. iHS (dual manhattan plot)
    3. P-value (manhattan plot)

## Summary Statistics Compilation
Since all results are splitted to populations and/or chromosomes, an additional run of scripts is necessary to condense all results into genome-wide files.

### Summary Statistics: Tables

#### Requirements

**Software**
* Unix OS: bash

**Data**
All files produced by correctly-completed analyses from the previous section, not moved nor renamed from their scripted location, and assuming 10,000 simulations.

#### Procedure

1. Change location to the base pipeline directory (`SelectionScan/`)
2. Set up the populations in the `#Parameters` section of these files:
    * `RUN_CONCAT_EMP.sh`
    * `RUN_CONCAT_SIM.sh`
3. Run `bash RUN_CONCAT_EMP.sh`. This will start concatenations for all empirical data (PopGenome and rEHH).
4. Run `bash RUN_CONCAT_SIM.sh`. This will start concatenations for all simulated data (PopGenome and rEHH).

#### Output
The tables generated have the same columns as the input files describe in the **Output** sections of **Selection Scan Pipeline**, except for those columns with only `NA` values.

* `PopGenome/results/`
    * `neut.results.adna.[pop].txt`
    * `SFS.results.adna.[all].txt`
    * `SFS.results.adna.[pop].txt`
    * `fst.results.adna.[all].txt`
    * `fst.results.adna.[pair].txt`
* `rEHH/results/`
    * `[pop].adna.ihs`
* `simulations/PopGenome/concats/`
    * `neut.sims.[pop].txt`
    * `SFS.sims.[all].txt`
    * `SFS.sims.[pop].txt`
    * `fst.sims.[all].txt`
    * `fst.sims.[pair].txt`
* `rEHH/concats/`
    * `[pop].sims.ihs`

### Summary Statistics: Plots

#### Requirements

**Software**
* Unix OS: bash
* [R 3.0+](http://cran.r-project.org/)
    * Packages: [data.table](http://cran.r-project.org/web/packages/data.table/index.html), [qqman](http://cran.r-project.org/web/packages/qqman/index.html), and dependencies.

**Data**
All files produced by correctly-completed concatenations from the previous section, not moved nor renamed from their scripted location.

#### Procedure

1. Set up the parameters in the `#<INPUT>` section of these files:
    * `PopGenome/NEUTplot.R`
    * `PopGenome/SFSplot.R`
    * `PopGenome/SFSplot.empxsim.R`
    * `PopGenome/FSTplot.R`
    * `PopGenome/pairFSTplot.R`
    * `rEHH/iHSplot.R`
    * `rEHH/XPEHHplot.R`
2. Run `bash RUN_PLOT_ALL.sh`. This will generate all plots for the concatenated empirical data.

#### Output
Plots and tables for Top 1% and Top 0.1% for all stats will be produced.

* `PopGenome/results/`
    * `[pop].[neutrality stat].plot.png`: neutrality statistic can be Tajima's D and Fu & Li's F or D.
    * `[all].SFS.bars.png`, `[all].SFS.lines.png`: site frequency spectrum for each population.
    * `SFS.EMPxSIM.abs.png`: Absolute SNP counts for present MAF bins.
    * `estimated.asc.bias`: Estimated ascertainment bias in percentage for each population, comparing the number of simulated versus empirical SNPs for each MAF bin.
    * `[all].Fst.plot.png`: Plot of region Fst value from an all-populations comparison.
    * `[pair].Fst.plot.png`: Plot of region Fst value from paired-populations comparisons.
    *`[stat].[pop/pair/all].top1.txt`,`[stat].[pop/pair/all].top01.txt`: Top 1% or 0.1% results for a stat (neutrality or Fst).

* `rEHH/results/`
    * `[pop].adna.ihs.png`: genome-wide *-log10(p-value)* Manhattan plot for each population.
    * `[pop].adna.ihs.qq.png`: Q-Q plot for empirical and expected p-values for each population.
    * `[pop].ihs.top1.txt`,`[pop].ihs.top01.txt`: Top 1% or 0.1% results for iHS values.
    * `[all].top1.venn.png`,`[all].top01.venn.png`: shared Top 1% or Top 0.1% SNPs for iHS values between populations.
    * `[pair].xpehh.all.txt`: XP-EHH statistics for all sites present in a population pair.
    * `[pair].xpehh.top1.txt`,`[pair].xpehh.top01.txt`: Top 1% or 0.1% results for iHS values.
    * `[pair].xpehh.png`: genome-wide raw XP-EHH score plot for each paired population.
    * `[pair].xpehh.logp.png`: genome-wide *-log10(p-value)* Manhattan plot for each paired population.


##### Sample XP-EHH Output

* `[pair].xpehh.*.txt`

|ID|CHR|POSITION|XPEHH|-log10(p-value) [bilateral]|
|---|---|---|---|---|
|rs###|##|#####|###|###|
|...|...|...|...|...|

---
