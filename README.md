# rpvg

Method for inferring path posterior probabilities and abundances from pangenome graph read alignments. For each paired-end read mapped to a [vg](https://github.com/vgteam/vg) pangenome graph (in xg format) the probability of it originating from each path (in a [GBWT](https://github.com/jltsiren/gbwt) index) is calculated. The mapping score and fragment length distribution is used when calculating this probability, and the mapping quality is converted into a separate "noise" probability. Furthermore, paths that share reads with positive probability are clustered into the same group. The read-path probabilities are used to calculate posterior probabilities and infer abundances for each group independently. The method supports mapped reads in both the standard vg alignment format (.gam) and the [vg mpmap](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-vg-mpmap) multipath alignment format (.gmap). 


### Compilation

*rpvg* requires that [CMake](https://cmake.org) (3.10 or higher), [protobuf](https://github.com/protocolbuffers/protobuf) (version 3) and OpenMP are installed before compilation. 

1. `git clone --recursive https://github.com/jonassibbesen/rpvg.git`
2. `cd rpvg && mkdir build && cd build`
3. `cmake ..`
4. `make -j <threads>` 

Compiling *rpvg* should take 5-10 minutes using 4 threads (`-j`). *rpvg* has been successfully built and tested on Linux (CentOS Linux 7 with GCC 8.1.0 and Ubuntu 18.04 with GCC 7.5.0) and Mac (macOS 10.14.6 with Clang 10.0.1). 

### Docker container

A docker container of the latest commit to master is available on [Docker Hub](https://hub.docker.com/repository/docker/jsibbesen/rpvg) and [Quay](https://quay.io/repository/jsibbesen/rpvg). 

### Running rpvg

*rpvg* requires the following five arguments:

```
rpvg -g graph.xg -p paths.gbwt -a alignments.gamp -o rpvg_results -i <inference-model>
```

The prefix used for all output files are given using `-o`. The number of threads can be given using `-t`. 

A small example dataset containing 36,120 haplotype-specific transcripts and 100,000 read pairs is available under [example](https://github.com/jonassibbesen/rpvg/tree/master/example). To infer the expression of the haplotype-specific transcripts in the small pantranscriptome using 4 threads use the following command within the *example* folder:

```
../bin/rpvg -t 4 -g graph.xg -p pantranscriptome.gbwt -f <(zcat pantranscriptome.txt.gz) -a mpmap_align.gamp -o rpvg --inference-model haplotype-transcripts
```

This should take less than a minute to run and will create two files: 

* *rpvg.txt*: Contains the estimated haplotype probability and transcript expression for each haplotype-specific transcript in the pantranscriptome.
* *rpvg_haps.txt*: Contains the estimated probability of each haplotype combination (e.g. diplotype) for all transcripts (only combinations with a probability above zero are shown). 

#### Paths:

The paths to be used for inference should be compressed and indexed using the [GBWT](https://github.com/jltsiren/gbwt). For transcriptome analyses a GBWT with transcript paths can be created using the `vg rna` subcommand in [vg](https://github.com/vgteam/vg). See [Transcriptomic analyses](https://github.com/vgteam/vg/wiki/Transcriptomic-analyses) wiki on the vg github for more information on how to use `vg rna`. 

To decrease the computation time of rpvg it is recommended that a [r-index](https://github.com/jltsiren/gbwt/wiki/Fast-Locate) of the paths is supplied together with the GBWT index. The `vg gbwt` subcommand in vg can be used to construct the r-index from a GBWT index (see [VG GBWT Subcommand](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand) wiki on the vg github). The name of the r-index should be the same as the GBWT index with an added *.ri* extension (e.g. *paths.gbwt.ri*).

#### Inference models:

The method currently contains four different inference models. Each model have been written with a particular path type and corresponding inference problem in mind:

* `haplotypes`: Infers posterior probabilities of haplotype combinations (e.g. diplotypes). For diploid ploidy it uses a branch and bound like algorithm to infer the most probable combination of haplotypes. A faster, less accurate Gibbs sampling scheme can be enabled using `--use-hap-gibbs`, which scales much better for higher ploidies. The maximum ploidy can be given using `-y`.

* `transcripts`: Infers abundances using a Expectation Maximization (EM) algorithm.

* `strains`: Infers abundances using a combination of weighted minimum path cover and EM. **Note that this algorithm has not yet been properly evaluated**.

* `haplotype-transcripts`: Infers abundances using a combination of haplotype sampling and EM. The most probable haplotype combinations are inferred using the same algorithm as used in the `haplotypes` inference model. By default, equivalent haplotypes are inferred for all clustered transcripts, however independent inference of haplotypes for each transcript can be enabled using the `--ind-hap-inference` option. The inference algorithm requires a file (`-f`) containing the haplotype and transcript origin of each path (`--write-info` output from `vg rna`). 

#### Alignment types:

* Use `-e` if the reads are from a strand-specific protocol. Supports the following library types: 
  * Use value `fr` when the first read is from the forward strand and the second read is from the reverse strand.
  * Use value `rf` when the first read is from the reverse stand and the second read is from the forward strand.
* Use `-u` if the input alignments are single-path (*.gam*) format. Default is multipath alignments from [vg mpmap](https://github.com/vgteam/vg/wiki/Multipath-alignments-and-vg-mpmap).
* Use `-s` for single-end reads. Note that the fragment length distribution will still be used for calculating the effective path length.
* Use `-l` for single-molecule long-reads. This is identical to the single-end mode (`-s`), but does not use effective path length normalization.

#### Fragment length distribution:

The fragment length distribution parameters are learned by rpvg. However, in order to learn this the maximum expected fragment length is needed. This is calculated from the expected fragment length distribution mean and standard deviation, which can be given using `-m` and `-d`, respectively. If these are not given the method will look for the parameters in the alignment file and pick the first values that it finds. The input parameters (`-m` and `-d`) are overwritten by the values estimated by rpvg when calculating the read-path probabilities. When the input is single-end reads (`-s`) the expected mean (`-m`) and standard deviation (`-d`) is required as it can not be estimated by rpvg and is needed for the effective path length calculation.
