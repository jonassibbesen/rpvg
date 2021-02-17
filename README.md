# rpvg

Method for inferring path posterior probabilities and abundances from variation graph read alignments. For each paired-end read mapped to a [vg](https://github.com/vgteam/vg) variation graph (in xg format) the probability of it originating from each path (in a [GBWT](https://github.com/jltsiren/gbwt) index) is calculated. The mapping score and fragment length distribution is used when calculating this probability, and the mapping quality is converted into a separate "noise" probability. Furthermore, paths that share reads with positive probability are clustered into the same group. The read-path probabilities are used to calculate posterior probabilities and infer abundances for each group independently. The method supports mapped paired-end reads in both the *vg map* [Alignment](https://github.com/vgteam/libvgio/blob/a369fb1f293545eccfdf2d6d3bd4a30b6f5ec664/deps/vg.proto#L111) format (.gam) and *vg mpmap* [MultipathAlignment](https://github.com/vgteam/libvgio/blob/a369fb1f293545eccfdf2d6d3bd4a30b6f5ec664/deps/vg.proto#L156) format (.gmap). 


### Compilation

*rpvg* requires that [protobuf](https://github.com/protocolbuffers/protobuf) and OpenMP are installed before compilation. 

1. `git clone --recursive https://github.com/jonassibbesen/rpvg.git`
2. `cd rpvg`
3. `mkdir build && cd build`
4. `cmake ..`
5. `make -j <threads>` 

### Docker container

A docker container of the latest commit to master is available on [Docker Hub](https://hub.docker.com/repository/docker/jsibbesen/rpvg) and [Quay](https://quay.io/repository/jsibbesen/rpvg). 

### Running rpvg

*rpvg* requires the following five arguments:

```
rpvg -g graph.xg -p paths.gbwt -a alignments.gam -o rpvg_results -i <inference-model>
```

The prefix used for all output files are given using `-o`. The number of threads can be given using `-t`. 

#### Paths:

The paths to be used for inference should be compressed and indexed using the [GBWT](https://github.com/jltsiren/gbwt). To decrease the computation time of *rpvg* it is recommeded that a [r-index](https://github.com/jltsiren/gbwt/wiki/Fast-Locate) of the paths is supplied together with the GBWT index. `vg gbwt` as part of the [vg toolkit](https://github.com/vgteam/vg) can be used to construct the r-index from a GBWT index. The name of the r-index should be the same as the GBWT index with an added *.ri* extension (e.g. *paths.gbwt.ri*).

#### Inference models:

The method currently contains four different inference models. Each model have been written with a particular path type and corresponding inference problem in mind:

* `haplotypes`: Infers haplotype/diplotype/... posterior probabilities. For diploid ploidy it uses a branch and bound like algorithm to infer the most probable combination of haplotypes. A faster, less accurate Gibbs sampling scheme can be enabled using `--use-hap-gibbs`, which scales much better for higher ploidies. The maximum ploidy can be given using `-y`.

* `transcripts`: Infers abundances using a Expectation Maximization (EM) algorithm.

* `strains`: Infers abundances using a combination of weighted minimum path cover and EM. **Note that this algorithm has not yet been properly evaluated**.

* `haplotype-transcripts`: Infers abundances using a combination of haplotype sampling and EM. The most probable haplotype combinations are inferred using the same algorithm as used in the `haplotypes` inference model. By default, equivalent haplotypes are inferred for all clustered transcripts, however independent inference of haplotypes for each transcript can be enabled using the `--ind-hap-inference` option. The inference algorithm requires a file (`-f`) containing the transcript origin of each path (`--write-info` output from *vg rna*). 

#### Alignment types:

* Use `-e` if the reads are from a strand-specific protocol. Supports the following library types: 
  * Use value `fr` when the first read is from the forward strand and the second read is from the reverse strand.
  * Use value `rf` when the first read is from the reverse stand and the second read is from the forward strand.
* Use `-u` if the input alignments are single-path (*.gam*) format. Default is multipath alignments from *vg mpmap*.
* Use `-s` for single-end reads. Note that the fragment length distribution will still be used for calculating the effective path length.
* Use `-l` for single-molecule long-reads. This is identical to the single-end mode (`-s`), but does not use effective path length normalization.

#### Fragment length distribution:

The fragment length distribution parameters are learned by *rpvg*. However, in order to learn this the maximum expected fragment length is needed. This is calculated from the expected fragment length distribution mean and standard deviation, which can be given using `-m` and `-d`, respectively. If these are not given the method will look for the parameters in the alignment file and pick the first values that it finds. The input parameters (`-m` and `-d`) are overwritten by the values estimated by *rpvg* when calculating the read-path probabilities. When the input is single-end reads (`-s`) the expected mean (`-m`) and standard deviation (`-d`) is required as it can not be estimated by *rpvg* and is needed for the effective path length calculation.
