# rpvg
Method for inferring path posterior probabilities and abundances from variation graph read aligments. For each paired-end read mapped to a [vg](https://github.com/vgteam/vg) variation graph in xg format the probability of it originating from each path in a [GBWT](https://github.com/jltsiren/gbwt) index is calculated. The mapping score and fragment length distribution is used when calculating this probability, and the mapping quality is converted into a seperate "noise" probability. Furthermore, paths that share nodes or reads with positive probability are clustered into the same group. The read-path probabilities are used to calculate posterior probabilities and infer abundances for each group indepedently. The method supports mapped paired-end reads in both the *vg map* [Alignment](https://github.com/vgteam/libvgio/blob/a369fb1f293545eccfdf2d6d3bd4a30b6f5ec664/deps/vg.proto#L111) format (.gam) and *vg mpmap* [MultipathAlignment](https://github.com/vgteam/libvgio/blob/a369fb1f293545eccfdf2d6d3bd4a30b6f5ec664/deps/vg.proto#L156) format (.gmap). 


### Compilation
*rpvg* requires that [protobuf](https://github.com/protocolbuffers/protobuf), [htslib](https://github.com/samtools/htslib) and OpenMP are installed before compilation. 

1. `git clone --recursive https://github.com/jonassibbesen/rpvg.git`
2. `cd rpvg`
3. `mkdir build && cd build`
4. `cmake ..`
5. `make -j <threads>`

### Docker container
A docker container of the latest commit to master is available on [Docker Hub](https://hub.docker.com/repository/docker/jsibbesen/rpvg).

### Running rpvg
*rpvg* requires the following four arguments:
```
rpvg -g graph.xg -p paths.gbwt -a alignments.gam -i <inference-model> > abundances.txt
```
The number of threads can be given using `-t`.

#### Inference models:
The method currently contains four different inference models. Each model have been written with a particurlar path type and corresponding inference problem in mind:

* `haplotypes`: Infers haplotype/diplotype/... posterior probabilities. By default it uses a Gibbs sampling scheme to infer the probabilities, however exact inference can be enabled using `-j`. The exact inference scales exponentially in the sample ploidy and it is therefore only recommended for when the number of haplotypes are low or the ploidy is 1 or 2. The ploidy can be given using `-y`.

* `transcripts`: Infers abundances using a Expectation Maximization (EM) algorithm.

* `strains`: Infers abundances using a combination of weighted minimim path cover and EM. **Note that this algorithm is work in progress and have therefore not been properly evaluated yet**.

* `haplotype-transcripts`: Infers abundances using a combination of haplotype sampling and EM. By default it uses Gibbs sampling to estimate the most probable haplotype combinations (see `haplotypes` for more information). The algorithm requires a file (`-f`) containing the transcript origin of each path (`--write-info` output from *vg rna*).

#### Alignment types:
* Use `-u` if the input alignment format is multipath (*.gamp*) from *vg mpmap*.
* Use `-s` for single-end reads. Note that the fragment length distribtion will still be used for calculating the effective path length.
* Use `-l` for single-molecule long-reads. This is identical to the single-end mode (`-s`), but does not use the effective path length.

#### Fragment length distribution:
Fragment length distribution mean and standard deviation can be given using `-m` and `-d`, respectively. If these are not given the method will look for the parameters in the alignment file and pick the first values that it finds. 

