# rpvg
Method for calculating read-path probabilities using variation graphs. For each paired-end read mapped to a [vg](https://github.com/vgteam/vg) variation graph the probability of it originating from each path in a [GBWT](https://github.com/jltsiren/gbwt) index is calculated. The mapping score and fragment length distribution is used when calculating this probability, and the mapping quality is converted into a seperate "noise" probability. Furthermore, paths that share reads with positive probability are clustered into the same group. The method supports mapped paired-end reads in both the *vg map* [Alignment](https://github.com/vgteam/libvgio/blob/a369fb1f293545eccfdf2d6d3bd4a30b6f5ec664/deps/vg.proto#L111) format (.gam) and *vg mpmap* [MultipathAlignment](https://github.com/vgteam/libvgio/blob/a369fb1f293545eccfdf2d6d3bd4a30b6f5ec664/deps/vg.proto#L156) format (.gmap). 


### Compilation
*rpvg* requires that [protobuf](https://github.com/protocolbuffers/protobuf), [htslib](https://github.com/samtools/htslib) and OpenMP are installed before compilation. 

1. `git clone --recursive https://github.com/jonassibbesen/rpvg.git`
2. `cd rpvg`
3. `mkdir build && cd build`
4. `cmake ..`
5. `make -j <threads>`


### Running rpvg
*rpvg* requires the following three arguments:
```
rpvg -g graph.xg -p paths.gbwt -a alignments.gam > probs.txt
```
Use `-m` if the input alignment format is multipath (.gamp). Fragment length distribution mean and standard deviation can be given using `-i` and `-d`, respectively. If these are not given the method will look for the parameters in the alignment file and pick the first values that it finds. 
