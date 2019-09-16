# rpvg
Method for calculating Read Probabilities using Variation Graphs.

### Compilation
*rpvg* requires that [protobuf](https://github.com/protocolbuffers/protobuf), [htslib](https://github.com/samtools/htslib) and OpenMP are installed before compilation. 

1. `git clone https://github.com/jonassibbesen/rpvg.git`
2. `cd rpvg`
3. `mkdir build && cd build`
4. `cmake ..`
5. `make -j <threads>`
