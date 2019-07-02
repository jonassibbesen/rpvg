# vgprob
Method for calculating variation graph read mapping probabilities given a set of paths in a GBWT index.

### Compilation
*vgprob* requires that [protobuf](https://github.com/protocolbuffers/protobuf) and OpenMP is installed before compilation. 

1. `git clone https://github.com/jonassibbesen/vgprob.git`
2. `cd vgprob`
3. `mkdir build && cd build`
4. `cmake ..`
5. `make -j <threads>`
