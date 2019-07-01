# vgprob
Method for calculating variation graph read mapping probabilities given a set of paths in a GBWT index.

### Compilation
vgprob depends on the following repositories:

* [cxxxopts](https://github.com/jarro2783/cxxopts)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
* [gbwt](https://github.com/jltsiren/gbwt)
* [protobuf](https://github.com/protocolbuffers/protobuf)
* [libvgio](https://github.com/vgteam/libvgio)
* [gssw](https://github.com/vgteam/gssw)

Remember to update the include and link directories in `CMakeLists.txt` so that they correspond to the location of depedencies on your machine. 
