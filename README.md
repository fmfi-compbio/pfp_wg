Original work:                              \
https://gitlab.com/manzai/Big-BWT           \
https://github.com/simongog/sdsl-lite       \
https://github.com/waYne1337/BWT-Tunneling  \

# How to install and run
1. Install system-wide SDSL according to the instructions
2. Run `git clone https://github.com/andynet/cds_project.git`
3. Open cds\_project/Big-BWT/makefile and set correct `include path-to-sdsl-repo/sdsl-lite/Make.helper`
4. Run `make`
5. Run `./bigbwt -w 4 -p 50 -i data/yeast.raw`

