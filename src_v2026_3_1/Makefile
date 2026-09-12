# Vanadis 3D v3.1 - Linux/WSL CUDA + GNU Fortran build
# Override CUDA_PATH or GPU_ARCH on the make command line if needed.

CUDA_PATH ?= /usr/local/cuda
FC        := gfortran
NVCC      := $(CUDA_PATH)/bin/nvcc

TARGET    ?= vanadis
FSRC      ?= VANADIS.FOR
CSRC      ?= cudatest.cu
FOBJ      := $(FSRC:.FOR=.o)
COBJ      := $(CSRC:.cu=.o)

# GTX 1060 = compute capability 6.1.  Change, e.g. GPU_ARCH=sm_86,
# when building for another NVIDIA GPU.
GPU_ARCH  ?= sm_61

FFLAGS    ?= -O3 -fopenmp -ffixed-line-length-none
NVCCFLAGS ?= -O3 -arch=$(GPU_ARCH)
LDLIBS    ?= -lgfortran

.PHONY: all run clean debug

all: $(TARGET)

$(FOBJ): $(FSRC)
	$(FC) $(FFLAGS) -c $< -o $@

$(COBJ): $(CSRC)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(TARGET): $(FOBJ) $(COBJ)
	$(NVCC) $(NVCCFLAGS) -o $@ $^ $(LDLIBS) -Xcompiler -fopenmp

run: $(TARGET)
	./$(TARGET)

# Runtime-check build of the Fortran source.  Link still uses CUDA object.
debug: FFLAGS = -O0 -g -fopenmp -ffixed-line-length-none \
                 -fcheck=all -ffpe-trap=invalid,zero,overflow \
                 -finit-real=snan -fbacktrace
debug: NVCCFLAGS = -O0 -g -G -arch=$(GPU_ARCH)
debug: clean all

clean:
	rm -f $(FOBJ) $(COBJ) $(TARGET)
