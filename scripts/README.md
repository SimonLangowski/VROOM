## GPU

## Installation requirements

Python dependencies: nvidia-cuda-toolkit
E.g using conda
`conda install -c conda-forge cuda-python`
`conda install -c "nvidia/label/cuda-12.2.1" cuda-toolkit`

Ensure CGBN is downloaded:
`git submodule update --init --recursive`

Ensure you have GPU drivers, nvcc installed
`nvidia-smi`
`nvcc --verison`

## Benchmarking

To run latency foccussed GPU benchmarks:
`python latency.py`

The python script performs the following tasks
1. Calculate various GPU configurations and their required parameters for CGBN and RNS Montgomery multiplication
2. The file `bench_cgbn.py` contains CGBN specific functions and `gpu.py` contains RNS Montgomery specific functions.  Each is devided into a `gen` func that creates data following the specific parameters, a `result` func that allocates data for results, and a `correct` func that checks the results.
3. The file `basic_benchmarker.py` calls the `gen` and `result` functions, and moves the corresponding data to the GPU.  It then runs the kernel and retrieves the results.  The results are passed to the `correct` function and data about performance is returned.

### Graph generation

Data can be found in the scripts/data/ folder

`python gpugraph.py` to generate graphs for the paper.