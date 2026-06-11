import argparse
import ctypes
import sys
import tempfile
import os
if os.name == "posix":
    from cuda import cuda, cudart, nvrtc
import os
import subprocess
import numpy as np
from copy import deepcopy


# test parameters
DEFAULT_BLOCKSIZE = 512
DEFAULT_BLOCKSPER = 1
DEFAULT_ITERATIONS = 100
DEFAULT_TOTAL = 512*32

path = "./"

# https://github.com/sdiehl/gpu-offload/blob/main/Minimal.ipynb

def _cudaGetErrorEnum(error):
    if isinstance(error, cuda.CUresult):
        err, name = cuda.cuGetErrorName(error)
        return name if err == cuda.CUresult.CUDA_SUCCESS else "<unknown>"
    elif isinstance(error, cudart.cudaError_t):
        return cudart.cudaGetErrorName(error)[1]
    elif isinstance(error, nvrtc.nvrtcResult):
        return nvrtc.nvrtcGetErrorString(error)[1]
    else:
        raise RuntimeError('Unknown error type: {}'.format(error))

def checkCudaErrors(result):
    if result[0].value:
        raise RuntimeError("CUDA error code={}({})".format(result[0].value, _cudaGetErrorEnum(result[0])))
    if len(result) == 1:
        return None
    elif len(result) == 2:
        return result[1]
    else:
        return result[1:]
    
# Based on Michael Sullivan's
# bolilerplate code to benchmark things 
  
class CKernelHelper:
    def __init__(self, args, code, devID, additionalDefines=[]):
        with open("sourceCode.cu", "w+") as f:
            f.write(code)
       
        CUDA_HOME = os.getenv('CUDA_HOME') if 'CUDA_HOME' in os.environ else '/usr/local/cuda'
        if CUDA_HOME == None:
            raise RuntimeError('Environment variable CUDA_HOME is not defined')
        include_dirs = os.path.join(CUDA_HOME, 'include')

        # Initialize CUDA
        checkCudaErrors(cudart.cudaFree(0))

        major = checkCudaErrors(cudart.cudaDeviceGetAttribute(cudart.cudaDeviceAttr.cudaDevAttrComputeCapabilityMajor, devID))
        minor = checkCudaErrors(cudart.cudaDeviceGetAttribute(cudart.cudaDeviceAttr.cudaDevAttrComputeCapabilityMinor, devID))
        prefix = 'sm'
        arch_arg = bytes(f'--gpu-architecture={prefix}_{major}{minor}', 'ascii')

        try:
            opts = [b'--fmad=true', arch_arg,
                    '--include-path={}'.format(include_dirs).encode('UTF-8'),
                    b'--std=c++20', b'-dopt=on', b'-O3']
            if (args.dump_sass):
                opts.append(b'-lineinfo')
            for d in additionalDefines:
                opts.extend(["--define-macro", d])
            command = ["nvcc", "sourceCode.cu", "-o", path + "sourceCode.cubin", "--cubin"]
            command.extend(opts)
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as err:
            raise SystemError('subprocess failed with error: {}'.format(err)) from None
        data = self.getCUBIN()
        self.module = checkCudaErrors(cuda.cuModuleLoadData(np.char.array(data)))

    def getCUBIN(self):
        with open(path + "sourceCode.cubin", "rb") as f:
            return f.read()

    def getFunction(self, name):
        return checkCudaErrors(cuda.cuModuleGetFunction(self.module, name))

def parse_args(custom_args=None, software_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--blocksize', type=int, default=DEFAULT_BLOCKSIZE, help='Threads per block (default: %s)' %
                        (DEFAULT_BLOCKSIZE,))
    parser.add_argument('--blocksper', type=int, default=DEFAULT_BLOCKSPER, help='Blocks per SM (default: %s)' %
                        (DEFAULT_BLOCKSPER,))
    # parser.add_argument('--repeat', type=int, default=DEFAULT_REPEAT, help='Repeats per thread (default: %s)' %
    #                     (DEFAULT_REPEAT,))
    parser.add_argument('--dump-sass', action='store_true', default=False, help='Dump sass to file')
    parser.add_argument('--dump-src', action='store_true', default=False, help='Dump replaced source code to file')
    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (default: 0)')
    parser.add_argument('--iterations', type=int, default=DEFAULT_ITERATIONS, help='Number of iterations to run for accurate timing (default: %s)'
                        % (DEFAULT_ITERATIONS,))
    parser.add_argument('--total', type=int, default=DEFAULT_TOTAL, help='Total number of parallel instantiations of the computation (default: %s)' %
                        (DEFAULT_TOTAL,))
    parser.add_argument('--skip-correct', type=bool, default=False, help="Skip correctness check")
    parser.add_argument('--debug', type=bool, default=False, help="set debug flag")
    if custom_args is not None:
        for custom_arg in custom_args:
            arg_flag, arg_type, arg_default, arg_help = custom_arg
            parser.add_argument(arg_flag, type=arg_type, default=arg_default, help=arg_help)
    return parser.parse_args(software_args)


def _get_cmd_output(cmd):
    import subprocess

    try:
        print('Executing:\n{}'.format(cmd))
        return subprocess.check_output(cmd.split()).decode(sys.stdout.encoding)
    except subprocess.CalledProcessError as e:
        raise SystemError('subprocess failed with error: {}'.format(e.output))

CACHE = [None, {}, {}, {}]

def destroy_cache():
    [device_data, device_output, stream, start, stop, host_raw] = CACHE[0]
    [checkCudaErrors(cudart.cudaFree(d_A)) for d_A in device_data]
    checkCudaErrors(cudart.cudaFree(device_output))
    checkCudaErrors(cudart.cudaStreamDestroy(stream))
    checkCudaErrors(cudart.cudaEventDestroy(start))
    checkCudaErrors(cudart.cudaEventDestroy(stop))
    CACHE[0] = None
    CACHE[1] = {}
    CACHE[2] = {}

# Use cached data if all other arguments are equal
CACHE_IGNORE = ["blocksper", "pow"]
def is_cached(args):
    args_dict = vars(args)
    CACHE[3] = deepcopy(args_dict)
    cached_dict = CACHE[2]
    if CACHE[0] == None:
        print("cache empty")
        return False
    if len(args_dict) != len(cached_dict):
        print(args_dict, cached_dict)
        return False
    for k, v in args_dict.items():
        if k not in CACHE_IGNORE:
            if k not in cached_dict or v != cached_dict[k]:
                print(f"Mismatch on key {k} {v}")
                return False
    return True

def store_cache(args, device_data, device_output, stream, start, stop, host_raw):
    CACHE[0] = [device_data, device_output, stream, start, stop, host_raw]
    CACHE[1] = deepcopy(vars(args))
    CACHE[2] = CACHE[3]

def load_cache(args):
    for k, v in CACHE[1].items():
        if k not in CACHE_IGNORE:
            setattr(args, k, v)
    return CACHE[0]

# Cuda does not allow 0 length arrays
def sanitize(arrs):
    for i in range(len(arrs)):
        if arrs[i].nbytes == 0:
            arrs[i] = np.zeros(1, dtype=np.uint32)
    return arrs

# display_name: displayed during printing
# file_name: cuda file name
# custom_args: additional passed to arg parse
# gen_func: takes in args and returns data for cuda as a list of numpy arrays, and data for checking
# result_func: takes in args and returns a numpy array capable of holding the function output, the required shared memory size, and the number of elements processed per block
# cuda_name: name of cuda function
# (optional) correct_func: function that takes data for checking and output of cuda to determine correctness
# (optional) software_args: for passing arguments directly from another python file (rather than reading command line)
# Return: time taken
def run_benchmark(display_name, file_name, cuda_name, custom_args, gen_func, result_func, correct_func=None, software_args=None, use_cache=False):
    args = parse_args(custom_args, software_args)
    checkCudaErrors(cudart.cudaSetDevice(args.gpu))
    properties = checkCudaErrors(cudart.cudaGetDeviceProperties(args.gpu))
    if use_cache and is_cached(args):
        print("Using cached data")
        [device_data, device_output, stream, start, stop, host_raw] = load_cache(args)
        host_output, shared_size, elements_per_block = result_func(args)
    else:
        if CACHE[0] is not None:
            destroy_cache()
        # Generate host side data
        print("Generating python data")
        host_data, host_raw = gen_func(args)
        host_data = sanitize(host_data)
        print("Moving data to GPU")
        # Allocate device data
        device_data = [checkCudaErrors(cudart.cudaMalloc(a.nbytes)) for a in host_data]
        host_output, shared_size, elements_per_block = result_func(args)
        device_output = checkCudaErrors(cudart.cudaMalloc(host_output.nbytes))

        # boilerplate
        stream = checkCudaErrors(cudart.cudaStreamCreateWithFlags(cudart.cudaStreamNonBlocking))
        start = checkCudaErrors(cudart.cudaEventCreate())
        stop = checkCudaErrors(cudart.cudaEventCreate())

        # Copy from host to device
        for (h_A, d_A) in zip(host_data, device_data):
            checkCudaErrors(cudart.cudaMemcpy(d_A, h_A, h_A.nbytes, cudart.cudaMemcpyKind.cudaMemcpyHostToDevice))
        store_cache(args, device_data, device_output, stream, start, stop, host_raw)

    checkCudaErrors(cudart.cudaMemcpy(device_output, host_output, host_output.nbytes, cudart.cudaMemcpyKind.cudaMemcpyHostToDevice))

    with open(file_name) as f:
        text = f.read()
    
    # To avoid trouble
    prefix = "PYREPLACE_"
    for (arg, val) in sorted(vars(args).items(), key=lambda x : len(x[0]), reverse=True):
        text = text.replace(prefix + arg.upper(), str(val))
    
    kernelHelper = CKernelHelper(args, text, args.gpu)
    if args.dump_src:
        with open('%s_GENERATED.cuh' % cuda_name, 'w') as g:
            g.write(cuda_name)

    if args.dump_sass:
        with tempfile.NamedTemporaryFile(suffix='.cubin') as f:
            f.write(kernelHelper.getCUBIN())
            f.flush()
            with open('%s.sass' % cuda_name, 'w') as g:
                print('Generating SASS file...')
                g.write(_get_cmd_output('nvdisasm -g {}'.format(f.name)))

    cudaFunc = kernelHelper.getFunction(cuda_name.encode('UTF-8'))

    # Expand shared memory size if needed
    if shared_size > 48*1024:
        checkCudaErrors(cuda.cuFuncSetAttribute(cudaFunc, cuda.CUfunction_attribute.CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES, shared_size))

    threads = cudart.dim3()
    threads.x = args.blocksize
    threads.y = 1
    threads.z = 1

    if (args.total % elements_per_block != 0):
        raise ValueError(f"Pad total to multiple of {elements_per_block} ({args.total - (args.total % elements_per_block) + elements_per_block})")
    block_chunks_required = args.total // elements_per_block

    # assert(args.total >= args.blocksper * args.blocksize)
    # assert(args.total % (args.blocksper * args.blocksize) == 0)
    num_blocks = min(args.blocksper * properties.multiProcessorCount, block_chunks_required)
    base_chunks_per_block = block_chunks_required // num_blocks

    grid = cudart.dim3()
    grid.x = num_blocks
    grid.y = 1
    grid.z = 1

    arguments = []
    argumentTypes = []
    arguments.extend(device_data)
    argumentTypes.extend([ctypes.c_void_p for _ in device_data])
    arguments.append(device_output)
    argumentTypes.append(ctypes.c_void_p)
    arguments.append(int(block_chunks_required))
    argumentTypes.append(ctypes.c_int)
    arguments.append(int(base_chunks_per_block))
    argumentTypes.append(ctypes.c_int)
    kernelArguments = (tuple(arguments), tuple(argumentTypes))

    print(f"Blocks: {threads} grid {grid}")

    checkCudaErrors(cudart.cudaEventRecord(start, stream))
    for _ in range(args.iterations):
        checkCudaErrors(cuda.cuLaunchKernel(cudaFunc, grid.x, grid.y, grid.z,
          threads.x, threads.y, threads.z, shared_size, stream, kernelArguments, 0))
    checkCudaErrors(cudart.cudaEventRecord(stop, stream))
    checkCudaErrors(cudart.cudaEventSynchronize(stop))
    msecTotal = checkCudaErrors(cudart.cudaEventElapsedTime(start, stop))
    print("{:s}: Run time: {:,.3f} ms for {:d} iterations of {:d} elements".format(display_name, msecTotal, args.iterations, args.total))

    # Copy result from device memory to host memory
    checkCudaErrors(cudart.cudaMemcpy(host_output, device_output, host_output.nbytes, cudart.cudaMemcpyKind.cudaMemcpyDeviceToHost))

    if correct_func and not args.skip_correct:
        print("checking results...")
        correct = correct_func(args, host_raw, host_output)
        print(correct)
    else:
        correct = "not checked"

    # Cleanup
    # Not using cache; free data
    if not use_cache:
        destroy_cache()
    return {"time": msecTotal,
            "name": display_name,
            "func": cuda_name,
            "shared_mem": shared_size,
            "iterations": args.iterations,
            "total": args.total,
            "blocksize": args.blocksize,
            "blocksper": args.blocksper,
            "args": args, # for task dependent arguments
            "correct": correct
            }