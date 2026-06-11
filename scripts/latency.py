from run_on_params import print_results, run_wrapped
from bench_cgbn import run_cgbn_benchmark
from gpu import run_rns_latency_benchmark, VECTOR
from create_graph import gpu_latency_plot

# From the AI google
def next_power_of_2(n):
    if n <= 0:
        return 1
    if n & (n - 1) == 0:
        return n
    return 1 << n.bit_length()

def run_threadwise_rns(modbits, tpi=32, stride=1):
    blocksize = max(32, tpi*stride)
    args = {
        "blocksize": blocksize,
        "blocksper": 1,
        "total": max(1, blocksize // tpi),
        "btpi": tpi,
        "modbits": modbits,
        "variation": VECTOR,
        "m": 1, # allow total 1?
        "stride": stride
    }
    run_wrapped(args, [], run_rns_latency_benchmark)

def run_latencies(func=run_cgbn_benchmark):
    blocksize = 32
    blocksper = 1
    results = []
    maxbits = 4096
    for modbits in range(64, maxbits + 31, 32):
        mod_words = modbits // 32
        word_needed = mod_words
        if func == run_rns_latency_benchmark:
            word_needed += 2
            modbits -= 4
        tpi = next_power_of_2(word_needed)
        # if func == run_rns_benchmark:
        #     tpi = 4
        if func == run_cgbn_benchmark:
            tpi = min(tpi, 32)
            if mod_words % tpi != 0:
                continue
        for stride in [1, 2]:
            if func == run_rns_latency_benchmark:
                blocksize = max(32, tpi*stride)
            min_total = max(1, blocksize // tpi)
            args = {
                "blocksize": blocksize,
                "blocksper": blocksper,
                "total": min_total,
                "btpi": tpi,
                "modbits": modbits,
                "pow": modbits, # simulate constant time exponentiation?
                "m": 8, # only used for tpi matrix stuff
            }
            if func == run_rns_latency_benchmark:
                args["m"] = 1
                args["variation"] = VECTOR
                args["stride"] = stride
            # elif func == run_rns_benchmark:
            #     args['reduction-method'] = 1
            
            if stride == 1 or func == run_rns_latency_benchmark:
                # if func == run_cgbn_benchmark and tpi == 1:
                #     args['variation'] == -1
                #     run_wrapped(args, results, run_rns_benchmark)
                # else:
                run_wrapped(args, results, func)
    
    return results


if __name__ == "__main__":
    # run_threadwise_rns(64, 8, 2)
    # run_threadwise_rns(64, 8, 4)

    
    # run_threadwise_rns(960, 64, 1)
    # run_threadwise_rns(1536, 64, 1)

    # results = run_latencies(run_rns_benchmark)

    results2 = run_latencies(run_cgbn_benchmark)
    results = run_latencies(run_rns_latency_benchmark)
    results.extend(results2)
    fn = print_results(results)
    gpu_latency_plot(fn)