import csv
import re
import sys
import glob
import matplotlib.pyplot as plt
import ast
import json
import numpy as np

def trim_dict(d):
    out = {}
    for (k, v) in d.items():
        out[k.strip()] = v.strip()
    return out

def pad_str(s, pad_length):
    return s + ' ' * (pad_length - len(s))

def read_csv(fn):
    outx = {}
    outy = {}
    with open(fn, 'r') as file:
        csvreader = csv.DictReader(file)
        best_results_per_variation = {}
        for row in csvreader:
            row = trim_dict(row)
            time = float(row['time']) #ms
            total = int(row['args_total'])
            iterations = int(row['args_iterations'])
            powbits = int(row['args_pow'])
            time_per_op = time / (total * iterations * powbits)
            ms_time_per_bilop = (10**9) * time_per_op
            modbits = int(row['args_modbits'])
            function = row['args_mult_class']
            key = (modbits, function)
            if key in best_results_per_variation:
                _, cur = best_results_per_variation[key]
                if ms_time_per_bilop < cur:
                    best_results_per_variation[key] = (row, ms_time_per_bilop)
            else:
                best_results_per_variation[key] = (row, ms_time_per_bilop)
        padding_length = 0
        for k in best_results_per_variation:
            _, name = k
            padding_length = max(len(name), padding_length)
        for k in sorted(best_results_per_variation.keys()):
            modbits, function = k
            row, ms_time_per_bilop = best_results_per_variation[k]
            blocksize = row['args_blocksize']
            blocksper = row['args_blocksper']
            powbits = int(row['args_pow'])
            print(f"{pad_str(function, padding_length)} {modbits}: {ms_time_per_bilop} ms for 1B (x{powbits}); ({blocksize} {blocksper})")
            name = row['args_variation_name'] + row['args_reduction_method_name']
            outx[name] = append_and_return(outx.get(name, []), modbits)
            outy[name] = append_and_return(outy.get(name, []), ms_time_per_bilop)
    return outx, outy

def read_latency_csv(fn, tpi_label=True):
    outx = {}
    outy = {}
    with open(fn, 'r') as file:
        csvreader = csv.DictReader(file)
        best_results_per_variation = {}
        for row in csvreader:
            row = trim_dict(row)
            time = float(row['time']) #ms
            iterations = int(row['args_iterations'])
            powbits = int(row['args_pow'])
            time_per_op = time / (iterations * powbits) # ms
            ns_time_per_op = 1000000 * time_per_op # ns
            modbits = int(row['args_modbits'])
            function = row['args_mult_class']
            tpi = row['args_tpi']
            # function += "tpi=" + str(tpi)
            key = (modbits, function)
            if key in best_results_per_variation:
                _, cur = best_results_per_variation[key]
                if ns_time_per_op < cur:
                    best_results_per_variation[key] = (row, ns_time_per_op)
            else:
                best_results_per_variation[key] = (row, ns_time_per_op)
        padding_length = 0
        for k in best_results_per_variation:
            _, name = k
            padding_length = max(len(name), padding_length)
        for k in sorted(best_results_per_variation.keys()):
            modbits, function = k
            row, ns_time_per_op = best_results_per_variation[k]
            blocksize = row['args_blocksize']
            blocksper = row['args_blocksper']
            powbits = int(row['args_pow'])
            print(f"{pad_str(function, padding_length)} {modbits}: {ns_time_per_op} ns for (x{powbits}); ({blocksize} {blocksper})")
            tpi = row['args_tpi']
            stride = row['args_stride']
            name = row['args_variation_name'] + row['args_reduction_method_name'] 
            if stride:
               name += "stride=" + str(stride)
            if tpi_label:
                name += "tpi=" + str(0 if int(tpi) <= 32 else 1)
            if row['correct'] == "False":
                continue
            # if int(row['local_mem']) > 10:
            #     continue
            # if ms_time_per_bilop > 0.005:
            #     continue
            outx[name] = append_and_return(outx.get(name, []), modbits)
            outy[name] = append_and_return(outy.get(name, []), ns_time_per_op)
    return outx, outy

def append_and_return(list, item):
    list.append(item)
    return list

def extend_and_return(d1, d2, relabel=""):
    outx = {}
    for key, value in d1.items():
        outx[key] = value
    for key, value in d2.items():
        outx[key + relabel] = value
    return outx

def read_dictionary(fn, outx, outy, power):
    with open(fn, "r") as file:
        d = ast.literal_eval(file.read())
        for (modbits, name), time in d.items():
            outx[name] = append_and_return(outx.get(name, []), modbits)
            outy[name] = append_and_return(outy.get(name, []), time / power)

def read_json(fn, outx, outy, name, modbits):
    with open(fn, "r") as f:
        j = json.load(f)
        time = j['mean']['point_estimate']
        outx[name] = append_and_return(outx.get(name, []), modbits)
        outy[name] = append_and_return(outy.get(name, []), time)


def plot_things(x_dict, y_dict, xlabel, ylabel, fig_name, title, y_ticks=None, x_ticks=None):
    print(x_dict, y_dict)
    xmin = 99999999999999
    xmax = 0
    ymin = xmin
    ymax = xmax
    for k in x_dict.keys():
        plt.plot(x_dict[k], y_dict[k], label=k, color=better_colors(k))
        xmin = min(xmin, min(x_dict[k]))
        xmax = max(xmax, max(x_dict[k]))
        ymin = min(ymin, min(y_dict[k]))
        ymax = max(ymax, max(y_dict[k]))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    step = 64
    max_ticks = 10
    if x_ticks is None:
        x_ticks = range(xmin, xmax+1, step)
        while len(x_ticks) > max_ticks:
            step += 64
            x_ticks = range(xmin, xmax+1, step)
    if y_ticks is None:
        y_ticks = np.linspace(ymin, ymax, max_ticks)
    plt.xticks(x_ticks)
    if y_ticks is not None:
        plt.yticks(y_ticks)
    plt.tight_layout()
    plt.legend()
    plt.savefig(f"{fig_name}.png", dpi=300)
    plt.figure()

def vs_baseline_acceleration_plot(baseline, x_dict, y_dict, xlabel, ylabel, title):
    base_pt = x_dict[baseline]
    base_time = y_dict[baseline]
    y_frac = {}
    for k, y_values in y_dict.items():
        for i, x in enumerate(x_dict[k]):
            idx = base_pt.index(x)
            if idx >= 0:
                bt = base_time[idx]
                diff = y_values[i] - bt
                if bt != 0:
                    frac_diff = diff / bt
                else:
                    frac_diff = 0 # some bug
                percent_diff = (frac_diff + 1) * 100
                y_frac[k] = append_and_return(y_frac.get(k, []), percent_diff)
    plt.plot([128, 4096], [25, 25], linestyle='--', color='grey')
    plot_things(x_dict, y_frac, xlabel, ylabel, title, "Percent change compared to Flint baseline", range(0, 201, 25))

def better_colors(lineid):
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    better = {
        "arkworks (CPU rust baseline)": 3,
        "flint (CPU c++ baseline)": 5,
        "AVX512 Opt Mont": 9,
        "CGBN (GPU baseline)": 2,
        "CGBN": 2,
        "GPU Opt Mont with TPI=t": 0,
        "TPI = t": 0,
        "TPI = t shared memory only": 4,
        "GPU Opt Mont with TPI=2t": 1,
        "TPI = 2t": 1,
        "TPI = 2t shared memory only": 6,
    }
    return colors[better[lineid]]

def better_titles(outx, outy):
    nx = {}
    ny = {}
    better = {
        "arkworks": "arkworks (CPU rust baseline)",
        "flint": "flint (CPU c++ baseline)",
        "int": "AVX512 Opt Mont",
        "CGBN": "CGBN (GPU baseline)",
        "TPI = t": "GPU Opt Mont with TPI=t",
        "TPI = 2t": "GPU Opt Mont with TPI=2t"
    }
    for k in outx.keys():
        if k in better:
            nx[better[k]] = outx[k]
            ny[better[k]] = outy[k]
        else:
            nx[k] = outx[k]
            ny[k] = outy[k]
    return nx, ny

def zoomed_plot_data(outx, outy, min, max):
    x_t = {}
    y_t = {}
    for k in outx.keys():
        for i, x in enumerate(outx[k]):
            if (x >= min) and (x < max):
                x_t[k] = append_and_return(x_t.get(k, []), x)
                y_t[k] = append_and_return(y_t.get(k, []), outy[k][i])
    return x_t, y_t

def trim_64(outx, outy):
    for k in outx.keys():
        fx = []
        fy = []
        for (i,x) in enumerate(outx[k]):
            if x % 64 == 0:
                fx.append(x)
                fy.append(outy[k][i])
        outx[k] = fx
        outy[k] = fy

def cpu_plot(power=10000):
    outx = {}
    outy = {}
    read_dictionary("./data/avx_metal.txt", outx, outy, power)
    trim_64(outx, outy)
    outx, outy = better_titles(outx, outy)
    
    vs_baseline_acceleration_plot("flint (CPU c++ baseline)", outx, outy, "Modulus size (bits)", "% vs Flint", "cpu-avx-baseline-cmp")
    rust_plot(outx, outy)
    outx, outy = better_titles(outx, outy)
    
    plot_things(outx, outy, "Modulus size (bits)", "Time (ns)", "cpu-avx-comparison", "Time per modular multiplication")

    plt.plot([300, 300], [5, 44], linestyle='--', color='grey')
    plt.plot([468, 468], [5, 55], linestyle='--', color='grey')
    xs, ys = zoomed_plot_data(outx, outy, 128, 768)
    plot_things(xs, ys, "Modulus size (bits)", "Time (ns)", "cpu-crossover", "Time per modular multiplication (CPU, small moduli)", range(0, 141, 20))
    return outx, outy

def rust_plot(outx, outy):
    # glob_str = "../cpu/baselines/modmul/target/criterion/*/new/estimates.json"
    glob_str = "./data/rust/cpu/baselines/modmul/target/criterion/*/new/estimates.json"
    json_files = glob.glob(glob_str)
    glob_idx = glob_str.find("*") 
    glob_idx2 = len(glob_str) - glob_idx - 1
    glob_idx += 4
    m_f = [(int(fn[glob_idx:-glob_idx2]), fn) for fn in json_files]
    for mf in sorted(m_f):
        modbits, fn = mf
        read_json(fn, outx, outy, "arkworks", modbits)

def get_latest():
    csv_files = glob.glob("data/*.csv")
    csv_files.sort()


def gpu_latency_plot(fn):
    # outx, outy = read_latency_csv(sys.argv[1])
    # fn = 'data/1748014251.8196764.csv'
    outx, outy = read_latency_csv(fn)
    # fn = 'data/1748014568.2274084.csv'
    # outx2, outy2 = read_latency_csv(fn)
    # outx, outy = read_latency_csv("cgbnalltpi.csv")
    # outx2, outy2 = read_latency_csv("latencygpu.csv")
    # outx, outy = extend_and_return(outx, outx2, "w"), extend_and_return(outy, outy2, "w")
    plot_things(outx, outy, "modbits", "time", "latency", "Latency GPU comparison")

if __name__ == "__main__":
    cpu_plot()
    # fn = "data/1748823599.0757022.csv"
    # gpu_latency_plot(fn)