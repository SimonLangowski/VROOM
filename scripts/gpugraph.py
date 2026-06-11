from create_graph import *


rns_mont_fn = "data/1748888074.0996492.csv" # "data/1748823599.0757022.csv"
rns_mont_fn2 = "data/1748892351.7703345.csv"

all_fn = "data/1748900378.1196628.csv"
# Graph 1: CGBN vs RNSMontgomery (modexp?)
def graph1():
    outx, outy = read_latency_csv(all_fn, False)
    parsed = plot_segments(outx, outy, "CGBN", "CGBN vs RNS Mont")
    all_graph(parsed)
    outx, outy = vs_baseline_acceleration_plot2("CGBNmontmul", outx, outy)
    plot_segments(outx, outy, "CGBNBaseline", "GPU baseline comparison", "% vs CGBN", range(64, 4097, 448), np.arange(0, 125, 12.5))

def all_graph(parsed):
    outx, outy = cpu_plot()
    for k, v in parsed.items():
        lineid, _ = k
        xs, ys = v
        outx[lineid] = xs
        outy[lineid] = ys
    outx, outy = better_titles(outx, outy)
    plot_things(outx, outy, "Modulus size (bits)", "Time (ns)", "ALL", "Time per modular multiplication", range(0, 7001, 700), range(64, 4097, 448))

# Graph 2: TPI and overparallelization
def graph2():
    outx, outy = read_latency_csv(rns_mont_fn)
    plot_segments(outx, outy, "TPI", "Overparallelization")

# Zoomed in, with and without warp shuffle instructions
def graph3():
    outx, outy = read_latency_csv(rns_mont_fn, False)
    outx3, outy3 = read_latency_csv(all_fn, False)
    outx2, outy2 = read_latency_csv(rns_mont_fn2, False)
    outx, outy = extend_and_return(outx, outx2, "shared"), extend_and_return(outy, outy2, "shared")
    cgbn = "CGBNmontmul"
    outx[cgbn] = outx3[cgbn]
    outy[cgbn] = outy3[cgbn]
    plot_segments(outx, outy, "TPIWS", "Time per modular multiplication (GPU, small moduli)", "Time (ns)", range(64, 1025, 96), range(0, 721, 72))

# Perhaps alternate solid and dotted for different TPI regions?
# Label with text
# Different color for warp shuffle, and stride

def plot_segments(outx, outy, fig_name, title, ylabel="Time (ns)", x_ticks=None, y_ticks=None):
    plt.title(title)
    plt.xlabel("Modulus size (bits)")
    plt.ylabel(ylabel)
    xmin = 99999999999999
    xmax = 0
    ymin = xmin
    ymax = xmax
    with_connected_tpi = {}
    linenames = {
        "CGBNmontmul" : "CGBN",
        "MontRNSReduceVecMont32stride=1" : "TPI = t",
        "MontRNSReduceVecMont32stride=2" : "TPI = 2t",
        "MontRNSReduceVecMont32stride=4" : "TPI = 4t",
        "MontRNSReduceVecMont32stride=1shared" : "TPI = t",
        "MontRNSReduceVecMont32stride=2shared" : "TPI = 2t",
        "MontRNSReduceVecMont32stride=4shared" : "TPI = 4t",
    }
    for k in outx.keys():
        x = outx[k]
        y = outy[k]
        parts = k.split("tpi=")
        lineid = linenames[parts[0]]
        if len(parts) > 1:
            tpi = parts[1]
        else:
            tpi = 0
        if "4" in lineid:
            continue
        if fig_name == "TPIWS" and "CGBN" not in lineid:
            if "shared" not in parts[0]:
                # lineid += " warp shuffle"
                # if tpi == '1':
                #     continue
                lineid = lineid
            else:
                lineid += " shared memory only"
                # tpi = tpi[0]
        if "TPI = t" and fig_name == "TPI":
            print([x[i] if y[i] > 0.002 else 0 for i in range(len(y))])
        kk = (lineid, tpi)
        xs, ys = with_connected_tpi.get(kk, ([], []))
        assert(len(x) == len(y))
        for i in range(len(x)):
            if fig_name != "TPIWS" or (x[i] < 1024 and y[i] < 720):
                if x[i] < 4032 or lineid == "CGBN":
                    xs.append(x[i])
                    ys.append(y[i])
        if len(xs) > 0:
            with_connected_tpi[kk] = (xs, ys)
    if fig_name == "TPI":
        styles = {
            "TPI = t": "solid",
            "TPI = 2t": "dashed",
            # "TPI = 4t": "dotted",
        }
    elif fig_name == "TPIWS":
        styles = {
            "TPI = t": "solid",
            "TPI = 2t": "dashed",
            "TPI = t shared memory only": "solid",
            "TPI = 2t shared memory only": "dashed",
            "CGBN": "solid",
        }
    elif fig_name == "CGBN" or fig_name == "CGBNBaseline":
        styles = {
            "TPI = t": "solid",
            "TPI = 2t": "dashed",
            "CGBN": "solid",
        }
    legend = sorted(list(styles.keys()))
    if fig_name == "TPIWS":
        legend.reverse()
        # # warp shuffle is not used beyond tpi=32 even if set in code
        # for i in range(0, 4, 2):
        #     x, y = with_connected_tpi[(legend[i], 0)]
        #     x2, y2 = with_connected_tpi[(legend[i+1], 0)]
        #     if len(x2) > len(x):
        #         ws_len = len(x)
        #         x2n = x2[:ws_len]
        #         x.extend(x2[ws_len:])
        #         y2n = y2[:ws_len]
        #         y.extend(y2[ws_len:])
        #         with_connected_tpi[(legend[i], 0)] = (x, y)
        #         with_connected_tpi[(legend[i+1], 0)] = (x2, y2)
    if fig_name == "TPIWS" or fig_name == "CGBN" or fig_name == "CGBNBaseline":
        wo = {}
        # Manual insert so legend in correct ordering
        for l in legend:
            kk = (l, '8')
            if kk not in with_connected_tpi:
                kk = (l, '16')
            if kk not in with_connected_tpi:
                kk = (l, 0)
            if kk not in with_connected_tpi:
                kk = (l, '0shared')
            if kk not in with_connected_tpi:
                kk = (l, '1shared')
            if kk not in with_connected_tpi:
                kk = (l, '0')
            wo[kk] = with_connected_tpi[kk]
        for k, v in with_connected_tpi.items():
            if k not in wo:
                wo[k] = v
        with_connected_tpi = wo

    tpicolors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    for k, v in with_connected_tpi.items():
        lineid, tpi = k
        if (int(tpi) > 256):
            continue
        id = list(styles.keys()).index(lineid)
        xs, ys = v
        # style = "dashed" if int(tpi).bit_length() % 2 == 1 else "solid"
        # color=colors[lineid]
        # color = tpicolors[tpi]
        # style = styles[lineid]
        color =  better_colors(lineid) #tpicolors[id]
        plt.plot(xs, ys, label=lineid, color=color)
        # Plot segment with dash/solid?
        
        # Label midpoint of segment
        xmid = (min(xs) + max(xs)) / 2
        ymid = (min(ys) + max(ys)) / 2
        # may need to offset upwards
       
        if fig_name == "TPI":
            spacing = 0.00015
        elif fig_name == "TPIWS":
            spacing = 0.00005
        elif fig_name == "CGBN":
            spacing = 0.00015
        elif fig_name == "CGBNBaseline":
            spacing = 0
        ymid += -(2*id - 1) * spacing # * (-1 if int(tpi) >= 128 else 1)
        if id < 2 and fig_name == "TPI":
            plt.text(xmid, ymid, tpi)

        xmin = min(xmin, min(xs))
        xmax = max(xmax, max(xs))
        ymin = min(ymin, min(ys))
        ymax = max(ymax, max(ys))
    if fig_name == "CGBNBaseline":
        plt.plot([64, 4096], [25, 25], linestyle='--', color='grey')

        plt.plot([1980, 1980], [12.5, 75], linestyle='--', color='grey')
        plt.text(2700, 15, "TPI=256", color=better_colors("TPI = 2t"))
        plt.text(2700, 62.5, "TPI=128", color=better_colors("TPI = t"))
        plt.plot([960, 960], [12.5, 75], linestyle='--', color='grey')
        plt.text(1200, 15, "TPI=128", color=better_colors("TPI = 2t"))
        plt.text(1200, 62.5, "TPI=64", color=better_colors("TPI = t"))
        plt.plot([480, 480], [12.5, 75], linestyle='--', color='grey')
        plt.text(500, 15, "TPI=64", color=better_colors("TPI = 2t"))
        plt.text(500, 62.5, "TPI=32", color=better_colors("TPI = t"))
        plt.text(00, 15, "TPI=32", color=better_colors("TPI = 2t"))
    if fig_name == "TPIWS":
        plt.plot([442, 442], [72, 576], linestyle='--', color='grey')
        plt.text(700, 200, "TPI=64", color=better_colors("TPI = 2t"))
        plt.text(700, 650, "TPI=32", color=better_colors("TPI = t"))
        plt.text(200, 200, "TPI=32", color=better_colors("TPI = 2t"))
    step = 64
    max_ticks = 10
    if x_ticks is None:
        x_ticks = range(xmin, xmax+1, step)
        while len(x_ticks) > max_ticks:
            step += 64
            x_ticks = range(xmin, xmax+1, step)
    if y_ticks is None:
        y_ticks = np.linspace(0, ymax, max_ticks)
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)
    plt.tight_layout()
    if fig_name == "TPIWS":
        plt.legend(legend, loc='upper left')
    else:
        plt.legend(legend)
    plt.savefig(f"{fig_name}.png", dpi=300)
    plt.figure()
    return with_connected_tpi


def vs_baseline_acceleration_plot2(baseline, x_dict, y_dict):
    base_pt = x_dict[baseline]
    base_time = y_dict[baseline]
    y_frac = {}
    for k, y_values in y_dict.items():
        for i, x in enumerate(x_dict[k]):
            # Because CGBN has so few points I need to do a weighted average
            idx = 0
            while base_pt[idx] < x:
                idx += 1
            if idx == 0:
                idx = 1
            x0 = base_pt[idx - 1]
            x1 = base_pt[idx]
            y0 = base_time[idx - 1]
            y1 = base_time[idx]
            weighted_average = y0 + (y1 - y0) * ((x - x0)/(x1 - x0))
            diff = y_values[i] - weighted_average
            if weighted_average != 0:
                frac_diff = diff / weighted_average
            else:
                frac_diff = 0 # some bug
            percent_diff = (frac_diff + 1) * 100
            y_frac[k] = append_and_return(y_frac.get(k, []), percent_diff)
    return x_dict, y_frac

if __name__ == "__main__":
    graph1()
    graph2()
    graph3()