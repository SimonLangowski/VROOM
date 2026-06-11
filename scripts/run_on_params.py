import sys
import time
from basic_benchmarker import destroy_cache

def params_to_args(args):
    # Take command line parameters as overrides
    arglist = sys.argv[1:]
    for (flag_name, value) in args.items():
        flag = "--" + str(flag_name)
        if flag not in arglist:
            arglist.append(flag)
            if value is not None:
                arglist.append(str(value))
    return arglist

def convert_to_csv(results):
    header = []
    result_rows = [[] for _ in range(len(results))]
    for result in results:
        for k in result.keys():
            if k != "args" and k not in header:
                header.append(k)
                for i, r in enumerate(results):
                    v = r.get(k, "")
                    sanitized = str(v).replace(",", ".").replace("\n", ";") # special characters in CSV
                    result_rows[i].append(sanitized)
            elif k == "args":
                arg_dict = vars(result['args'])
                for key in arg_dict.keys():
                    k = "args_" + key
                    if k not in header:
                        header.append(k)
                        for i, r in enumerate(results):
                            a = vars(r.get("args", None))
                            if a is not None:
                                v = a.get(key, "")
                                sanitized = str(v).replace(",", ".").replace("\n", ";") # special characters in CSV
                                result_rows[i].append(sanitized)
        
    csv_strings = [", ".join(header)]
    csv_strings.extend([", ".join(res_values) for res_values in result_rows])
    return "\n".join(csv_strings)

def print_results(results):
    csv_str = convert_to_csv(results)
    print(csv_str) # For output capturing
    fn = "data/" + str(time.time()) + ".csv"
    f = open(fn, "x")
    print(csv_str, file=f)
    return fn

def run_wrapped(args, results, func):
    print(args)
    try:
        result = func(params_to_args(args), True)
        results.append(result)
        return result
    except KeyboardInterrupt as e:
        # Store current results before crashing
        print_results(results)
        raise e
    except RuntimeError as e:
        print(e)
        return None
    except SystemError as e:
        print_results(results)
        raise e
