
# Could unconditionally add 10 * (2**52 - t) + ((30 * p) % (2*52 - t)) for correction, then full reduce.
# Could multiply by (-)2/3 to make linear term into quadratic term, then full reduce, to avoid blowup

q = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab


def parse_polys(poly_str):
    polys = []
    poly_lines = poly_str.split('\n')
    for line in poly_lines:
        if 'Polynomial' in line:
            coeff_map = {}
            parts = [x.strip() for x in line.split(":")[1].strip().split(" ")]
            assert(len(parts) % 2 == 0)
            for i in range(0, len(parts), 2):
                coeff = int(parts[i])
                term = parts[i+1]
                if "y" in term and "x" in term:
                    # x*y term
                    x = int(term.split("x")[1].split("y")[0])
                    y = int(term.split("y")[1])
                    coeff_map[(x, ("y", y))] = coeff
                elif term.startswith("y") and term.count("y") == 2:
                    # y*y term
                    parts_y = term.split("y")
                    y1 = int(parts_y[1])
                    y2 = int(parts_y[2])
                    coeff_map[(("y", y1), ("y", y2))] = coeff
                elif len(term.split("x")) == 3:
                    # x*x term
                    x1 = int(term.split("x")[1])
                    x2 = int(term.split("x")[2])
                    coeff_map[(x1, ("x", x2))] = coeff
                elif term.startswith("x"):
                    # linear x term
                    x = int(term.split("x")[1])
                    coeff_map[(x, -1)] = coeff
                elif term.startswith("y"):
                    # linear y term
                    y = int(term.split("y")[1])
                    coeff_map[(("y", y), -1)] = coeff
            polys.append(coeff_map)
    return polys

def group_pos_neg_terms(poly0, poly1):
    p0_unique = {}
    p1_unique = {}
    pos_terms = {}
    neg_terms = {}
    for term in poly0.keys():
        if term in poly1:
            if poly0[term] == poly1[term]:
                pos_terms[term] = poly0[term]
            elif poly0[term] == -poly1[term]:
                neg_terms[term] = -poly0[term]
            else:
                raise ValueError(f"Coefficients do not match for term {term}")
        else:
            p0_unique[term] = poly0[term]
    for term in poly1.keys():
        if term not in poly0:
            p1_unique[term] = poly1[term]
    return p0_unique, p1_unique, pos_terms, neg_terms

def relinearize(terms):
    new_terms = {}
    for term in terms.keys():
        if term[1] == -1:
            if terms[term] == 2:
                new_terms[(term[0], 0.67)] = 3
            elif terms[term] == -2:
                new_terms[(term[0], 0.67)] = -3
            elif terms[term] == 1:
                new_terms[(term[0], 0.33)] = 3
            elif terms[term] == -1:
                new_terms[(term[0], 0.33)] = -3
            else:
                # Preserve linear terms that aren't 1, -1, 2, or -2
                new_terms[term] = terms[term]
        else:
            new_terms[term] = terms[term]
    return new_terms

def factor_common_coeffcient(terms, factor):
    terms_with_factor = {}
    terms_without_factor = {}
    for term in terms.keys():
        if terms[term] % factor == 0:
            terms_with_factor[term] = terms[term] // factor
        else:
            terms_without_factor[term] = terms[term]
    return terms_with_factor, terms_without_factor

def apply_strategies(pos_terms, neg_terms, p0_unique, p1_unique):
    # 1 negative term, make it unique
    # 2 multiplies, versus 1 multiply, 1 add, 1 subtract
    if len(neg_terms) == 1:
        neg_term = list(neg_terms.keys())[0]
        p0_unique[neg_term] = -neg_terms[neg_term]
        p1_unique[neg_term] = neg_terms[neg_term]
        neg_terms = {}
    # Change the linear terms to quadratic term
    pos_terms = relinearize(pos_terms)
    neg_terms = relinearize(neg_terms)
    p0_unique = relinearize(p0_unique)
    p1_unique = relinearize(p1_unique)
    
    common_factor = 1
    pos_f3, neg_f3, p0_f3, p1_f3 = factor_common_coeffcient(pos_terms, 3), factor_common_coeffcient(neg_terms, 3), factor_common_coeffcient(p0_unique, 3), factor_common_coeffcient(p1_unique, 3)
    if len(pos_f3[1]) == 0 and len(neg_f3[1]) == 0 and len(p0_f3[1]) == 0 and len(p1_f3[1]) == 0:
        pos_terms = pos_f3[0]
        neg_terms = neg_f3[0]
        p0_unique = p0_f3[0]
        p1_unique = p1_f3[0]
        common_factor *= 3

    pos_f2, neg_f2, p0_f2, p1_f2 = factor_common_coeffcient(pos_terms, 2), factor_common_coeffcient(neg_terms, 2), factor_common_coeffcient(p0_unique, 2), factor_common_coeffcient(p1_unique, 2)
    if len(pos_f2[1]) == 0 and len(neg_f2[1]) == 0 and len(p0_f2[1]) == 0 and len(p1_f2[1]) == 0:
        pos_terms = pos_f2[0]
        neg_terms = neg_f2[0]
        p0_unique = p0_f2[0]
        p1_unique = p1_f2[0]
        common_factor *= 2

    return pos_terms, neg_terms, p0_unique, p1_unique, common_factor

def print_poly(terms):
    string = ""
    if len(terms) == 0:
        return "{}"
    for i, term in enumerate(terms.keys()):
        if term[1] == -1:
            # Linear term
            if isinstance(term[0], tuple):
                # y variable
                y_idx = term[0][1]
                if terms[term] == 1:
                    string += f"y[{y_idx}]"
                elif terms[term] == -1:
                    string += f"minus_y{y_idx}"
                elif terms[term] < 0:
                    string += f"{-terms[term]} * minus_y{y_idx}"
                else:
                    string += f"{terms[term]} * y[{y_idx}]"
            else:
                # x variable
                if terms[term] == 1:
                    string += f"x[{term[0]}]"
                elif terms[term] == -1:
                    string += f"minus_x{term[0]}"
                elif terms[term] < 0:
                    string += f"{-terms[term]} * minus_x{term[0]}"
                else:
                    string += f"{terms[term]} * x[{term[0]}]"
        elif isinstance(term[1], tuple):
            # Second variable is a tuple (type, index)
            var1_type, var1_idx = term[1]
            if isinstance(term[0], tuple):
                # y * y term
                var0_type, var0_idx = term[0]
                if terms[term] == 1:
                    string += f"y[{var0_idx}] * y[{var1_idx}]"
                elif terms[term] == -1:
                    # Pull out negate: y * minus_y (one negate beforehand)
                    string += f"y[{var0_idx}] * minus_y{var1_idx}"
                elif terms[term] < 0:
                    string += f"{-terms[term]} * y[{var0_idx}] * minus_y{var1_idx}"
                else:
                    string += f"{terms[term]} * y[{var0_idx}] * y[{var1_idx}]"
            else:
                # x * y or x * x term
                if var1_type == "y":
                    # x * y term
                    if terms[term] == 1:
                        string += f"x[{term[0]}] * y[{var1_idx}]"
                    elif terms[term] == -1:
                        # Negative coefficient: prefer negating y term
                        string += f"x[{term[0]}] * minus_y{var1_idx}"
                    elif terms[term] < 0:
                        # Negative coefficient: prefer negating y term
                        string += f"{-terms[term]} * x[{term[0]}] * minus_y{var1_idx}"
                    else:
                        string += f"{terms[term]} * x[{term[0]}] * y[{var1_idx}]"
                else:
                    # x * x term: pull out negate as x * minus_y (one negate beforehand; minus_y -> minus_x for sqr/inverse)
                    if terms[term] == 1:
                        string += f"x[{term[0]}] * x[{var1_idx}]"
                    elif terms[term] == -1:
                        string += f"x[{term[0]}] * minus_y{var1_idx}"
                    elif terms[term] < 0:
                        string += f"{-terms[term]} * x[{term[0]}] * minus_y{var1_idx}"
                    else:
                        string += f"{terms[term]} * x[{term[0]}] * x[{var1_idx}]"
        else:
            # Both are integers - x * x term (old format): pull out negate as x * minus_y
            if terms[term] == 1:
                string += f"x[{term[0]}] * x[{term[1]}]"
            elif terms[term] == -1:
                string += f"x[{term[0]}] * minus_y{term[1]}"
            elif terms[term] < 0:
                string += f"{-terms[term]} * x[{term[0]}] * minus_y{term[1]}"
            else:
                string += f"{terms[term]} * x[{term[0]}] * x[{term[1]}]"
        if i < len(terms) - 1:
            string += " + "
    return string

def print_split_poly(terms, factor=2):
    terms_with_factor, terms_without_factor = factor_common_coeffcient(terms, factor)
    if len(terms_with_factor) == 0:
        # No terms have factor
        return f"{print_poly(terms_without_factor)}"
    elif len(terms_without_factor) == 0:
        # All terms have factor
        return f"({print_poly(terms_with_factor)}) * {factor}"
    else:
        return f"(({print_poly(terms_with_factor)}) * {factor}) + {print_poly(terms_without_factor)}"

def print_negatives(string, poly_len):
    negatives_needed = []
    for i in range(poly_len):
        if f"minus_y{i}" in string:
            negatives_needed.append(i)
    negatives_string = ""
    for i in range(poly_len):
        if i in negatives_needed:
            negatives_string += f"\tauto minus_y{i} = ring.negate(y[{i}]);\n"
    return negatives_string + string


def print_function(name,poly_str):
    string = ""
    polys = parse_polys(poly_str)
    for i in range(0, len(polys), 2):
        p0 = polys[i]
        if i + 1 >= len(polys):
            # Odd number of polynomials - handle last one separately
            p0_unique = p0
            p1_unique = {}
            pos_terms = {}
            neg_terms = {}
            pos_terms, neg_terms, p0_unique, p1_unique, common_factor = apply_strategies(pos_terms, neg_terms, p0_unique, p1_unique)
            stringi = f"\tauto poly_{i}_unique = "
            if common_factor != 1:
                stringi += f"("
            if len(p0_unique) > 0:
                stringi += ""
                stringi += f"{print_split_poly(p0_unique)}"
            if common_factor != 1:
                stringi += f") * {common_factor}"
            stringi += f";"
            string += stringi + "\n"
            string += f"\treturn ring.template batch_reduce_expand(" + ",".join([f"poly_{j}_unique" for j in range(len(polys))]) + ");\n"
            string += "}"
            if "0.67" in string or "0.33" in string:
                string = string.replace("x[0.67]", "two_thirds")
                string = string.replace("y[0.67]", "two_thirds")
                string = string.replace("minus_x0.67", "neg_two_thirds")
                string = string.replace("minus_y0.67", "neg_two_thirds")
                string = string.replace("x[0.33]", "one_third")
                string = string.replace("y[0.33]", "one_third")
                string = string.replace("minus_x0.33", "neg_one_third")
                string = string.replace("minus_y0.33", "neg_one_third")
                two_thirds = (2 * pow(3, -1, q)) % q
                neg_two_thirds = (q - two_thirds) % q
                c = f"\tauto two_thirds = ring.wrap(ring.template from_hex(\"{hex(two_thirds)}\"));\n"
                c += f"\tauto neg_two_thirds = ring.wrap(ring.template from_hex(\"{hex(neg_two_thirds)}\"));\n"
                string = c + string
            string = print_negatives(string, 12)
            if "sqr" in name or "inverse" in name:
                # For squaring and inversion functions, all variables are x
                string = string.replace("y[", "x[")
                string = string.replace("minus_y", "minus_x")
            string = f"/* {poly_str} */\n" +  name + " {\n" + string
            print(string)
            return string
        p1 = polys[i+1]
        p0_unique, p1_unique, pos_terms, neg_terms = group_pos_neg_terms(p0, p1)
        string += f"\t// {i} {i+1}\n"
        pos_terms, neg_terms, p0_unique, p1_unique, common_factor = apply_strategies(pos_terms, neg_terms, p0_unique, p1_unique)
        if len(pos_terms) > 0:
            string += f"\tauto poly_{i}_pos_shared = offset + {print_split_poly(pos_terms)};\n"
        if len(neg_terms) > 0:
            string += f"\tauto poly_{i}_neg_shared = {print_split_poly(neg_terms)};\n"
        stringi = f"\tauto poly_{i}_unique = "
        stringip1 = f"\tauto poly_{i+1}_unique = "
        if common_factor != 1:
            stringi += f"("
            stringip1 += f"("
        started = False
        if len(pos_terms) > 0:
            stringi += f"poly_{i}_pos_shared"
            stringip1 += f"poly_{i}_pos_shared"
            started = True
        if len(neg_terms) > 0:
            stringi += f" - poly_{i}_neg_shared"
            stringip1 += f" + poly_{i}_neg_shared"
            started = True
        if len(p0_unique) > 0:
            if started:
                stringi += f" + "
            else:
                stringi += ""
            stringi += f"{print_split_poly(p0_unique)}"
        if len(p1_unique) > 0:
            if started:
                stringip1 += f" + "
            else:
                stringip1 += ""
            stringip1 += f"{print_split_poly(p1_unique)}"
        if common_factor != 1:
            stringi += f") * {common_factor}"
            stringip1 += f") * {common_factor}"
        stringi += f";"
        stringip1 += f";"
        string += stringi + "\n" + stringip1 + "\n"
    string += f"\treturn ring.template batch_reduce_expand(" + ",".join([f"poly_{i}_unique" for i in range(len(polys))]) + ");\n"
    string += "}"
    if "0.67" in string or "0.33" in string:
        string = string.replace("x[0.67]", "two_thirds")
        string = string.replace("y[0.67]", "two_thirds")
        string = string.replace("minus_x0.67", "neg_two_thirds")
        string = string.replace("minus_y0.67", "neg_two_thirds")
        string = string.replace("x[0.33]", "one_third")
        string = string.replace("y[0.33]", "one_third")
        string = string.replace("minus_x0.33", "neg_one_third")
        string = string.replace("minus_y0.33", "neg_one_third")
        two_thirds = (2 * pow(3, -1, q)) % q
        neg_two_thirds = (q - two_thirds) % q
        c = f"\tauto two_thirds = ring.wrap(ring.template from_hex(\"{hex(two_thirds)}\"));\n"
        c += f"\tauto neg_two_thirds = ring.wrap(ring.template from_hex(\"{hex(neg_two_thirds)}\"));\n"
        if "0.33" in string:
            one_third = pow(3, -1, q) % q
            neg_one_third = (q - one_third) % q
            c += f"\tauto one_third = ring.wrap(ring.template from_hex(\"{hex(one_third)}\"));\n"
            c += f"\tauto neg_one_third = ring.wrap(ring.template from_hex(\"{hex(neg_one_third)}\"));\n"
        string = c + string
        string = string
    string = print_negatives(string, 12)
    if "sqr" in name or "inverse" in name:
        string = string.replace("y[", "x[")
        string = string.replace("minus_y", "minus_x")

    string = f"/* {poly_str} */\n" +  name + " {\n" + string
    print(string)
    return string

fp12 = '''
Polynomial 0(22) : 1 x0y0 -1 x1y1 +1 x2y4 -1 x2y5 -1 x3y4 -1 x3y5 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3 +1 x6y10 -1 x6y11 -1 x7y10 -1 x7y11 +1 x8y8 -1 x8y9 -1 x9y8 -1 x9y9 +1 x10y6 -1 x10y7 -1 x11y6 -1 x11y7

Polynomial 1(22) : 1 x0y1 +1 x1y0 +1 x2y4 +1 x2y5 +1 x3y4 -1 x3y5 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3 +1 x6y10 +1 x6y11 +1 x7y10 -1 x7y11 +1 x8y8 +1 x8y9 +1 x9y8 -1 x9y9 +1 x10y6 +1 x10y7 +1 x11y6 -1 x11y7

Polynomial 2(18) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5 +1 x6y6 -1 x7y7 +1 x8y10 -1 x8y11 -1 x9y10 -1 x9y11 +1 x10y8 -1 x10y9 -1 x11y8 -1 x11y9

Polynomial 3(18) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5 +1 x6y7 +1 x7y6 +1 x8y10 +1 x8y11 +1 x9y10 -1 x9y11 +1 x10y8 +1 x10y9 +1 x11y8 -1 x11y9

Polynomial 4(14) : 1 x0y4 -1 x1y5 +1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1 +1 x6y8 -1 x7y9 +1 x8y6 -1 x9y7 +1 x10y10 -1 x10y11 -1 x11y10 -1 x11y11

Polynomial 5(14) : 1 x0y5 +1 x1y4 +1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0 +1 x6y9 +1 x7y8 +1 x8y7 +1 x9y6 +1 x10y10 +1 x10y11 +1 x11y10 -1 x11y11

Polynomial 6(20) : 1 x0y6 -1 x1y7 +1 x2y10 -1 x2y11 -1 x3y10 -1 x3y11 +1 x4y8 -1 x4y9 -1 x5y8 -1 x5y9 +1 x6y0 -1 x7y1 +1 x8y4 -1 x8y5 -1 x9y4 -1 x9y5 +1 x10y2 -1 x10y3 -1 x11y2 -1 x11y3

Polynomial 7(20) : 1 x0y7 +1 x1y6 +1 x2y10 +1 x2y11 +1 x3y10 -1 x3y11 +1 x4y8 +1 x4y9 +1 x5y8 -1 x5y9 +1 x6y1 +1 x7y0 +1 x8y4 +1 x8y5 +1 x9y4 -1 x9y5 +1 x10y2 +1 x10y3 +1 x11y2 -1 x11y3

Polynomial 8(16) : 1 x0y8 -1 x1y9 +1 x2y6 -1 x3y7 +1 x4y10 -1 x4y11 -1 x5y10 -1 x5y11 +1 x6y2 -1 x7y3 +1 x8y0 -1 x9y1 +1 x10y4 -1 x10y5 -1 x11y4 -1 x11y5

Polynomial 9(16) : 1 x0y9 +1 x1y8 +1 x2y7 +1 x3y6 +1 x4y10 +1 x4y11 +1 x5y10 -1 x5y11 +1 x6y3 +1 x7y2 +1 x8y1 +1 x9y0 +1 x10y4 +1 x10y5 +1 x11y4 -1 x11y5

Polynomial 10(12) : 1 x0y10 -1 x1y11 +1 x2y8 -1 x3y9 +1 x4y6 -1 x5y7 +1 x6y4 -1 x7y5 +1 x8y2 -1 x9y3 +1 x10y0 -1 x11y1

Polynomial 11(12) : 1 x0y11 +1 x1y10 +1 x2y9 +1 x3y8 +1 x4y7 +1 x5y6 +1 x6y5 +1 x7y4 +1 x8y3 +1 x9y2 +1 x10y1 +1 x11y0
'''

mulx0y0z00 = '''
Polynomial 0(10) : 1 x0y0 -1 x1y1 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3 +1 x8y4 -1 x8y5 -1 x9y4 -1 x9y5

Polynomial 1(10) : 1 x0y1 +1 x1y0 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3 +1 x8y4 +1 x8y5 +1 x9y4 -1 x9y5

Polynomial 2(8) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x10y4 -1 x10y5 -1 x11y4 -1 x11y5

Polynomial 3(8) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x10y4 +1 x10y5 +1 x11y4 -1 x11y5

Polynomial 4(6) : 1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1 +1 x6y4 -1 x7y5

Polynomial 5(6) : 1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0 +1 x6y5 +1 x7y4

Polynomial 6(10) : 1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5 +1 x6y0 -1 x7y1 +1 x10y2 -1 x10y3 -1 x11y2 -1 x11y3

Polynomial 7(10) : 1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5 +1 x6y1 +1 x7y0 +1 x10y2 +1 x10y3 +1 x11y2 -1 x11y3

Polynomial 8(6) : 1 x0y4 -1 x1y5 +1 x6y2 -1 x7y3 +1 x8y0 -1 x9y1

Polynomial 9(6) : 1 x0y5 +1 x1y4 +1 x6y3 +1 x7y2 +1 x8y1 +1 x9y0

Polynomial 10(6) : 1 x2y4 -1 x3y5 +1 x8y2 -1 x9y3 +1 x10y0 -1 x11y1

Polynomial 11(6) : 1 x2y5 +1 x3y4 +1 x8y3 +1 x9y2 +1 x10y1 +1 x11y0
'''

inv_fp12_to_fp6 = '''
Polynomial 0(13) : 1 x0x0 -1 x1x1 -1 x8x8 +1 x9x9 +2 x2x4 -2 x2x5 -2 x3x4 -2 x3x5 -2 x6x10 +2 x6x11 +2 x7x10 +2 x7x11 +2 x8x9

Polynomial 1(12) : -1 x8x8 +1 x9x9 +2 x0x1 +2 x2x4 +2 x2x5 +2 x3x4 -2 x3x5 -2 x6x10 -2 x6x11 -2 x7x10 +2 x7x11 -2 x8x9

Polynomial 2(11) : 1 x4x4 -1 x5x5 -1 x6x6 +1 x7x7 +2 x0x2 -2 x1x3 -2 x4x5 -2 x8x10 +2 x8x11 +2 x9x10 +2 x9x11

Polynomial 3(10) : 1 x4x4 -1 x5x5 +2 x0x3 +2 x1x2 +2 x4x5 -2 x6x7 -2 x8x10 -2 x8x11 -2 x9x10 +2 x9x11

Polynomial 4(9) : 1 x2x2 -1 x3x3 -1 x10x10 +1 x11x11 +2 x0x4 -2 x1x5 -2 x6x8 +2 x7x9 +2 x10x11

Polynomial 5(8) : -1 x10x10 +1 x11x11 +2 x0x5 +2 x1x4 +2 x2x3 -2 x6x9 -2 x7x8 -2 x10x11
'''
inv_fp6_to_fp12 = '''
Polynomial 0(10) : 1 x0y0 -1 x1y1 +1 x2y4 -1 x2y5 -1 x3y4 -1 x3y5 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3

Polynomial 1(10) : 1 x0y1 +1 x1y0 +1 x2y4 +1 x2y5 +1 x3y4 -1 x3y5 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3

Polynomial 2(8) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5

Polynomial 3(8) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5

Polynomial 4(6) : 1 x0y4 -1 x1y5 +1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1

Polynomial 5(6) : 1 x0y5 +1 x1y4 +1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0

Polynomial 6(10) : -1 x6y0 +1 x7y1 -1 x8y4 +1 x8y5 +1 x9y4 +1 x9y5 -1 x10y2 +1 x10y3 +1 x11y2 +1 x11y3

Polynomial 7(10) : -1 x6y1 -1 x7y0 -1 x8y4 -1 x8y5 -1 x9y4 +1 x9y5 -1 x10y2 -1 x10y3 -1 x11y2 +1 x11y3

Polynomial 8(8) : -1 x6y2 +1 x7y3 -1 x8y0 +1 x9y1 -1 x10y4 +1 x10y5 +1 x11y4 +1 x11y5

Polynomial 9(8) : -1 x6y3 -1 x7y2 -1 x8y1 -1 x9y0 -1 x10y4 -1 x10y5 -1 x11y4 +1 x11y5

Polynomial 10(6) : -1 x6y4 +1 x7y5 -1 x8y2 +1 x9y3 -1 x10y0 +1 x11y1

Polynomial 11(6) : -1 x6y5 -1 x7y4 -1 x8y3 -1 x9y2 -1 x10y1 -1 x11y0
'''
inv_fp6_to_fp2_c0_c1_c2 = '''
Polynomial 0(6) : 1 x0x0 -1 x1x1 -1 x2x4 +1 x2x5 +1 x3x4 +1 x3x5

Polynomial 1(5) : 2 x0x1 -1 x2x4 -1 x2x5 -1 x3x4 +1 x3x5

Polynomial 2(5) : 1 x4x4 -1 x5x5 -1 x0x2 +1 x1x3 -2 x4x5

Polynomial 3(5) : 1 x4x4 -1 x5x5 -1 x0x3 -1 x1x2 +2 x4x5

Polynomial 4(4) : 1 x2x2 -1 x3x3 -1 x0x4 +1 x1x5

Polynomial 5(3) : -1 x0x5 -1 x1x4 +2 x2x3
'''
inv_fp6_to_fp2_ret = '''
Polynomial 0(10) : 1 x0x6 -1 x1x7 +1 x2x10 -1 x2x11 -1 x3x10 -1 x3x11 +1 x4x8 -1 x4x9 -1 x5x8 -1 x5x9

Polynomial 1(10) : 1 x0x7 +1 x1x6 +1 x2x10 +1 x2x11 +1 x3x10 -1 x3x11 +1 x4x8 +1 x4x9 +1 x5x8 -1 x5x9
'''

inv_fp2_to_fp6 = '''
Polynomial 0(2) : 1 x0x2 -1 x1x3

Polynomial 1(2) : 1 x0x3 +1 x1x2

Polynomial 2(2) : 1 x0x4 -1 x1x5

Polynomial 3(2) : 1 x0x5 +1 x1x4

Polynomial 4(2) : 1 x0y0 -1 x1y1

Polynomial 5(2) : 1 x0y1 +1 x1y0
'''

inv_fp2_to_fp = '''
Polynomial 0(2) : 1 x0x0 +1 x1x1
'''
inv_fp_to_fp2 = '''
Polynomial 0(1) : 1 x0y0

Polynomial 1(1) : -1 x1y0
''' 

fp12sqr = '''
Polynomial 0(13) : 1 x0x0 -1 x1x1 +1 x8x8 -1 x9x9 +2 x2x4 -2 x2x5 -2 x3x4 -2 x3x5 +2 x6x10 -2 x6x11 -2 x7x10 -2 x7x11 -2 x8x9

Polynomial 1(12) : 1 x8x8 -1 x9x9 +2 x0x1 +2 x2x4 +2 x2x5 +2 x3x4 -2 x3x5 +2 x6x10 +2 x6x11 +2 x7x10 -2 x7x11 +2 x8x9

Polynomial 2(11) : 1 x4x4 -1 x5x5 +1 x6x6 -1 x7x7 +2 x0x2 -2 x1x3 -2 x4x5 +2 x8x10 -2 x8x11 -2 x9x10 -2 x9x11

Polynomial 3(10) : 1 x4x4 -1 x5x5 +2 x0x3 +2 x1x2 +2 x4x5 +2 x6x7 +2 x8x10 +2 x8x11 +2 x9x10 -2 x9x11

Polynomial 4(9) : 1 x2x2 -1 x3x3 +1 x10x10 -1 x11x11 +2 x0x4 -2 x1x5 +2 x6x8 -2 x7x9 -2 x10x11

Polynomial 5(8) : 1 x10x10 -1 x11x11 +2 x0x5 +2 x1x4 +2 x2x3 +2 x6x9 +2 x7x8 +2 x10x11

Polynomial 6(10) : 2 x0x6 -2 x1x7 +2 x2x10 -2 x2x11 -2 x3x10 -2 x3x11 +2 x4x8 -2 x4x9 -2 x5x8 -2 x5x9

Polynomial 7(10) : 2 x0x7 +2 x1x6 +2 x2x10 +2 x2x11 +2 x3x10 -2 x3x11 +2 x4x8 +2 x4x9 +2 x5x8 -2 x5x9

Polynomial 8(8) : 2 x0x8 -2 x1x9 +2 x2x6 -2 x3x7 +2 x4x10 -2 x4x11 -2 x5x10 -2 x5x11

Polynomial 9(8) : 2 x0x9 +2 x1x8 +2 x2x7 +2 x3x6 +2 x4x10 +2 x4x11 +2 x5x10 -2 x5x11

Polynomial 10(6) : 2 x0x10 -2 x1x11 +2 x2x8 -2 x3x9 +2 x4x6 -2 x5x7

Polynomial 11(6) : 2 x0x11 +2 x1x10 +2 x2x9 +2 x3x8 +2 x4x7 +2 x5x6
'''

cyclotomic_sqr_fp12 = '''
Polynomial 0(6) : 3 x0x0 -2 x0 -3 x1x1 +3 x8x8 -3 x9x9 -6 x8x9

Polynomial 1(5) : -2 x1 +3 x8x8 -3 x9x9 +6 x0x1 +6 x8x9

Polynomial 2(6) : -2 x2 +3 x4x4 -3 x5x5 +3 x6x6 -3 x7x7 -6 x4x5

Polynomial 3(5) : -2 x3 +3 x4x4 -3 x5x5 +6 x4x5 +6 x6x7

Polynomial 4(6) : 3 x2x2 -3 x3x3 -2 x4 +3 x10x10 -3 x11x11 -6 x10x11

Polynomial 5(5) : -2 x5 +3 x10x10 -3 x11x11 +6 x2x3 +6 x10x11

Polynomial 6(5) : 2 x6 +6 x2x10 -6 x2x11 -6 x3x10 -6 x3x11

Polynomial 7(5) : 2 x7 +6 x2x10 +6 x2x11 +6 x3x10 -6 x3x11

Polynomial 8(3) : 2 x8 +6 x0x8 -6 x1x9

Polynomial 9(3) : 2 x9 +6 x0x9 +6 x1x8

Polynomial 10(3) : 2 x10 +6 x4x6 -6 x5x7

Polynomial 11(3) : 2 x11 +6 x4x7 +6 x5x6
'''

frob3 = '''
Polynomial 0(1) : 0x1 x0

Polynomial 1(1) : 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa x1

Polynomial 2(1) : 0x1 x3

Polynomial 3(1) : 0x1 x2

Polynomial 4(1) : 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa x4

Polynomial 5(1) : 0x1 x5

Polynomial 6(2) : 0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x6 +0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x7

Polynomial 7(2) : 0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x6 +0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x7

Polynomial 8(2) : 0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x8 +0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x9

Polynomial 9(2) : 0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x8 +0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x9

Polynomial 10(2) : 0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x10 +0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x11

Polynomial 11(2) : 0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x10 +0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x11
'''

frob2 = '''
Polynomial 0(1) : 0x1 x0

Polynomial 1(1) : 0x1 x1

Polynomial 2(1) : 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe x2

Polynomial 3(1) : 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe x3

Polynomial 4(1) : 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac x4

Polynomial 5(1) : 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac x5

Polynomial 6(1) : 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff x6

Polynomial 7(1) : 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff x7

Polynomial 8(1) : 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa x8

Polynomial 9(1) : 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa x9

Polynomial 10(1) : 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad x10

Polynomial 11(1) : 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad x11
'''

frob1 = '''
Polynomial 0(1) : 0x1 x0

Polynomial 1(1) : 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa x1

Polynomial 2(1) : 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac x3

Polynomial 3(1) : 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac x2

Polynomial 4(1) : 0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad x4

Polynomial 5(1) : 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe x5

Polynomial 6(2) : 0x1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8 x6 +0xfc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3 x7

Polynomial 7(2) : 0xfc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3 x6 +0xfc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3 x7

Polynomial 8(2) : 0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x8 +0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x9

Polynomial 9(2) : 0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 x8 +0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 x9

Polynomial 10(2) : 0x5b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116 x10 +0x144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995 x11

Polynomial 11(2) : 0x144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995 x10 +0x144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995 x11
'''

def to_hex_limbs(x):
    return [hex(x >> (64*i) & (2**64 - 1)) for i in range(6)]

def print_frobenius(name, frobenius_str, frob_num):
    polynomial_idx = 0
    constants = {}
    for line in frobenius_str.split('\n'):
        if 'Polynomial' in line:
            parts = [x.strip() for x in line.split(":")[1].strip().split(" ")]
            assert len(parts) % 2 == 0
            for i in range(0, len(parts), 2):
                coeff = int(parts[i].split("x")[1], base=16)
                term = int(parts[i+1].split("x")[1])
                c = constants.get(coeff, [])
                c.append((polynomial_idx, term))
                constants[coeff] = c
            polynomial_idx += 1
    string = name + " {\n"
    polynomials = [[] for _ in range(polynomial_idx)]
    for i, (coeff, terms) in enumerate(constants.items()):
        if coeff == 1:
            for term in terms:
                polynomials[term[0]].append(f"x[{term[1]}]")
        elif coeff - q == -1:
            for term in terms:
                polynomials[term[0]].append(f"ring.negate(x[{term[1]}])")
        else:
            val = 'ring.wrap(ring.from_hex("' + hex(coeff) + '"))'
            constant_name = f"frob{frob_num}_{i}"
            string += f"\tauto {constant_name} = {val};\n"
            unique_terms = set([t[1] for t in terms])
            for u_term in unique_terms:
                term_name = f"frob{frob_num}_{i}_{u_term}"
                string += f"\tauto {term_name} = x[{u_term}] * {constant_name};\n"
                for term in terms:
                    if term[1] == u_term:
                        polynomials[term[0]].append(term_name)
    requires_reduction = []
    for poly_idx, polynoms in enumerate(polynomials):
        if "frob" in polynoms[0]:
            requires_reduction.append(poly_idx)
        else:
            assert(len(polynoms) == 1)
            continue
        if (len(polynoms) == 1):
            string += f"\tauto poly{poly_idx}_unique = {polynoms[0]};\n"
        else:
           if poly_idx % 2 == 0:
                for i,p in enumerate(polynoms):
                    for j, p2 in enumerate(polynomials[poly_idx + 1]):
                        if p == p2:
                            shared_term = "shared_" + p
                            string += f"\tauto {shared_term} = {p};\n"
                            polynoms[i] = shared_term
                            polynomials[poly_idx+1][j] = shared_term
    for poly_idx, polynoms in enumerate(polynomials):
        if len(polynoms) > 1:
            polynoms.sort(reverse=True)
            assert(polynoms[0].startswith("shared_"))
            string += f"\tauto poly{poly_idx}_unique = " + ' + '.join(polynoms) + ";\n"
    reduced_names = ", ".join([f"reduced_poly{i}" for i in requires_reduction])
    reduction_names = ", ".join([f"poly{i}_unique" for i in requires_reduction])
    string += f"\tauto [{reduced_names}] = ring.template batch_reduce_expand({reduction_names});\n"
    string += f"\treturn std::array<typename Ring::StandardElement, {len(polynomials)}>" + "{" + ", ".join([f"reduced_poly{i}" if i in requires_reduction else polynomials[i][0] for i in range(len(polynomials))]) + "};\n"
    string += "}"
    print(string)
    return string


if __name__ == "__main__":
    print('''#pragma once
#include <array>
#include <cstdint>''')
    print_function("template<class Ring>\ninline auto fp12_flat_mul(const std::array<typename Ring::StandardElement, 12> &x, const std::array<typename Ring::StandardElement, 12> &y, Ring &ring)", fp12)
    print_function("template<class Ring>\ninline auto flat_mulxy0z00(const std::array<typename Ring::StandardElement, 12> &x, const std::array<typename Ring::StandardElement, 6> &y, Ring &ring)", mulx0y0z00)
    print_function("template<class Ring>\ninline auto fp12_flat_sqr(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring)", fp12sqr)
    print_function("template<class Ring>\ninline auto fp12_flat_cyclotomic_sqr(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring)", cyclotomic_sqr_fp12)
 

    print_frobenius("template<class Ring>\ninline auto fp12_flat_frobenius_map1(Ring &ring, const std::array<typename Ring::StandardElement, 12> &x)", frob1, 1)
    print_frobenius("template<class Ring>\ninline auto fp12_flat_frobenius_map2(Ring &ring, const std::array<typename Ring::StandardElement, 12> &x)", frob2, 2)
    print_frobenius("template<class Ring>\ninline auto fp12_flat_frobenius_map3(Ring &ring, const std::array<typename Ring::StandardElement, 12> &x)", frob3, 3)

    print('''template<class Ring>
inline auto fp12_flat_conjugate(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring) {
	return std::array<typename Ring::StandardElement, 12>{x[0], x[1], x[2], x[3], x[4], x[5], ring.negate(x[6]), ring.negate(x[7]), ring.negate(x[8]), ring.negate(x[9]), ring.negate(x[10]), ring.negate(x[11])};
}''')

    print_function("template<class Ring>\ninline auto fp12_flat_inverse_to_fp6(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring)", inv_fp12_to_fp6)
    print_function("template<class Ring>\ninline auto fp6_flat_inverse_to_fp12(const std::array<typename Ring::StandardElement, 12> &x, const std::array<typename Ring::StandardElement, 6> &y, Ring &ring)", inv_fp6_to_fp12)
    print_function("template<class Ring>\ninline auto fp6_flat_inverse_to_fp2_c0_c1_c2(const std::array<typename Ring::StandardElement, 6> &x, Ring &ring)", inv_fp6_to_fp2_c0_c1_c2)
    print_function("template<class Ring>\ninline auto fp6_flat_inverse_to_fp2_ret(const std::array<typename Ring::StandardElement, 6> &x, Ring &ring)", inv_fp6_to_fp2_ret)
    print_function("template<class Ring>\ninline auto fp2_flat_inverse_to_fp6(const std::array<typename Ring::StandardElement, 6> &x, const std::array<typename Ring::StandardElement, 2> &y, Ring &ring)", inv_fp2_to_fp6)
    print_function("template<class Ring>\ninline auto fp2_flat_inverse_to_fp(const std::array<typename Ring::StandardElement, 2> &x, Ring &ring)", inv_fp2_to_fp)
    print_function("template<class Ring>\ninline auto fp_flat_inverse_to_fp2(const std::array<typename Ring::StandardElement, 2> &x, const std::array<typename Ring::StandardElement, 1> &y, Ring &ring)", inv_fp_to_fp2)