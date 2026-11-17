#!/usr/bin/env python3
"""Verify the finite QRAQM exponent using embedded exact parameters.

Only Python's standard library and mpmath are required. No input or output
files are used. These calculations establish arithmetic feasibility under
the algorithmic assumptions stated in the paper; they do not establish a
coherent update implementation or its hard workspace bound.
"""

from fractions import Fraction
import sys

from mpmath import iv, mp

DEPTH = 9
ALPHABET = (-2, -1, 0, 1, 2)
LABELS = 3
PRECISION = 100
TARGET = Fraction(99, 500)
RESERVE = Fraction(1, 2**26)
ROOT_PAIRS = (
    (-2, 2),
    (-1, 1),
    (-1, 2),
    (0, 0),
    (0, 1),
    (1, -1),
    (1, 0),
    (2, -2),
    (2, -1),
)
SUBTREE_PAIRS = (
    (-2, 0),
    (-2, 1),
    (-2, 2),
    (-1, -1),
    (-1, 0),
    (-1, 1),
    (-1, 2),
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 1),
)

# The complete joint law has 185 coordinates: the ordered root, followed
# by eight unordered split laws in each of two equally weighted classes.
ZERO_MASSES = (
    33, 44, 54, 55, 62, 64, 66, 74,
    76, 81, 84, 87, 88, 89, 92, 95,
)
FREE_MASSES = {
    14: "74351584103623/3345121733088535",
    15: "326435434859/2572126947144344",
    16: "7061302622099719/9007199254740992",
    17: "218009962699575/1989419827000804",
    18: "973837751885/9409585759127149",
    19: "5574209779933/8797377194893120",
    25: "33913490737275/1717615190777411",
    26: "290813615392/4806176347421393",
    27: "7750493605115287/9007199254740992",
    28: "599276143839297/9007199254740992",
    29: "71682336273/9922239252011858",
    30: "691053270091/6395589227219087",
    36: "27716150300655/5777563520834824",
    38: "4122652133214493/4503599627370496",
    39: "134371336383299/3105584783790389",
    40: "32834542685/8993296115547776",
    41: "242894100769/7300978517228440",
    47: "12762870959913/9399856781098550",
    49: "4296038487234695/4503599627370496",
    50: "164599802422071/6838852439371928",
    51: "2202665956/8890741517622841",
    52: "19099393322/4073895064672851",
    60: "4394912572384763/4503599627370496",
    61: "97666919881593/7679436652566938",
    63: "2656131203/7585188376345836",
    71: "8895187559894721/9007199254740992",
    72: "60307568070703/9190375785952953",
    73: "3556763/9369366399943945",
    82: "8949722529002611/9007199254740992",
    83: "30921732022837/9191844911552865",
    85: "13225083/5340720181476124",
    91: "161963133764/9285239460467627",
    93: "2244485557837801/2251799813685248",
    94: "13190611124200/7659341718034267",
    96: "540/2178710988894757",
    102: "1139687601818521/9007199254740992",
    103: "115407128185009/6704545726080549",
    104: "3762583978967139/9863248685506357",
    105: "501466446121559/1125899906842624",
    106: "9231448124969/2875657291009983",
    107: "35095797385465/3625038125263774",
    113: "1702535879304672/9530889178959953",
    114: "170490366063078/8772878170250251",
    115: "4221546977109024/9791189183662277",
    116: "2516333532581977/9100710090375066",
    117: "24095755972055/7671615294361569",
    118: "58583415372669/7844992619914190",
    124: "2014126776440233/9007199254740992",
    125: "190441923613269/9480505603481857",
    126: "2290527670677904/6054467375371503",
    127: "1596659152807967/7341898704619154",
    128: "42996816849539/8561135444799341",
    129: "26498965095567/3298174893778235",
    135: "106434112791329/562949953421312",
    136: "138504738960158/8444316251518803",
    137: "3242063469577955/9007199254740992",
    138: "387952348599657/1780748697772061",
    139: "21604600286383/3770393606566225",
    140: "85978657889398/8602719205885591",
    146: "669952358163997/4503599627370496",
    147: "84845569317416/9253273242198809",
    148: "1863800325446746/4542916348481077",
    149: "2091300223239923/9919496486470068",
    150: "20504706926459/4827160390490344",
    151: "39129787646721/4419699167321573",
    157: "572224191991897/5414517760545567",
    158: "34210220275226/8675606061316393",
    159: "2284430969923291/4503599627370496",
    160: "1584594685518181/8392857292535982",
    161: "20978368874456/9789119993022817",
    162: "37116469377453/6944788024591040",
    168: "135715371626485/2147223172012823",
    169: "9506728701557/7803383742365862",
    170: "2830245927826321/4503599627370496",
    171: "345024448219331/2251799813685248",
    172: "5069547209455/6460790402172627",
    173: "17701246418463/7133791043100226",
    175: "76872426426/5174924672762479",
    176: "270098363/9954446624136088",
    177: "6718701101658/6893004903879167",
    178: "820566379722488/7620930756880425",
    179: "76911513830261/3339207832093092",
    180: "61537199689/4252170950793916",
    181: "1702279017264259/2251799813685248",
    182: "1038324299352974/9331793789327569",
    183: "113912721175/6514259466754709",
    184: "6541958049637/6388157105660263",
}

# Each row contains the first two label probabilities as numerators over
# 2**60. The third is the exact complement. Only depths 2 through 5 need
# tables: all subsequent filters use the constant label zero.
FILTER_DENOMINATOR = 2**60
FILTER_NUMERATORS = (
    (
        (  # Depth 2.
            (1152521763494107600, 0),
            (1152278915004699936, 418301859374826),
            (1151495948966236903, 207719537),
            (1151639969565205274, 1281535041641702),
            (1567648312673821, 86367093298),
            (2788424367644832, 1150132989885100467),
            (1152878254996345166, 11539498862629),
            (1152916018880799897, 6639512119),
            (1152919194550999153, 1839835188067),
            (1152895968019840631, 10074310),
            (1152919194816764817, 2309503776393),
            (1151980508327357127, 98116334),
            (1152624314662437496, 297178839563450),
            (1152919315594071301, 420825046419),
            (1152913689818737452, 33369876),
            (1152919427103395578, 1548088046645),
            (1152912936597160218, 710293),
            (1152621842572334067, 299662034512909),
            (1152882041355872434, 704408901405),
        ),
        (  # Depth 3.
            (1142381958991315238, 2722549606023),
            (1056248508541634047, 96637276330354432),
            (654399725906264592, 105047561828),
            (277831340369167622, 875090051130056238),
            (7979328214573927, 0),
            (13706858932444276, 1139214362300910330),
            (1152450137186392723, 307195195589748),
            (1152513647222912474, 1109670281562),
            (1150030429615485971, 2886320005003150),
            (1082518425443699181, 38515873844),
            (1058873161798558084, 94048342808288892),
            (251036358570498421, 76891181378),
            (503291788349278278, 649629629088197384),
            (1152895595334168978, 14206233879846),
            (1150747419178267569, 0),
            (1152420788776214046, 497636260437122),
            (1046239729072314170, 0),
            (1139535458072574982, 13386046534271994),
            (1152394058522718259, 18853193160962),
        ),
        (  # Depth 4.
            (377023165648350613, 356643723080355412),
            (332020895117249883, 430532564995432912),
            (281609631068699527, 167220936858979),
            (302745854800750161, 848275550597101487),
            (1152921176871786135, 109534746729),
            (1152921283018245641, 132503131639),
            (1102445658006947034, 49187922768835110),
            (1143370724833444842, 84390759627930),
            (1137135051783365088, 15738181233153568),
            (967439014328716141, 5061954963471),
            (1038513128895464789, 114395542155129003),
            (253862323464228598, 692696570850862),
            (233175002557065165, 919618036505324083),
            (1151361085200879449, 598355768512167),
            (1129678873801038550, 12531542319916),
            (1147134292273985992, 5628794926518584),
            (489729963983520357, 214290615281775732),
            (553972887264443666, 395571263229977055),
            (1136049573802764342, 730998626369345),
        ),
        (  # Depth 5.
            (1149998552913979460, 2067617302307388),
            (1150110938165193053, 1295457627080099),
            (1150021166008875925, 1107507306877803),
            (1150114675103884013, 1307366184499219),
            (1150833966921601291, 1000359082324980),
            (1150664129115786590, 891197153781410),
            (1150279394326608285, 1244178874926436),
            (1149345029046419635, 2087981863458636),
            (1148434198201349391, 2872954545493746),
            (1150762712554870299, 1256135631102565),
            (1150642535603995977, 1298669458820151),
            (1150545259222439569, 1198180765572207),
            (1150990386417590427, 990840097783653),
            (1150636241528645962, 1157416712844854),
            (1150028429776604577, 1480037944755679),
            (1150964677291687723, 941328718554965),
            (1150681608670818329, 1174308670127591),
            (1149318569991340929, 1757960268806782),
            (1150194132864054734, 1402717223550386),
        ),
    ),
    (
        (  # Depth 2.
            (288931608209568384, 0),
            (310020799880770182, 837285466208504080),
            (44982783260972013, 89463732953),
            (62432919787927860, 1090488584818919116),
            (1949255690261867, 61903767691),
            (7818445204564484, 1145102735968856096),
            (1115739495025014786, 33279137260467392),
            (771565617995894108, 1513178184792886),
            (931898791809524192, 220056256082536375),
            (235998882530139431, 332549484732),
            (452878934431134939, 700042459841517511),
            (11803504835261509, 188898166539),
            (82546351008785524, 1070375016537764575),
            (1043814077465250592, 19608260950096157),
            (687393322771310921, 575743689435782),
            (1014913782157507061, 124775774247276335),
            (64652002818902159, 84888804999),
            (567466185017790592, 585455319589056384),
            (1040195484841243746, 5705917548176209),
        ),
        (  # Depth 3.
            (730108200210392080, 159441954078183),
            (392638819905060083, 760192849814233192),
            (168434129419093454, 117354289690),
            (91450255207879724, 1061471100678589646),
            (17716140335859714, 0),
            (19587669420471108, 1133333666179240444),
            (1110520473326942930, 35224341124260563),
            (1015132579941025231, 2531935173216071),
            (955877771075084893, 195791326308532671),
            (491859479814360428, 263068955482),
            (515999861257372067, 636921643349474909),
            (73279429236745919, 164170159563),
            (161180175618398926, 991741180850346897),
            (1106308381003646049, 17893614540513987),
            (913700109330457319, 0),
            (1045996411492900099, 101953426958436205),
            (297984384850579968, 0),
            (775990840757942827, 376930663848904149),
            (1100095390610657542, 88746116013216),
        ),
        (  # Depth 4.
            (350918135602153260, 273459348721260977),
            (212463369094726464, 557612578482325056),
            (195717931448519204, 200152681627081),
            (121630058347395172, 1027071040752441501),
            (10657692384228699, 47813880384800674),
            (7273149399424030, 1093733866591056995),
            (1136394350552775603, 15005589510709069),
            (1118834199281343518, 319142883951581),
            (1000570975665404288, 152187875674303104),
            (638599683120104029, 0),
            (633035995682403584, 519871709129782528),
            (147759366343155539, 868005867902127),
            (414575242992666025, 738051446699632983),
            (1128347511993028522, 16367037684826201),
            (1010626674924012825, 103985095283681),
            (1061528890401181494, 90775675542731082),
            (577448803018957160, 84563100512252908),
            (984515545909317182, 98141076572289734),
            (1132090173012219044, 927455010647121),
        ),
        (  # Depth 5.
            (88814469206136483, 591504493776185744),
            (13849359729254830, 652611626655104468),
            (64928031168612514, 317865683888926025),
            (19157273949384773, 832059741249029291),
            (20584897295634352, 460055412741799239),
            (3323116302524347, 631817198600876931),
            (181982456409305607, 635823618597530408),
            (869474221087367662, 265667323446076179),
            (492893465700489691, 656642076483143461),
            (717052641696440197, 65582034154549496),
            (292856970612004096, 857814509250392832),
            (12941182151119945, 346468774259095441),
            (19639627618289340, 921020479338216761),
            (741944553988688512, 409880224852115456),
            (852618319652409344, 180938058284937856),
            (597389575320845312, 554827418094508032),
            (33054032939148835, 518589627358454512),
            (43347278248105841, 751755454697115146),
            (548329651870416530, 327365179854655845),
        ),
    ),
)


def require(condition, message):
    """Make verification failures explicit, including under Python -O."""
    if not condition:
        raise ValueError(message)


def marginal(law, positions):
    result = {}
    for entry, probability in law.items():
        key = tuple(entry[position] for position in positions)
        result[key] = result.get(key, Fraction(0)) + probability
    return result


def equal_laws(first, second):
    return all(first.get(key, 0) == second.get(key, 0)
               for key in first.keys() | second.keys())


def solve_linear_system(matrix, right):
    """Solve a square nonsingular system by exact rational elimination."""
    size = len(right)
    require(len(matrix) == size and all(len(row) == size for row in matrix),
            "Invalid linear-system dimensions")
    rows = [list(row) + [value] for row, value in zip(matrix, right)]
    for column in range(size):
        pivot = next((r for r in range(column, size) if rows[r][column]), None)
        require(pivot is not None, "The specified free masses are not unique")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        divisor = rows[column][column]
        rows[column] = [entry / divisor for entry in rows[column]]
        for r in range(column + 1, size):
            multiplier = rows[r][column]
            if multiplier:
                rows[r] = [a - multiplier * b
                           for a, b in zip(rows[r], rows[column])]
    solution = [Fraction(0)] * size
    for r in range(size - 1, -1, -1):
        solution[r] = rows[r][-1] - sum(
            rows[r][j] * solution[j] for j in range(r + 1, size))
    return solution


def mass_index(side, level, pair_index):
    """Index an unordered split mass; levels run from 2 through DEPTH."""
    return len(ROOT_PAIRS) + ((side * (DEPTH - 1) + level - 2)
                             * len(SUBTREE_PAIRS) + pair_index)


def reconstruct_masses():
    """Complete the 87 free masses using all exact parent/child equations."""
    count = len(ROOT_PAIRS) + 2 * (DEPTH - 1) * len(SUBTREE_PAIRS)
    require(count == 185 and len(ZERO_MASSES) == 16 and len(FREE_MASSES) == 87,
            "Incorrect probability architecture")
    equations, right = [], []
    for total in (0, 1):
        row = [Fraction(0)] * count
        for index, pair in enumerate(ROOT_PAIRS):
            row[index] = Fraction(sum(pair) == total)
        equations.append(row)
        right.append(Fraction(1, 2))
    for side in range(2):
        for level in range(2, DEPTH + 1):
            for digit in ALPHABET:
                row = [Fraction(0)] * count
                for index, pair in enumerate(SUBTREE_PAIRS):
                    row[mass_index(side, level, index)] += (sum(pair) == digit)
                if level == 2:
                    for index, pair in enumerate(ROOT_PAIRS):
                        row[index] -= (pair[side] == digit)
                else:
                    for index, (a, b) in enumerate(SUBTREE_PAIRS):
                        row[mass_index(side, level - 1, index)] -= Fraction(
                            int(a == digit) + int(b == digit), 2)
                equations.append(row)
                right.append(Fraction(0))
    fixed = {index: Fraction(0) for index in ZERO_MASSES}
    require(not (fixed.keys() & FREE_MASSES.keys()), "Conflicting mass entries")
    fixed.update({index: Fraction(value) for index, value in FREE_MASSES.items()})
    unknown = [index for index in range(count) if index not in fixed]
    require(len(equations) == len(unknown) == 82, "Invalid type-flow rank")
    matrix = [[row[index] for index in unknown] for row in equations]
    residual = [value - sum(row[index] * mass for index, mass in fixed.items())
                for row, value in zip(equations, right)]
    fixed.update(zip(unknown, solve_linear_system(matrix, residual)))
    values = [fixed[index] for index in range(count)]
    require(min(values) >= 0, "A reconstructed mass is negative")
    for row, value in zip(equations, right):
        require(sum(a * b for a, b in zip(row, values)) == value,
                "An exact type-flow equation fails")
    return values


def load_laws():
    """Build and check every joint type and conditional filter channel."""
    masses = reconstruct_masses()
    root = dict(zip(ROOT_PAIRS, masses[:len(ROOT_PAIRS)]))
    require(sum(root.values()) == 1, "The root law is not normalized")
    require(all(a in ALPHABET and b in ALPHABET and a + b in (0, 1)
                for a, b in root), "Invalid root support")
    root_types = [marginal(root, (side,)) for side in range(2)]
    types, joints, extensions = [], [], []
    for side in range(2):
        parent = root_types[side]
        side_types, side_joints, side_extensions = [parent], [], []
        for level in range(2, DEPTH + 1):
            joint, extension = {}, {}
            ordered_index = 0
            for index, (a, b) in enumerate(SUBTREE_PAIRS):
                require(a <= b and a in ALPHABET and b in ALPHABET
                        and a + b in ALPHABET, "Invalid split support")
                mass = masses[mass_index(side, level, index)]
                ordered = ((a, b, mass),) if a == b else (
                    (a, b, mass / 2), (b, a, mass / 2))
                for x, y, probability in ordered:
                    joint[x, y] = probability
                    if level <= 5:
                        first, second = FILTER_NUMERATORS[side][level - 2][ordered_index]
                        channel = [Fraction(first, FILTER_DENOMINATOR),
                                   Fraction(second, FILTER_DENOMINATOR)]
                        channel.append(1 - sum(channel))
                    else:
                        channel = [Fraction(1), Fraction(0), Fraction(0)]
                    require(min(channel) >= 0 and sum(channel) == 1,
                            "An invalid filter channel was specified")
                    for label, value in enumerate(channel):
                        extension[x, y, label] = probability * value
                    ordered_index += 1
            require(ordered_index == 19 and sum(joint.values()) == 1,
                    "An invalid split law was specified")
            sums = {}
            for (a, b), probability in joint.items():
                sums[a + b,] = sums.get((a + b,), Fraction(0)) + probability
            require(equal_laws(sums, parent), "A parent sum law does not match")
            parent = marginal(joint, (0,))
            require(equal_laws(parent, marginal(joint, (1,))),
                    "The two child marginal laws differ")
            require(equal_laws(marginal(extension, (0, 1)), joint),
                    "A filter extension does not preserve its useful law")
            side_types.append(parent)
            side_joints.append(joint)
            side_extensions.append(extension)
        types.append(side_types)
        joints.append(side_joints)
        extensions.append(side_extensions)
    return root, root_types, types, joints, extensions


def interval(value):
    value = Fraction(value)
    return iv.mpf(value.numerator) / value.denominator


def entropy(law):
    return -sum((interval(p) * iv.log(interval(p), 2)
                 for p in law.values() if p), interval(0))


def conditional_information(law, first, second, given=()):
    """Evaluate conditional mutual information directly from exact ratios."""
    joint = marginal(law, tuple(first) + tuple(second) + tuple(given))
    left = marginal(law, tuple(first) + tuple(given))
    right = marginal(law, tuple(second) + tuple(given))
    conditioning = marginal(law, tuple(given))
    first_size, second_size = len(first), len(second)
    result = interval(0)
    for key, probability in joint.items():
        if not probability:
            continue
        x = key[:first_size]
        y = key[first_size:first_size + second_size]
        z = key[first_size + second_size:]
        ratio = probability * conditioning[z] / (left[x + z] * right[y + z])
        if ratio != 1:
            result += interval(probability) * iv.log(interval(ratio), 2)
    return result


def feature_values(laws):
    """Compute entropy features and independently cross-check information laws."""
    root, root_types, types, joints, extensions = laws
    result = {"one": interval(1)}
    for level in range(1, DEPTH + 1):
        result[f"h{level}"] = sum(entropy(types[side][level - 1])
                                  for side in range(2)) / 2
    result["f1"] = entropy(root) - sum(entropy(p) for p in root_types)
    root_information = conditional_information(root, (0,), (1,))
    require((result["f1"] + root_information).a <= 0 <=
            (result["f1"] + root_information).b, "Root information identity failed")
    for level in range(2, DEPTH + 1):
        values = {name: [] for name in ("f", "d", "q", "r")}
        for side in range(2):
            joint = joints[side][level - 2]
            extension = extensions[side][level - 2]
            child = types[side][level - 1]
            f = entropy(joint) - 2 * entropy(child)
            alpha = conditional_information(extension, (0,), (2,))
            beta = conditional_information(extension, (1,), (2,))
            kappa = conditional_information(extension, (0, 1), (2,))
            d = conditional_information(extension, (0,), (2,), (1,))
            q = conditional_information(extension, (1,), (2,), (0,))
            r = conditional_information(extension, (0,), (1,), (2,))
            for difference in (d - (kappa - beta), q - (kappa - alpha),
                               r + f - (kappa - alpha - beta)):
                require(difference.a <= 0 <= difference.b,
                        "A mutual-information chain-rule check failed")
            for name, value in zip(("f", "d", "q", "r"), (f, d, q, r)):
                values[name].append(value)
        for name, entries in values.items():
            result[f"{name}{level}"] = sum(entries) / 2
    return result


FEATURE_NAMES = (["one"] + [f"h{j}" for j in range(1, DEPTH + 1)]
                 + [f"f{j}" for j in range(1, DEPTH + 1)]
                 + [f"{name}{j}" for name in ("d", "q", "r")
                    for j in range(2, DEPTH + 1)])


def basis(name):
    return tuple(Fraction(entry == name) for entry in FEATURE_NAMES)


def add(*vectors):
    return tuple(sum(entries) for entries in zip(*vectors))


def scale(value, vector):
    return tuple(Fraction(value) * entry for entry in vector)


def subtract(first, second):
    return add(first, scale(-1, second))


def resource_prescription():
    """Construct the common-height rates as exact symbolic expressions."""
    one = basis("one")
    l0, u0 = scale(-1, basis("f2")), basis("d2")
    t0 = add(scale(-1, basis("f3")), basis("d3"), basis("q3"),
             scale(-1, basis("r3")))
    height = add(l0, scale(RESERVE, one))
    marked = add(scale(-1, l0), scale(-2 * RESERVE, one))
    update = add(u0, scale(RESERVE, one))
    total = add(t0, scale(3 * RESERVE, one))
    moduli = [one]
    for level in range(1, DEPTH + 1):
        moduli.append(add(one, scale(-(level + 2), l0),
                          *(scale(-1, basis(f"f{j}"))
                            for j in range(1, level + 1)),
                          scale(-(level + 3) * RESERVE, one)))
    lists = [None] * (DEPTH + 1)
    lists[DEPTH] = height
    for level in range(DEPTH - 1, -1, -1):
        lists[level] = add(scale(2, lists[level + 1]),
                           scale(-1, moduli[level]), moduli[level + 1],
                           basis(f"f{level + 1}"))
    require(lists[0] == marked and all(value == height for value in lists[1:]),
            "The exact common-height recurrence failed")
    require(subtract(moduli[1], moduli[2]) == scale(RESERVE, one),
            "The first modulus gap is not exactly the prescribed reserve")
    return moduli, height, marked, update, total


def interval_maximum(values):
    return iv.mpf([max(value.a for value in values),
                   max(value.b for value in values)])


def verify():
    """Verify every numerical source, setup, update, capacity and budget term."""
    iv.dps = PRECISION
    laws = load_laws()
    features = feature_values(laws)
    moduli, height, marked, update, total = resource_prescription()
    one = basis("one")
    rows = {}

    def at_most(name, budget, cost):
        rows[name] = subtract(budget, cost)

    def positive(name, expression):
        rows[name] = expression

    positive("height", height)
    positive("negative marked exponent", scale(-1, marked))
    positive("update budget", update)
    positive("total budget", total)
    at_most("root modulus", one, moduli[1])
    positive("bottom modulus", moduli[DEPTH])
    setup_costs = [moduli[DEPTH], height, scale(Fraction(1, 2), height),
                   scale(Fraction(1, 2), add(scale(2, height),
                                            scale(-1, one), moduli[1]))]
    storage_costs = [moduli[DEPTH], height]
    update_costs = [scale(0, one), scale(Fraction(1, 2),
                                       add(height, scale(-1, one), moduli[1]))]
    for level in range(2, DEPTH + 1):
        d, q, r = (basis(f"{name}{level}") for name in ("d", "q", "r"))
        family = add(scale(-1, basis(f"f{level}")), d, q, scale(-1, r))
        setup_costs.extend((family, add(height, d), height,
                            add(height, scale(Fraction(1, 2), q)),
                            add(height, scale(Fraction(1, 2), r))))
        storage_costs.extend((family, add(height, d), height))
        update_costs.extend((d, scale(Fraction(1, 2), q),
                             scale(Fraction(1, 2), r)))
        positive(f"capacity {level}",
                 add(basis(f"h{level}"), scale(-1, moduli[level]),
                     scale(Fraction(1, 2**level), marked), scale(-1, height)))
    positive("capacity 1", add(basis("h1"), scale(-1, moduli[1]),
                               scale(Fraction(1, 2), marked), scale(-1, height)))
    for level in range(1, DEPTH):
        positive(f"modulus gap {level}", subtract(moduli[level], moduli[level + 1]))

    # The complete left-source carrier is charged explicitly, even though
    # it is cheaper than the fixed right source for these parameters.
    carrier = subtract(basis(f"h{DEPTH}"), moduli[DEPTH])
    setup_costs.append(carrier)
    storage_costs.append(carrier)
    positive("left source contains sampled subset", subtract(carrier, height))
    positive("right source dominates carrier", subtract(moduli[DEPTH], carrier))
    positive("carrier below total budget", subtract(total, carrier))
    for index, cost in enumerate(setup_costs):
        at_most(f"setup charge {index}", total, cost)
    for index, cost in enumerate(storage_costs):
        at_most(f"storage charge {index}", total, cost)
    for index, cost in enumerate(update_costs):
        at_most(f"update charge {index}", update, cost)
    walk = add(scale(Fraction(1, 2), subtract(height, marked)), update)
    at_most("declared walk charge", total, walk)

    def evaluate(expression):
        return sum((interval(coefficient) * features[name]
                    for name, coefficient in zip(FEATURE_NAMES, expression)
                    if coefficient), interval(0))

    slacks = {}
    for name, expression in rows.items():
        value = evaluate(expression)
        require(value.a > 0, f"Nonpositive resource margin: {name}")
        slacks[name] = value
    declared = evaluate(total)
    require(declared.b < interval(TARGET).a,
            "The declared total exponent does not improve the target")
    sigma = interval_maximum([evaluate(cost) for cost in setup_costs])
    storage = interval_maximum([evaluate(cost) for cost in storage_costs])
    source_fraction = evaluate(moduli[DEPTH]) / features[f"h{DEPTH}"]
    require(source_fraction.a > 0 and source_fraction.b < 1,
            "The source-support proportion lies outside (0,1)")
    return {"declared": declared, "setup": sigma, "storage": storage,
            "walk": evaluate(walk), "update": evaluate(update),
            "carrier": evaluate(carrier), "target_margin": interval(TARGET) - declared,
            "minimum_slack": min((value.a for value in slacks.values())),
            "checked_resource_inequalities": len(slacks)}


def main():
    result = verify()
    print("PASS: exact joint types, filter channels and common-height identities")
    print(f"PASS: {result['checked_resource_inequalities']} strict resource inequalities "
          f"at {PRECISION}-digit interval precision")
    for label, name in (("Declared time and conditional total-space budget", "declared"),
                        ("Setup exponent", "setup"),
                        ("Persistent storage exponent", "storage"),
                        ("Declared walk exponent", "walk"),
                        ("Strict target margin", "target_margin"),
                        ("Minimum strict resource margin", "minimum_slack")):
        with mp.workdps(35):
            approximate = mp.nstr(mp.mpf(result[name].mid), 25)
        print(f"{label} (approx.): {approximate}")
    print("Conditional exponent bound: 0.198")
    print("The list, coherent-access/update and hard-workspace assumptions remain separate.")


if __name__ == "__main__":
    try:
        main()
    except (ArithmeticError, ValueError) as error:
        sys.exit(f"FAIL: {error}")
