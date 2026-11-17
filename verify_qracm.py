#!/usr/bin/env python3
"""Verify the finite QRACM subset-sum exponent using exact probabilities.

Run with Python 3 and mpmath 1.3.0. All parameters are embedded below.
Only arithmetic feasibility is checked. The list-distribution, survival,
occurrence and coherent-access assumptions are those stated in the paper.
No optimizer, external data, network access or output files are used.
"""

from fractions import Fraction
from mpmath import iv


ALPHABET = (-1, 0, 1, 2)
DEPTH = 4
SCALE = 10**14
TARGET = Fraction(2099947, 10000000)
BOTTOM_HEIGHT = Fraction(1123259302947, 5349021095909)
SOURCE_FRACTION = Fraction(1, 3)
PRECISION = 100

# For each coordinate class, the root row contains five nonfinal split
# numerators. Subsequent rows, at nodes 1 through 14, contain eight each.
# Split groups are ordered by their sum, then by the left digit. The final
# probability in each group is one minus the sum of the listed values.
ROOT_SPLITS = [[2088997030214, 94752084228903, 53595732734, 26214026038727, 73191864108997],
 [8397162296039, 88403226034021, 653898122655, 60222999079699, 38755663385201]]

INTERNAL_SPLITS = [[[36656712640765, 1319022067060, 95869717173907, 37670565982, 29285950304091, 68479764589823,
   245360854299, 94138602581846],
  [48876284005284, 541199187594, 98937200925245, 2870601076, 52153306997296, 47831823871489,
   1075205778226, 97994154033954],
  [52510856030195, 481765703089, 98786468700965, 5438882533, 33072750938878, 66895650279652, 139924195800,
   90200917527859],
  [54550140947734, 2011995317063, 96652701415188, 117485414695, 52456315650390, 47376632663419,
   9141488762612, 85866905006901],
  [49095077370540, 386758288624, 99220426351414, 2025740681, 51318764458062, 48669926392334,
   1024456938341, 96375096800367],
  [45915153286503, 396642692126, 99023631275808, 1444927292, 51090476490593, 48866798612897, 698340004933,
   95345570930165],
  [12235955118313, 1512523170, 99676386067339, 0, 3258097654775, 96741902345225, 0, 100000000000000],
  [70923138470502, 119148857716, 99860904910603, 0, 71660321495769, 28339678504231, 0, 100000000000000],
  [49999989027034, 181058253673, 99637893352315, 0, 49999998846332, 50000001153668, 1970833,
   99999998029167],
  [49999971546351, 164275323622, 99671439499887, 610, 50000022757959, 49999977241431, 0, 99999998732612],
  [50001460555776, 80720152922, 99838565727861, 806437, 50000123014567, 49999876178996, 0,
   99999557088546],
  [49998982361848, 79246810384, 99841494000852, 0, 49999894098881, 50000105755741, 36364990,
   99995211069510],
  [50001199729300, 90640470167, 99818723178603, 8338553, 50000082264273, 49999906126929, 17607019518,
   99977213164875],
  [50000694025308, 107942530589, 99784086538752, 5816942, 50000133516631, 49999860666427, 14024784251,
   99984013572452]],
 [[53213618565699, 1245922405989, 96913794007014, 16801251696, 37120607831320, 61545447804640,
   1910889861518, 89960485588905],
  [41545412029281, 3263081496861, 88710956776299, 345976482315, 42217559290332, 56204541393108,
   4138530745687, 88577241628477],
  [60639897769686, 504423527409, 98973193714720, 606960378, 38498550185679, 61424436662711, 142674047033,
   72569582249645],
  [53454156204330, 1554768912069, 97336877847045, 88716571748, 51987707764227, 47879152113718,
   9135641183370, 85459693630871],
  [44810650704569, 2622407342508, 90748888866945, 95586722904, 44297983932102, 54613016747676,
   5063194795077, 82973841999395],
  [40444044701293, 2385251146766, 90284989091066, 67342095167, 40532727568499, 58376991668335,
   3563845415855, 81270221518496],
  [62264821020653, 166444529333, 99801560700621, 0, 75919813975578, 24080186024422, 0, 100000000000000],
  [55156908889518, 29299356807, 99905783198604, 0, 27025632097590, 72974367902410, 0, 100000000000000],
  [50000004091812, 132915550146, 99734162861462, 0, 50000005663217, 49999994335365, 0, 100000000000000],
  [50000154033121, 125391624711, 99749215675371, 0, 50000017663939, 49999982336061, 0, 99999960798301],
  [49999553762921, 1076733871758, 97846515776929, 6721044477, 49993311553182, 49993247183148,
   1191742667177, 97615633476494],
  [50000167284952, 1227048488586, 97545915018800, 898243009, 49999084936093, 49999114654753, 142961571436,
   99713579037389],
  [49999770567089, 1158947658253, 97682096437450, 7459033502, 49992567424586, 49992510673347,
   1226480351173, 97546859631355],
  [50000022679492, 1230285978991, 97539435651513, 1849823619, 49998208731760, 49998088694303,
   288346539085, 99419403751882]]]

# Filter rows use lexicographic ordered pairs with sum in ALPHABET.
# At nodes 1 through 6, each row gives the first two label numerators.
# Label 2 has the complementary probability. Other filters are constant.
FILTER_NUMERATORS = [[[[2976984277474, 2066997], [49932639067787, 13983], [99204838596881, 161990], [915091621, 10777520275],
   [28533089512, 37050746], [919371953374, 101995], [53793659058549, 29562191], [448745, 7967866201296],
   [30327486, 29762751819], [986188953, 11666079], [19370, 99416942940575], [145506666, 37137252238450]],
  [[1758164, 48354982511751], [1734275, 48350896254210], [0, 50712087708874],
   [9949274230, 48312232977364], [28190894, 48373120829153], [905386, 48371039938329],
   [50507145, 50793146555760], [533009377194, 48089784078013], [1520685376, 48375913487239],
   [41351831, 48466971580253], [78548267676589, 8743184570834], [16930217586447, 27091790153504]],
  [[53175377861364, 178576901569], [43356736724976, 18588748737121], [4725895037665, 91585161972272],
   [53344321257389, 25617846], [53297126909439, 2298861963], [53132619185743, 292487716991],
   [49428356216839, 10011943505847], [53352636620451, 321664], [53305534256697, 18738440],
   [53295453939021, 2391417632], [53643447251466, 1542619], [53417619953508, 2505371]],
  [[367999698504, 52167340404045], [15356369160671, 50095086804466], [22979262127313, 35531834308388],
   [2485083818447, 47221741572769], [39622517784, 51245748723891], [5023773813632, 51633074345036],
   [22837677969251, 35529033262532], [4359327210384, 46824771519622], [330070627828, 50989142326233],
   [1797086544070, 52724161058035], [3241250696251, 47145857805061], [72401411136, 4129745669733]],
  [[17951861, 44323238458011], [8632238, 44338550942534], [1509180288, 84178077516735],
   [95074792033, 44269339040645], [64397444, 44308753468081], [3756082, 44321044861218],
   [0, 86978434850037], [3750002170399, 42650439233069], [2632584938, 44303727896507],
   [1155744882, 44236773199715], [53591281771937, 20248946345306], [84260664120, 43324681531118]],
  [[43192778821877, 56807205554758], [43200423178646, 56799571537868], [51349379866020, 48648729341109],
   [42575549115663, 56045458308745], [43183005592482, 56816169434322], [43206960879903, 56792958484053],
   [58283623679209, 41633062666656], [32495101730661, 42837209669515], [43152543491509, 56828169713495],
   [43190313595540, 56806883346481], [5297424697932, 7174522622372], [43315199772728, 56309117627664]]],
 [[[48519216454117, 52196], [48525167533772, 17372], [48524259870217, 0], [48499356970517, 31054207405],
   [48516725255547, 140537748], [48519373216300, 329128], [48527208743098, 5828258],
   [46362660384999, 4427622443748], [48503958003668, 20943228946], [48508860971394, 15005610],
   [1664668251027, 96556236703900], [42862657733953, 11212146513957]],
  [[17876336403, 96097817410872], [49743918, 52263280473959], [14014814, 1071678480042],
   [10826551553478, 89165737583864], [336621312448, 99407376122569], [1640116100, 94500620198658],
   [1459829972, 14589265083911], [67113310736816, 32886615371607], [5384519839162, 94609127897763],
   [29103331116, 99819478613096], [98585425872113, 1414573992491], [66100946551190, 33898912697486]],
  [[51073262209630, 20753318017109], [51078876720845, 20752930813310], [33341770788578, 36661996395756],
   [51072062482880, 20750451952578], [51073736401297, 20753047988312], [51075132302737, 20754129172632],
   [39420888473414, 30684864893188], [51070896568407, 20750255250785], [51074815424040, 20751414503313],
   [51076735890840, 20750522134814], [51088239972553, 20739003198769], [51082863035480, 20744173476133]],
  [[122887885122, 51868784668371], [14795410443043, 49781835880147], [21075941882277, 31709067504921],
   [3334910433602, 47076729050547], [4741346166, 51110067633285], [6060512486784, 51039171576356],
   [23494354451055, 30683849219767], [5367253896858, 46620248542324], [21503216343, 50881199759416],
   [963544530504, 53055228185437], [4189145007656, 46771150617267], [82472974551, 3993616290122]],
  [[10314078673, 53566876651116], [286426822, 53553760960753], [207270223, 53803028500573],
   [8685677522889, 48950215360038], [130227698248, 53518420733450], [3627461024, 53571764824088],
   [211149169, 53927650978498], [52695041851808, 25360687753317], [1504025427312, 52787253800796],
   [42461570389, 53554434370627], [84660445785571, 8246527391937], [7040498688672, 50095233165892]],
  [[44999684670434, 54971450004605], [45031628048354, 54966884962081], [47597626058995, 52402048314701],
   [39796211887703, 48673379124175], [44886053989240, 54865663265171], [45013737604637, 54973440387434],
   [48324711493724, 51674257577399], [24373087568190, 29814823557042], [44274794159002, 54136172495658],
   [44964606356522, 54952465623126], [10828995379487, 13122704685732],
   [42298285679133, 51894220542132]]]]

# K[v] is the cumulative modulus shared by the two children of v.
CHILD_MODULI = {
    1: Fraction(4438413308952, 8663522859127),
    2: Fraction(805845594465, 1431843611537),
    3: Fraction(2474210743293, 8111807237746),
    4: Fraction(3126790979971, 9359904369534),
    5: Fraction(3638186018853, 8788687160918),
    6: Fraction(4105473057359, 9921465150964),
    7: Fraction(246098898929, 2589967553085),
    8: Fraction(491812915853, 5115678264976),
    9: Fraction(508431299265, 4050813590344),
    10: Fraction(603544620740, 4808780232153),
    11: Fraction(1083233820766, 5159855929833),
    12: Fraction(1462667204314, 6965277142263),
    13: Fraction(1931067219149, 9199765101091),
    14: Fraction(1095319666365, 5215954155814),
}


def require(condition, message):
    """Raise an explicit error even if Python assertions are disabled."""
    if not condition:
        raise ArithmeticError(message)


def normalized(values, message):
    require(all(p >= 0 for p in values), message + ": negative probability")
    require(sum(values) == 1, message + ": probabilities do not sum to one")


def marginal(law, positions):
    result = {}
    for symbols, probability in law.items():
        key = tuple(symbols[position] for position in positions)
        result[key] = result.get(key, Fraction(0)) + probability
    return result


def reconstruct_probabilities():
    """Reconstruct every exact joint law, marginal and filter extension.

    The two physical coordinate classes have weight one half. Root sums
    0 and 1 each have probability one half. All other parent types follow
    from the displayed conditional splits and their exact child marginals.
    """
    require(len(ROOT_SPLITS) == len(INTERNAL_SPLITS) == len(FILTER_NUMERATORS) == 2,
            "Expected two coordinate classes")
    internal_count = 2**DEPTH - 1
    node_count = 2**(DEPTH + 1) - 1
    pair_order = tuple((x, y) for x in ALPHABET for y in ALPHABET
                       if x + y in ALPHABET)
    all_joints, all_types, all_filters = [], [], []
    for coordinate_class in range(2):
        require(len(ROOT_SPLITS[coordinate_class]) == 5, "Root split dimension")
        require(len(INTERNAL_SPLITS[coordinate_class]) == 14,
                "Internal split dimension")
        require(len(FILTER_NUMERATORS[coordinate_class]) == 6,
                "Filter dimension")
        types = [{0: Fraction(1, 2), 1: Fraction(1, 2)}] + [None] * (node_count - 1)
        joints, filters = [], []
        for node in range(internal_count):
            row = (ROOT_SPLITS[coordinate_class] if node == 0
                   else INTERNAL_SPLITS[coordinate_class][node - 1])
            require(len(row) == (5 if node == 0 else 8), "Split row dimension")
            require(all(isinstance(x, int) and 0 <= x <= SCALE for x in row),
                    "Invalid split numerator")
            cursor, joint = 0, {}
            parent_sums = (0, 1) if node == 0 else ALPHABET
            for parent_sum in parent_sums:
                pairs = [(x, parent_sum - x) for x in ALPHABET
                         if parent_sum - x in ALPHABET]
                count = len(pairs) - 1
                probabilities = [Fraction(x, SCALE) for x in row[cursor:cursor + count]]
                cursor += count
                probabilities.append(1 - sum(probabilities))
                normalized(probabilities, "Conditional split")
                for pair, probability in zip(pairs, probabilities):
                    joint[pair] = types[node].get(parent_sum, Fraction(0)) * probability
            require(cursor == len(row), "Unused split parameters")
            normalized(list(joint.values()), "Joint type")
            sum_type = {}
            for (x, y), probability in joint.items():
                sum_type[x + y] = sum_type.get(x + y, Fraction(0)) + probability
            require(sum_type == types[node], "Parent sum type mismatch")
            left_type = {x: sum(p for (a, b), p in joint.items() if a == x)
                         for x in ALPHABET}
            right_type = {y: sum(p for (a, b), p in joint.items() if b == y)
                          for y in ALPHABET}
            types[2 * node + 1], types[2 * node + 2] = left_type, right_type
            normalized(list(left_type.values()), "Left marginal")
            normalized(list(right_type.values()), "Right marginal")
            channel = {}
            if 1 <= node <= 6:
                filter_rows = FILTER_NUMERATORS[coordinate_class][node - 1]
                require(len(filter_rows) == len(pair_order), "Filter row count")
                require(set(joint) == set(pair_order), "Filter support mismatch")
                for pair, numerators in zip(pair_order, filter_rows):
                    require(len(numerators) == 2 and all(isinstance(x, int) for x in numerators),
                            "Filter row dimension")
                    first = tuple(Fraction(x, SCALE) for x in numerators)
                    channel[pair] = first + (1 - sum(first),)
            else:
                channel = {pair: (Fraction(1), Fraction(0), Fraction(0))
                           for pair in joint}
            extended = {}
            for pair, probabilities in channel.items():
                normalized(probabilities, "Filter channel")
                for label, probability in enumerate(probabilities):
                    extended[pair + (label,)] = joint[pair] * probability
            require(marginal(extended, (0, 1)) == joint, "Filter extension mismatch")
            joints.append(joint)
            filters.append(extended)
        all_joints.append(joints)
        all_types.append(types)
        all_filters.append(filters)
    return all_joints, all_types, all_filters


def interval(value):
    value = Fraction(value)
    return iv.mpf(value.numerator) / value.denominator


def entropy(values):
    """Directed interval enclosure of the entropy of an exact rational law."""
    return -sum((interval(p) * iv.log(interval(p), 2) for p in values if p), interval(0))


def interval_max(values):
    values = list(values)
    return iv.mpf([max(value.a for value in values), max(value.b for value in values)])


def information_features(joints, types, filters):
    """Average the classwise information quantities before any cost maximum."""
    node_count, internal_count = len(types[0]), len(joints[0])
    h = [interval(0) for _ in range(node_count)]
    f = [interval(0) for _ in range(internal_count)]
    features = [(interval(0),) * 3 for _ in range(internal_count)]
    for coordinate_class in range(2):
        for node, distribution in enumerate(types[coordinate_class]):
            h[node] += entropy(distribution.values()) / 2
        for node, joint in enumerate(joints[coordinate_class]):
            left, right = types[coordinate_class][2 * node + 1:2 * node + 3]
            joint_entropy = entropy(joint.values())
            f[node] += (joint_entropy - entropy(left.values()) - entropy(right.values())) / 2
            if node == 0 or node >= 7:
                # The channel is constant; its three information values are exactly zero.
                require(all(p == 0 for (x, y, z), p in filters[coordinate_class][node].items() if z != 0),
                        "Constant filter expected")
                continue
            extended = filters[coordinate_class][node]
            center_entropy = entropy(marginal(extended, (2,)).values())
            alpha = entropy(left.values()) + center_entropy - entropy(marginal(extended, (0, 2)).values())
            beta = entropy(right.values()) + center_entropy - entropy(marginal(extended, (1, 2)).values())
            kappa = joint_entropy + center_entropy - entropy(extended.values())
            features[node] = tuple(old + new / 2 for old, new in zip(features[node], (alpha, beta, kappa)))
    return h, f, features


# Affine expressions have exact rational coefficients in the entropy features.
# They establish the root and full-source identities before interval evaluation.
def affine_sum(*expressions):
    result = {}
    for expression in expressions:
        for key, coefficient in expression.items():
            result[key] = result.get(key, Fraction(0)) + coefficient
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def affine_scale(multiplier, expression):
    return {key: Fraction(multiplier) * coefficient for key, coefficient in expression.items()
            if multiplier * coefficient}


def constant(value):
    return {"one": Fraction(value)} if value else {}


def evaluate(expression, features):
    return sum((interval(coefficient) * features[key] for key, coefficient in expression.items()), interval(0))


def reconstruct_rates(h, f):
    """Construct every list/modulus rate, preserving exact symbolic identities."""
    require(set(CHILD_MODULI) == set(range(1, 15)), "Modulus table dimension")
    require(0 < SOURCE_FRACTION < 1 and SOURCE_FRACTION == Fraction(1, 3),
            "Full-source split must be one third")
    require(BOTTOM_HEIGHT == Fraction(1123259302947, 5349021095909),
            "Common stored bottom height mismatch")
    feature_values = {"one": interval(1), "h15": h[15]}
    feature_values.update({"f" + str(node): value for node, value in enumerate(f)})
    child = {node: constant(value) for node, value in CHILD_MODULI.items()}
    lists = {node: constant(BOTTOM_HEIGHT) for node in range(16, 31)}
    lists[15] = affine_sum({"h15": Fraction(1)}, affine_scale(-1, child[7]))
    child[0] = affine_sum(*lists.values(), constant(-1),
                          *[affine_scale(-1, child[node]) for node in range(1, 7)],
                          *[child[node] for node in range(7, 15)],
                          *[{"f" + str(node): Fraction(1)} for node in range(15)])
    moduli = {0: constant(1)}
    moduli.update({node: child[(node - 1) // 2] for node in range(1, 31)})
    for node in reversed(range(15)):
        lists[node] = affine_sum(lists[2 * node + 1], lists[2 * node + 2],
                                 affine_scale(-1, moduli[node]), child[node],
                                 {"f" + str(node): Fraction(1)})
    require(lists[0] == {}, "Root list recurrence is not identically zero")
    require(affine_sum({"h15": Fraction(1)}, affine_scale(-1, moduli[15]),
                       affine_scale(-1, lists[15])) == {},
            "Implicit bottom is not the complete modular type domain")
    require(all(lists[node] == constant(BOTTOM_HEIGHT) for node in range(16, 31)),
            "Stored bottom heights differ")
    return ({node: evaluate(value, feature_values) for node, value in lists.items()},
            {node: evaluate(value, feature_values) for node, value in moduli.items()},
            {node: evaluate(value, feature_values) for node, value in child.items()})


def decimal_bound(value, places=15, upper=True):
    """Format an endpoint with decimal rounding away from the interval."""
    sign, mantissa, exponent, _ = value._mpi_[1 if upper else 0]
    endpoint = Fraction((-1 if sign else 1) * mantissa) * Fraction(2) ** exponent
    scaled = endpoint * 10**places
    integer = (-((-scaled.numerator) // scaled.denominator) if upper
               else scaled.numerator // scaled.denominator)
    prefix = "-" if integer < 0 else ""
    integer = abs(integer)
    return prefix + str(integer // 10**places) + "." + str(integer % 10**places).zfill(places)


def verify():
    """Check all finite constraints and complete resource bounds in the paper."""
    iv.dps = PRECISION
    require(DEPTH == 4 and ALPHABET == (-1, 0, 1, 2), "Architecture mismatch")
    require(TARGET == Fraction(2099947, 10000000), "Target mismatch")
    joints, types, filters = reconstruct_probabilities()
    h, f, filter_features = information_features(joints, types, filters)
    lists, moduli, child_moduli = reconstruct_rates(h, f)
    implicit = (0, 1, 3, 7, 15)
    stored = [node for node in range(1, 31) if node not in implicit]
    bottom_stored = [node for node in stored if node >= 15]
    internal_stored = [node for node in stored if node < 15]
    require(len(bottom_stored) == 15 and len(internal_stored) == 11,
            "Stored source or internal merge missing")
    target = interval(TARGET)
    positive_checks = {}

    def positive(name, value):
        require(value.a > 0, name + ": lower interval endpoint is not positive")
        positive_checks[name] = value

    # Structural zero rates are established symbolically in reconstruct_rates.
    for node in range(1, 31):
        positive("list rate " + str(node), lists[node])
        positive("modulus rate " + str(node), moduli[node])
        positive("modulus upper bound " + str(node), 1 - moduli[node])
        if node != 15:
            positive("type capacity " + str(node), h[node] - moduli[node] - lists[node])
    for node in range(15):
        positive("nested modulus " + str(node), moduli[node] - child_moduli[node])
        positive("useful pair information " + str(node), -f[node])
    for node in range(1, 7):
        alpha, beta, kappa = filter_features[node]
        positive("query information " + str(node), alpha)
        positive("data information " + str(node), beta)
        positive("query incidence " + str(node), kappa - alpha)
        positive("data incidence " + str(node), kappa - beta)
        positive("filtered pair information " + str(node), -f[node] + kappa - alpha - beta)

    # Full right-source storage and every one of the fifteen stored sources.
    source = interval(SOURCE_FRACTION) * h[15]
    positive("right source domain", source)
    positive("complementary source domain", (1 - interval(SOURCE_FRACTION)) * h[15])
    static_costs = {"implicit right-source construction": source}
    for node in stored:
        static_costs["raw list " + str(node)] = lists[node]
    for node in bottom_stored:
        support = (lists[node] + moduli[node]) / 2
        positive("source half-domain " + str(node), h[node] / 2 - support)
        static_costs["stored source " + str(node)] = support

    # Stored merges report every output; both candidate visits and empty
    # query-address headers are counted. All index writes are charged fully.
    for node in internal_stored:
        left, right = lists[2 * node + 1], lists[2 * node + 2]
        alpha, beta, kappa = filter_features[node]
        delta = moduli[node] - child_moduli[node]
        query_tree = left + kappa - alpha
        candidate_tree = left + right - delta + kappa - alpha - beta
        report = (interval_max((query_tree, candidate_tree)) + lists[node]) / 2
        static_costs["stored merge " + str(node)] = interval_max((
            left, right, lists[node], kappa, right + kappa - beta, report))
    for node in implicit[:-1]:
        alpha, beta, kappa = filter_features[node]
        static_costs["implicit family " + str(node)] = kappa
        static_costs["implicit right index " + str(node)] = lists[2 * node + 2] + kappa - beta

    # Bottom-to-root nested amplitude amplification, using the actual
    # occurrence preparation recurrence rather than a free penalty budget.
    generator = interval_max((interval(0), moduli[15] - source)) / 2
    for node in reversed(implicit[:-1]):
        alpha, beta, kappa = filter_features[node]
        delta = moduli[node] - child_moduli[node]
        right = lists[2 * node + 2]
        address_penalty = kappa - alpha - right + delta - f[node]
        candidate_penalty = kappa - alpha - beta - f[node]
        positive("implicit pair information " + str(node), candidate_penalty)
        generator += interval_max((address_penalty, candidate_penalty)) / 2
    static = interval_max(static_costs.values())
    time = interval_max((static, generator))
    for name, cost in static_costs.items():
        positive("target gap: " + name, target - cost)
    positive("target gap: implicit preparation", target - generator)
    positive("target gap: total time", target - time)
    positive("target gap: total space", target - static)
    return {
        "time": time,
        "space": static,
        "implicit_preparation": generator,
        "target_margin": target - time,
        "minimum_positive_margin": interval_max((-value for value in positive_checks.values())) * -1,
        "checked_inequalities": len(positive_checks),
    }


def main():
    result = verify()
    print("PASS: QRACM exact probability laws and all finite resource inequalities")
    print("Interval precision:", PRECISION, "decimal digits")
    print("Time exponent upper bound:", decimal_bound(result["time"]))
    print("Total-space exponent upper bound:", decimal_bound(result["space"]))
    print("Implicit preparation upper bound:", decimal_bound(result["implicit_preparation"]))
    print("Target exponent: 0.2099947")
    print("Target margin lower bound:", decimal_bound(result["target_margin"], 18, upper=False))
    print("Minimum positive margin lower bound:",
          decimal_bound(result["minimum_positive_margin"], 18, upper=False))
    print("Checked strict inequalities:", result["checked_inequalities"])
    print("Algorithmic assumptions remain conditional, as stated in the paper.")


if __name__ == "__main__":
    main()
