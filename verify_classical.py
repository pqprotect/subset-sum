#!/usr/bin/env python3
"""Verify the classical filtered representation-tree exponent.

All exact probability parameters are embedded below. This program checks
arithmetic feasibility of the finite construction; the recursive-list,
survival and access assumptions remain the hypotheses stated in the paper.

Requires Python 3 and mpmath 1.3.x. Run: python verify_classical.py
"""

from fractions import Fraction

from mpmath import iv

ALPHABET = (-2, -1, 0, 1, 2, 3)
DEPTH = 7
LABELS = 5
SCALE = 10**15
TARGET = Fraction("0.26339")
RESERVE = Fraction(1, 10**6)
BALANCE_DEPTH = 4
LABEL_TRANSPOSE = (0, 2, 1, 4, 3)

# Root rows: pair, joint-mass numerator over 2*SCALE, then the first
# four conditional label numerators over SCALE. The fifth is a complement.
ROOT_ROWS = (((-2, 2), 4436919484560, 4185991151348, 1522221320, 995483161894952, 1849943847),
 ((-2, 3), 189526284474, 4779875043, 5574073857, 999981567502445, 4048850034),
 ((-1, 1), 199979517485737, 742082631365041, 39357241953, 225800305143587, 43920776658),
 ((-1, 2), 48959863253981, 190663387363448, 617573189, 799324807418168, 621753432),
 ((0, 0), 732804534963726, 939304602099249, 11399998602615, 2936958300032, 13753361681069),
 ((0, 1), 671650700752187, 932046389833331, 1235997834866, 29632332185830, 1446665340509),
 ((1, -1), 62386275473247, 320236507470620, 297318766782183, 381380267, 374188754394347),
 ((1, 0), 272491410622780, 803503538386278, 76790922978095, 91445586219, 94757128869695),
 ((2, -2), 392752592730, 35787095389, 432362558172471, 2452291752, 567595490136191),
 ((2, -1), 6701745176734, 15873367727329, 426593168103984, 892213042, 557221317526731),
 ((3, -2), 6753909844, 147878627, 460841342508963, 1343393851, 539156825687638))

# Nonfinal unordered pairs, grouped by their integer sum.
SPLIT_PAIRS = ((-2, 0), (-2, 1), (-2, 2), (-1, 1), (-2, 3), (-1, 2), (-1, 3), (0, 2), (0, 3))

# Conditional split numerators over SCALE, by class and levels 2 through 7.
SPLIT_NUMERATORS = (((184348484396233,
   8089894165553,
   235120790126,
   82750880059507,
   2047272955,
   7886194324375,
   170993628126,
   167062478292296,
   4063269623478),
  (121120930424827,
   10533326640505,
   487019824286,
   69086066070496,
   1940593551,
   10210274143683,
   59863914415,
   99541888393479,
   33112612306116),
  (185159516204042,
   16167525219316,
   378523397862,
   63959380189629,
   19244179428,
   15542934098376,
   2341546654727,
   199874164154793,
   125156271353542),
  (271159053097506,
   1805277370780,
   15186787446,
   21020937032971,
   26387183709,
   1604671330563,
   4405473489863,
   255475105167061,
   433543716976359),
  (208392903263455,
   420474149753,
   1300071640,
   5406214171379,
   19745758376,
   394893125553,
   8211460905132,
   198012583856420,
   573942504565147),
  (179351725466068,
   168184831207,
   1433447915,
   676366820940,
   60308248297,
   77331576386,
   20673942102366,
   102496615108023,
   510407926244898)),
 ((339168650035370,
   53601714093576,
   2839281605257,
   290037419965728,
   54696942416,
   50825543781179,
   2937717663537,
   307232369423112,
   59755465813742),
  (330533560288195,
   106237384522682,
   16266042691828,
   323959352011482,
   927653224695,
   106527789462302,
   15422725214866,
   321907643943434,
   63833940770786),
  (352329545828143,
   90809723969208,
   16169931481868,
   315709534752207,
   1645298615510,
   91284484066193,
   16586488212595,
   351800777859379,
   62213564319892),
  (325869206918941,
   52695009355090,
   6452005095252,
   306504345887341,
   551393156167,
   49660063675559,
   7876482573148,
   320644754899453,
   83912630297782),
  (232633230470882,
   15138366124501,
   543407845672,
   153905384356408,
   15032454076,
   14571622130379,
   767597545865,
   232909238017270,
   38759171703205),
  (4681660857586,
   49935567602,
   1262849741,
   32604874304564,
   3564343150,
   52349401713,
   724096106474,
   4942580547610,
   667048800728620)))

# Ordered lexicographic list of unordered child pairs.
FILTER_PAIRS = ((-2, 0), (-2, 1), (-2, 2), (-2, 3), (-1, -1), (-1, 0), (-1, 1), (-1, 2), (-1, 3), (0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2))

# Channel rows use four numerators over SCALE off the diagonal.
# A diagonal row (u,v) means probabilities
# (1-(u+v)/SCALE, u/(2*SCALE), u/(2*SCALE), v/(2*SCALE), v/(2*SCALE)).
# Unlisted levels have the constant label 0.
FILTER_NUMERATORS = {(0, 2): ((770267563189443, 7649628375, 78769521916127, 69474924653372),
          (271949909309274, 1039160198, 668889420442321, 27976989431361),
          (22804111446363, 1756385876, 970928215246990, 2616604124768),
          (2559847466438, 205802111, 997384035421485, 23912631608),
          (1007036686601, 151333459623206),
          (842355060340882, 24196194874, 7546939903166, 74696054388468),
          (735526833377093, 1245459953, 133082427721327, 65121551553847),
          (219714115615066, 1046913247, 726774175138020, 24186485244337),
          (10996509510089, 114528758, 988626807354867, 234115764454),
          (751046911671, 150792152205243),
          (842471114365643, 18504647393, 7522991454284, 74684516836751),
          (720707573253274, 11313571561, 112348917375957, 72920260719389),
          (995694245711240, 2039947461, 3578783568641, 665516414068),
          (759277585585, 151553461185732),
          (881358093529031, 373810300102, 5481971941995, 53337780291947)),
 (0, 3): ((787076914899163, 7750939483, 170409847303124, 18374663138994),
          (133057750075055, 739917401, 858826947392049, 3763996309999),
          (7349378999630, 1690439933, 992159682691750, 207214408027),
          (117516288536882, 243220364, 882399485239452, 52197134135),
          (1373694701768, 54753073667366),
          (937233665746585, 19760652211, 9168401046024, 26749003844322),
          (718525546758532, 943466585, 240300474528627, 20527305552107),
          (111019194773933, 713946913, 881685718890867, 3330983135848),
          (972611705604487, 459829396, 26976596757169, 303121022773),
          (543940105437, 53961262005834),
          (937207073227343, 16611757468, 9204723231413, 26749324449246),
          (693787949793102, 11911298097, 259952688877227, 19330120448845),
          (995424395176655, 13786252874, 3624834773282, 655454577926),
          (1168447516525, 55102699977640),
          (964306249902945, 1319308206184, 17992935139304, 7251189985949)),
 (0, 4): ((774423689664714, 6474825835, 182917843409383, 20622145143673),
          (89637722704743, 493785811, 904467103384879, 3070820357726),
          (24024627373229, 257102357, 975066589029183, 517726067468),
          (1345815281873, 8531144514, 998299319637197, 235772979779),
          (6473033849893, 50264423981160),
          (915251379890347, 50903063104, 18955439494165, 32309341053668),
          (436819011751122, 4810038988, 530646135521136, 15939365499295),
          (130792575746864, 2591662442, 862419017303855, 3606425705874),
          (6370600089168, 6525714278, 993055465655378, 334590676706),
          (685057649377, 67124046032499),
          (915671289246036, 44263671349, 18660375723781, 32433620516201),
          (827841266981875, 138211193338, 128101244631453, 21786372744971),
          (988487197805406, 321146452521, 8841721038754, 1660241893451),
          (5964668685058, 49166198957433),
          (255624149937372, 573186023135618, 127218101800800, 14631372217548)),
 (0, 5): ((444625048051973, 252952363467820, 252940245616349, 477362726357),
          (126083669124679, 72241258560477, 72302099742900, 663816502736),
          (20258000841989, 11847216944089, 11743450804358, 437928292441),
          (68889581763740, 42866674694803, 42846114640012, 954056880416),
          (532339316003668, 5563987645857),
          (443762680167366, 244944698751424, 244959141473364, 11468222153444),
          (357600943419621, 213212554716712, 213302597567084, 8219776327511),
          (107959284441198, 62381183859711, 62457853959452, 1281402888751),
          (189068271253998, 127043586411656, 127050763586880, 476981470924),
          (520939424240989, 21831154349561),
          (443897128543181, 244995203013755, 245129753128848, 10906544927828),
          (439901302235907, 247422171363208, 247473461616976, 4574856644017),
          (457984374520909, 263533513207592, 263711966800117, 349815164263),
          (531763731376697, 5731520944728),
          (358628269554793, 213962576207590, 213859536154271, 207732377692475)),
 (1, 2): ((269212504955722, 100862912088, 710293729229660, 236541098370),
          (34756962329080, 2416366786, 954790623134019, 12125675655),
          (3161615345088, 2014136607, 993565923569926, 1882572731),
          (1535018908870, 858302306, 997566046488478, 1434409898),
          (63297658820251, 16429552413085),
          (659099486669728, 3909571011890, 303222376495750, 2691535686255),
          (215455745287141, 111991261668, 754391797818945, 258450452439),
          (25586427354787, 2140598579, 962783080703595, 11017234098),
          (11040605899643, 1525596920, 985955532381594, 1086524355),
          (94173319980451, 25200367857326),
          (684197482770449, 3229574213565, 280855558540277, 2900240536550),
          (180610284639269, 119771557307, 794184002913461, 269394198032),
          (332999238613341, 7995169821, 662003492613697, 9914677325),
          (66137739927976, 22733028422084),
          (788397989552068, 2577384394763, 187393503133625, 2236629099369)),
 (1, 3): ((498846860825223, 151573884771, 493861326339143, 3676424846735),
          (98449917139406, 4103792994, 900087857124926, 741184362815),
          (13962553404977, 473754576, 985830866366774, 104949003108),
          (2779457410393, 1686000411, 997151818987466, 27250168606),
          (39760772555985, 14001821847165),
          (852570430581256, 2213217181466, 132638978036713, 6306859393972),
          (443109360129201, 144507835915, 550204646875649, 3286594095558),
          (93494396751318, 4170442841, 905119227075591, 697883007700),
          (20905168212682, 33072786836, 978707910031335, 162987909778),
          (37003501451855, 13997502222496),
          (852652583108671, 2093320527305, 132673293179680, 6303434220075),
          (449330834369087, 149670432837, 543877537663371, 3345242075158),
          (510424941702364, 43924103602681, 436488565324717, 3827611330014),
          (37243235519937, 14016200076291),
          (879527226438399, 2583615987688, 103773819149409, 7399305509410)),
 (1, 4): ((599341547867036, 654555276241, 276052438982021, 66082876982511),
          (209364058915787, 68157100876, 739938935265880, 27915427133328),
          (60474018662078, 6537618430, 924317241405343, 8878891610129),
          (4710282083987, 482496840, 993658372106996, 985146710373),
          (44799289607198, 150423661629919),
          (757572463263941, 3813191812173, 88224007418363, 77388598055636),
          (547163484738279, 692545309745, 345057695007943, 56540346204487),
          (227924681233066, 67801978873, 721336364510857, 28455446527605),
          (10185177566369, 5563755816, 985344645065981, 2508290567758),
          (33395586736711, 162130231800722),
          (760372868212727, 3932938431705, 86555325507596, 76009202820028),
          (606667212131393, 620228407042, 272739716148846, 66927018726488),
          (407195506653763, 550352188101, 510668042260786, 46968298466862),
          (47957650870000, 150911741371287),
          (786168506155031, 4424894472642, 43274107848818, 85707759696841)),
 (1, 5): ((305825151191802, 234050250027133, 233815328435505, 5107424959765),
          (109008773624268, 85993479234352, 85982035712105, 2538397682719),
          (33290946597559, 26758135627375, 26674873833377, 408897040039),
          (16434870592569, 12565411372845, 12582432180348, 330310518853),
          (581059463663195, 43744038358320),
          (345418315868605, 261727588294715, 261577142263363, 9936287535623),
          (290748883797630, 228715311891303, 228626460638281, 2693385899134),
          (122678327705488, 96803218332298, 96690703131241, 1979567113123),
          (68957505229979, 53097673710414, 53135027185839, 587527621085),
          (571401225372489, 60725386524796),
          (346912538696042, 263673782328158, 263633542174412, 11145764751637),
          (314924853582482, 240525820741672, 240398531669224, 5532820309156),
          (275926596203382, 209584043836074, 209710672312700, 486603351239),
          (577868262182676, 50235488618046),
          (375377270655821, 287681490149810, 288005527953084, 5989019328395)),
 (1, 6): ((317095442115168, 316303155230604, 316143890822106, 40363952402086),
          (199667781537555, 199907136362464, 199812045705732, 386885371425410),
          (64256712404568, 63692972446889, 63546127279755, 806585774957468),
          (104224619164491, 103946126817938, 103958415594829, 687329697479645),
          (653168724984252, 24775091598874),
          (324155269697516, 324122523364238, 324104904678561, 24355052236158),
          (312869083959722, 312383183545531, 312351918457252, 61031976488944),
          (210192312672633, 210003258415539, 209910751231226, 369083228937715),
          (254008545847064, 253266107782777, 253562804939282, 238619009023845),
          (657954753668537, 14953814348979),
          (324228465954632, 324291949190265, 324368345159805, 23568577876951),
          (320327381970373, 319342839864319, 319238972405577, 40004766898518),
          (327119212459756, 326525539697524, 326423259068900, 19330136567369),
          (652631068078631, 26300680446067),
          (331067772743666, 330508161073329, 330328565320810, 4516751407632))}

def require(condition, message):
    """Raise on a failed check, including when Python runs with -O."""
    if not condition:
        raise ValueError(message)


def marginal(distribution, coordinates):
    result = {}
    for symbols, probability in distribution.items():
        key = tuple(symbols[i] for i in coordinates)
        result[key] = result.get(key, Fraction(0)) + probability
    return result


def check_distribution(distribution, name):
    require(bool(distribution), name + ": empty distribution")
    require(min(distribution.values()) >= 0, name + ": negative mass")
    require(sum(distribution.values()) == 1, name + ": total mass is not one")


def label_row(numerators, diagonal=False):
    if diagonal:
        require(len(numerators) == 2, "Invalid diagonal channel row")
        u, v = (Fraction(n, 2 * SCALE) for n in numerators)
        row = (1 - 2 * u - 2 * v, u, u, v, v)
    else:
        require(len(numerators) == LABELS - 1, "Invalid channel row")
        first = tuple(Fraction(n, SCALE) for n in numerators)
        row = first + (1 - sum(first),)
    require(min(row) >= 0 and sum(row) == 1, "Invalid label probabilities")
    return row


def extend(joint, channels):
    extension = {}
    for (a, b), mass in joint.items():
        row = channels[a, b]
        for z, probability in enumerate(row):
            extension[a, b, z] = mass * probability
    check_distribution(extension, "Filter extension")
    require(marginal(extension, (0, 1)) == joint, "Incorrect extension marginal")
    return extension


def reconstruct_laws():
    """Reconstruct and check the exact rational laws from the integer tables."""
    require(DEPTH == 7 and LABELS == 5, "Incorrect architecture")
    require(SCALE == 10**15 and BALANCE_DEPTH == 4, "Incorrect exact scales")
    require(RESERVE == Fraction(1, 10**6), "Incorrect capacity reserve")
    require(len(SPLIT_NUMERATORS) == 2, "Two coordinate classes are required")
    all_pairs = tuple((a, b) for a in ALPHABET for b in ALPHABET
                      if a <= b and a + b in ALPHABET)
    require(FILTER_PAIRS == all_pairs, "Incorrect subtree pair support")
    groups = {}
    for pair in all_pairs:
        groups.setdefault(sum(pair), []).append(pair)
    nonfinal = {pair for group in groups.values() for pair in group[:-1]}
    require(set(SPLIT_PAIRS) == nonfinal, "Incorrect nonfinal split support")

    root = {}
    root_channels = {}
    for pair, mass, *numerators in ROOT_ROWS:
        require(pair not in root, "Repeated root pair")
        require(pair[0] in ALPHABET and pair[1] in ALPHABET,
                "Root symbol outside alphabet")
        require(sum(pair) in (0, 1), "Nonbinary root sum")
        root[pair] = Fraction(mass, 2 * SCALE)
        root_channels[pair] = label_row(numerators)
    check_distribution(root, "Root law")
    for s in (0, 1):
        require(sum(p for (a, b), p in root.items() if a + b == s)
                == Fraction(1, 2), "Root sum is not balanced")
    root_extension = extend(root, root_channels)
    root_marginals = [marginal(root, (side,)) for side in range(2)]
    profiles = []

    for side in range(2):
        require(len(SPLIT_NUMERATORS[side]) == DEPTH - 1,
                "Incorrect number of subtree levels")
        parent = root_marginals[side]
        levels = []
        for level, numerators in enumerate(SPLIT_NUMERATORS[side], 2):
            require(len(numerators) == len(SPLIT_PAIRS), "Incomplete split row")
            given = dict(zip(SPLIT_PAIRS,
                             (Fraction(n, SCALE) for n in numerators)))
            unordered = {}
            for s, pairs in groups.items():
                probabilities = [given[pair] for pair in pairs[:-1]]
                probabilities.append(1 - sum(probabilities))
                require(min(probabilities) >= 0, "Negative split probability")
                for pair, probability in zip(pairs, probabilities):
                    unordered[pair] = parent.get((s,), Fraction(0)) * probability

            channel_parameters = FILTER_NUMERATORS.get((side, level))
            constant = channel_parameters is None
            if not constant:
                require(len(channel_parameters) == len(all_pairs),
                        "Incomplete filter channel")
            joint = {}
            channels = {}
            for i, pair in enumerate(all_pairs):
                a, b = pair
                row = ((Fraction(1),) + (Fraction(0),) * (LABELS - 1)
                       if constant else
                       label_row(channel_parameters[i], diagonal=(a == b)))
                if a == b:
                    joint[a, b] = unordered[pair]
                    channels[a, b] = row
                    require(all(row[z] == row[LABEL_TRANSPOSE[z]]
                                for z in range(LABELS)),
                            "Diagonal channel violates label symmetry")
                else:
                    joint[a, b] = joint[b, a] = unordered[pair] / 2
                    channels[a, b] = row
                    channels[b, a] = tuple(row[z] for z in LABEL_TRANSPOSE)
            check_distribution(joint, "Subtree joint law")
            child = marginal(joint, (0,))
            require(child == marginal(joint, (1,)), "Unequal child marginals")
            summed = {s: sum(p for (a, b), p in joint.items() if a + b == s)
                      for s in ALPHABET}
            require(all(summed[s] == parent.get((s,), Fraction(0))
                        for s in ALPHABET), "Incorrect parent sum law")
            extension = extend(joint, channels)
            require(all(extension[a, b, z]
                        == extension[b, a, LABEL_TRANSPOSE[z]]
                        for a, b, z in extension), "Incorrect channel orientation")
            levels.append((joint, extension, child))
            parent = child
        profiles.append(levels)

    require(set(FILTER_NUMERATORS) == ({(0, j) for j in range(2, 6)}
                                     | {(1, j) for j in range(2, 7)}),
            "Incorrect constant-filter architecture")
    return root, root_extension, root_marginals, profiles


def interval_number(number):
    number = Fraction(number)
    return iv.mpf(number.numerator) / number.denominator


def entropy(distribution):
    return -sum((interval_number(p) * iv.log(interval_number(p), 2)
                 for p in distribution.values() if p), iv.mpf(0))


def information_rates(joint, extension):
    """Return I(X;Z), I(Y;Z), and I(XY;Z) with directed intervals."""
    labels = {z for (_, _, z), p in extension.items() if p}
    if len(labels) == 1:
        return (iv.mpf(0),) * 3
    hz = entropy(marginal(extension, (2,)))
    alpha = entropy(marginal(joint, (0,))) + hz - entropy(marginal(extension, (0, 2)))
    beta = entropy(marginal(joint, (1,))) + hz - entropy(marginal(extension, (1, 2)))
    kappa = entropy(joint) + hz - entropy(extension)
    return alpha, beta, kappa


def entropy_features(laws):
    """Average the two coordinate classes before any resource maximum."""
    root, extension, marginals, profiles = laws
    a, b, k = information_rates(root, extension)
    h = [(entropy(marginals[0]) + entropy(marginals[1])) / 2]
    f = [entropy(root) - entropy(marginals[0]) - entropy(marginals[1])]
    alpha, kappa = [(a + b) / 2], [k]
    for i in range(DEPTH - 1):
        class_features = []
        for side in range(2):
            joint, ext, child = profiles[side][i]
            a, b, k = information_rates(joint, ext)
            # The exact label-transposition symmetry checked above implies
            # I(X;Z)=I(Y;Z), including at a diagonal child pair.
            class_features.append((entropy(child), entropy(joint) - 2 * entropy(child), a, k))
        values = [(u + v) / 2 for u, v in zip(*class_features)]
        for target, value in zip((h, f, alpha, kappa), values):
            target.append(value)
    return [iv.mpf(1)] + h + f + alpha + kappa


class LinearForm:
    """An exact rational linear combination of the entropy features."""
    size = 1 + 4 * DEPTH

    def __init__(self, coefficients):
        self.coefficients = tuple(Fraction(x) for x in coefficients)
        require(len(self.coefficients) == self.size, "Invalid linear form")

    @classmethod
    def basis(cls, index):
        return cls(int(i == index) for i in range(cls.size))

    def __add__(self, other):
        return LinearForm(a + b for a, b in zip(self.coefficients, other.coefficients))

    def __neg__(self):
        return self.scale(-1)

    def __sub__(self, other):
        return self + (-other)

    def scale(self, factor):
        return LinearForm(Fraction(factor) * x for x in self.coefficients)

    def is_zero(self):
        return not any(self.coefficients)

    def evaluate(self, features):
        return sum((interval_number(c) * x for c, x in zip(self.coefficients, features)
                    if c), iv.mpf(0))


def check_resources(features):
    """Check every symbolic recurrence and every finite resource inequality."""
    zero = LinearForm([0] * LinearForm.size)
    one = LinearForm.basis(0)
    h = [None] + [LinearForm.basis(j) for j in range(1, DEPTH + 1)]
    f = [None] + [LinearForm.basis(DEPTH + j) for j in range(1, DEPTH + 1)]
    a = [None] + [LinearForm.basis(2 * DEPTH + j) for j in range(1, DEPTH + 1)]
    k = [None] + [LinearForm.basis(3 * DEPTH + j) for j in range(1, DEPTH + 1)]
    # Both deepest filters are exactly constant, as checked by reconstruction.
    a[DEPTH] = k[DEPTH] = zero
    incidence = [None] + [k[j] - a[j] for j in range(1, DEPTH + 1)]
    m = BALANCE_DEPTH
    t0 = (one - h[m] - sum(f[1:m + 1], zero)
          + sum(incidence[1:m + 1], zero)).scale(Fraction(1, m))
    time = t0 + one.scale(RESERVE / m)
    lists = [zero] + [time - incidence[j] for j in range(1, DEPTH + 1)]
    moduli = [one]
    for j in range(1, DEPTH + 1):
        moduli.append(moduli[-1] + lists[j - 1] - lists[j].scale(2) - f[j])
    source = (lists[DEPTH] + moduli[DEPTH]).scale(Fraction(1, 2))

    for j in range(1, DEPTH + 1):
        residual = lists[j - 1] - (lists[j].scale(2) - moduli[j - 1] + moduli[j] + f[j])
        require(residual.is_zero(), "List recurrence is not an exact identity")
        require((lists[j] + k[j] - a[j] - time).is_zero(),
                "Incidence-table costs are not exactly balanced")
    require((lists[DEPTH] - source.scale(2) + moduli[DEPTH]).is_zero(),
            "Bottom-source recurrence is not exact")
    require((h[m] - moduli[m] - lists[m] - one.scale(RESERVE)).is_zero(),
            "Capacity reserve is not exact")

    structural = {"positive time": time,
                  "root modulus upper": one - moduli[1],
                  "source domain": h[DEPTH].scale(Fraction(1, 2)) - source}
    costs = {"bottom source": source}
    for j in range(1, DEPTH + 1):
        structural[f"list {j}"] = lists[j]
        structural[f"modulus {j}"] = moduli[j]
        structural[f"capacity {j}"] = h[j] - moduli[j] - lists[j]
        if j < DEPTH:
            structural[f"nested modulus {j}"] = moduli[j] - moduli[j + 1]
        costs[f"input scan {j}"] = lists[j]
        costs[f"filter family {j}"] = k[j]
        costs[f"incidence table {j}"] = lists[j] + k[j] - a[j]
        costs[f"candidate visits {j}"] = (lists[j - 1] - f[j] + k[j] - a[j].scale(2))

    positive_margins = []
    for name, form in structural.items():
        value = form.evaluate(features)
        require(value.a > 0, name + ": nonpositive structural interval")
        positive_margins.append(value)
    active = []
    for name, cost in costs.items():
        gap = time - cost
        if gap.is_zero():
            active.append(name)
        else:
            value = gap.evaluate(features)
            require(value.a > 0, name + ": cost exceeds finite time budget")
            positive_margins.append(value)
        require((interval_number(TARGET) - cost.evaluate(features)).a > 0,
                name + ": cost reaches headline bound")
    require(all(f"incidence table {j}" in active for j in range(1, DEPTH + 1)),
            "A balancing identity was not recognized")
    finite = time.evaluate(features)
    margin = interval_number(TARGET) - finite
    require(margin.a > 0, "Finite time does not improve headline bound")
    return {"finite_exponent": finite,
            "headline_margin": margin,
            "minimum_positive_slack": min(positive_margins, key=lambda x: x.a),
            "exact_active_costs": active,
            "checked_costs": len(costs),
            "checked_structural_inequalities": len(structural)}


def verify():
    """Return the checked resource intervals; no files are read or written."""
    iv.dps = 100
    return check_resources(entropy_features(reconstruct_laws()))


def main():
    result = verify()
    print("PASS: exact rational type flows, channels and symbolic rate identities")
    print("PASS: all directed interval resource comparisons at 100 digits")
    print("Model: classical")
    print("Finite time and total-space exponent (approx.): "
          f"{float(result['finite_exponent'].a):.15f}")
    print("Certified time and total-space headline: < 0.26339")
    print("Headline margin (approx.): "
          f"{float(result['headline_margin'].a):.15g}")
    print("Minimum positive resource slack (approx.): "
          f"{float(result['minimum_positive_slack'].a):.15g}")
    print("Arithmetic feasibility does not prove the algorithmic assumptions.")


if __name__ == "__main__":
    main()
