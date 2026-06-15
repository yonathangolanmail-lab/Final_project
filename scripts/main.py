import matplotlib.pyplot as plt

# 1. DATA CONFIGURATION

file= open("../data/proteins.fa", "r")
data= file.read()
file.close()


LOCATIONS = ["Nucleus", "Mitochondria", "ER", "Cell Membrane", "Secreted"]

PROTEIN_GROUPS = {
    "Cell Membrane": ["STE2_YEAST","CAN1_YEAST","ITR1_YEAST","MEP2_YEAST","KCH1_YEAST","BOR1_YEAST","PDR12_YEAST","SHO1_YEAST"],
    "Mitochondria": ["ODO1_YEAST","CCPR_YEAST","PREP_YEAST","MTF1_YEAST","ADH3_YEAST","EFTU_YEAST","ODPA_YEAST","IDH2_YEAST"],
    "ER": ["BIP_YEAST","PDI_YEAST","EUG1_YEAST","MPD1_YEAST","SIL1_YEAST","LHS1_YEAST","MNL1_YEAST","SCJ1_YEAST"],
    "Nucleus": ["RAP1_YEAST","RAD50_YEAST","RAD9_YEAST","SNF6_YEAST","MED15_YEAST","PRP4_YEAST","TAF13_YEAST","SPT2_YEAST"],
    "Secreted": ["SAG1_YEAST","HS150_YEAST","UTH1_YEAST","PRY1_YEAST","SUN4_YEAST","PIR3_YEAST","PPA5_YEAST","CHIT_YEAST"]
}


# 2. TRUTH LABELS

true_locations = {}

for location, protein_list in PROTEIN_GROUPS.items():
    for protein in protein_list:
        true_locations[protein] = location

# 3. LOAD FASTA DATA

records = data.split(">")
proteins= []
correct_predictions = 0
total_predictions = 0

confusion = {}
for true_loc in LOCATIONS:
    confusion[true_loc] = {}

    for pred_loc in LOCATIONS:
        confusion[true_loc][pred_loc] = 0

# amino acid property groups (extendable if we want to add new motifs)
HYDROPHOBIC = set("AILMFWVPG")
BASIC = set("KR")
MITO_POSITIVE = set("RKH")

location_correct = {}
location_total = {}


for location in LOCATIONS:
    location_correct[location] = 0
    location_total[location] = 0

for record in records[1:]:
    lines= record.split("\n")
    header= lines[0]
    sequence= ""
    for i in range(1, len(lines)):
        sequence= sequence+ lines[i].strip()
    proteins.append((header, sequence))


# 4. FEATURE DETECTION

# Detect N-terminal signal peptide
# returns True if signal peptide detected in first 30 amino acids
def has_signal_peptide(seq):
    first30 = seq[:30]
    count = 0
    for aa in first30:
        if aa in HYDROPHOBIC:
            count += 1
    return count >= 15

# Detect ER retention motif
def has_er_retention(seq):
    return (
        seq.endswith("KDEL") or
        seq.endswith("HDEL") or
        seq.endswith("RDEL")
    )

# Detect Nuclear Localization Signal (NLS)
def has_nls(seq):
    for i in range(len(seq)-4):
        window = seq[i:i+5]
        count = 0
        for aa in window:
            if aa in BASIC:
                count += 1
        if count >= 4:
            return True
    return False

# Detect mitochondrial targeting signal
def has_mito_signal(seq):
    first30 =seq[:30]
    count =0
    for aa in first30:
        if aa in MITO_POSITIVE:
            count += 1
    return count >= 5

# Detect transmembrane helix
def has_transmembrane_helix(seq):
    for i in range(len(seq)-20):
        window = seq[i:i+20]
        count = 0
        for aa in window:
            if aa in HYDROPHOBIC:
                count += 1
        if count >= 15:
            return True
    return False


# 5. CLASSIFIER

# rule-based classifier using biological motif hierarchy
def predict_location(seq):
    if has_er_retention(seq):
        return "ER"
    if has_mito_signal(seq):
        return "Mitochondria"
    if has_nls(seq):
        return "Nucleus"
    if has_signal_peptide(seq):
        return "Secreted"
    if has_transmembrane_helix(seq):
        return "Cell Membrane"
    return "Cell Membrane"


# 6. EVALUATION
# confusion matrix + precision + recall per class
results = []
for p in proteins:
    header = p[0]
    sequence = p[1]
    protein_name = header.split("|")[2].split()[0]
    true_location = true_locations[protein_name]
    location_total[true_location] += 1
    prediction = predict_location(sequence)
    confusion[true_location][prediction] += 1
    correct = (prediction == true_location)

    results.append([protein_name,true_location,prediction,correct])
    
    if correct:
        correct_predictions += 1
        location_correct[true_location] += 1
    total_predictions += 1
    
    print(protein_name)
    print("True: ", true_location)
    print("predicted: ", prediction)
    print("Correct: ", correct)
    print()

accuracy = correct_predictions / total_predictions

print("VALIDATION TABLE")
print()

for row in results:
    print(row)

print("Accuracy: ", accuracy)
print()

print("PRECISION AND RECALL")
print()

metrics = []

print("PRECISION AND RECALL\n")

for location in LOCATIONS:

    tp = confusion[location][location]

    fn = sum(confusion[location][pred] for pred in LOCATIONS if pred != location)

    fp = sum(confusion[true][location] for true in LOCATIONS if true != location)

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    metrics.append([location, round(precision, 3), round(recall, 3)])

    print(location)
    print("Precision:", round(precision, 3))
    print("Recall:", round(recall, 3))
    print()


print("CONFUSION MATRIX")
print()

print("\t", end="")

for loc in LOCATIONS:
    print(loc[:4], end="\t")

print()

for true_loc in LOCATIONS:

    print(true_loc[:4], end="\t")

    for pred_loc in LOCATIONS:
        print(confusion[true_loc][pred_loc], end="\t")

    print()


# 6. OUTPUT + VISUALIZATION

print("SUMMARY TABLE")
print()

print("Location\tPrecision\tRecall")

for row in metrics:
    print(row[0], "\t", row[1], "\t\t", row[2])

locations = []
recalls = []

for row in metrics:
    locations.append(row[0])
    recalls.append(row[2])

plt.bar(locations, recalls)

plt.title("Recall by Cellular Location")

plt.ylabel("Recall")

plt.show()