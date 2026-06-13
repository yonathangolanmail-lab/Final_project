file= open("../data/proteins.fa", "r")
data= file.read()
file.close()

true_locations = {

    "STE2_YEAST":"Cell Membrane",
    "CAN1_YEAST":"Cell Membrane",
    "ITR1_YEAST":"Cell Membrane",
    "MEP2_YEAST":"Cell Membrane",
    "KCH1_YEAST":"Cell Membrane",
    "BOR1_YEAST":"Cell Membrane",
    "PDR12_YEAST":"Cell Membrane",
    "SHO1_YEAST":"Cell Membrane",

    "ODO1_YEAST":"Mitochondria",
    "CCPR_YEAST":"Mitochondria",
    "PREP_YEAST":"Mitochondria",
    "MTF1_YEAST":"Mitochondria",
    "ADH3_YEAST":"Mitochondria",
    "EFTU_YEAST":"Mitochondria",
    "ODPA_YEAST":"Mitochondria",
    "IDH2_YEAST":"Mitochondria",

    "BIP_YEAST":"ER",
    "PDI_YEAST":"ER",
    "EUG1_YEAST":"ER",
    "MPD1_YEAST":"ER",
    "SIL1_YEAST":"ER",
    "LHS1_YEAST":"ER",
    "MNL1_YEAST":"ER",
    "SCJ1_YEAST":"ER",

    "RAP1_YEAST":"Nucleus",
    "RAD50_YEAST":"Nucleus",
    "RAD9_YEAST":"Nucleus",
    "SNF6_YEAST":"Nucleus",
    "MED15_YEAST":"Nucleus",
    "PRP4_YEAST":"Nucleus",
    "TAF13_YEAST":"Nucleus",
    "SPT2_YEAST":"Nucleus",

    "SAG1_YEAST":"Secreted",
    "HS150_YEAST":"Secreted",
    "UTH1_YEAST":"Secreted",
    "PRY1_YEAST":"Secreted",
    "SUN4_YEAST":"Secreted",
    "PIR3_YEAST":"Secreted",
    "PPA5_YEAST":"Secreted",
    "CHIT_YEAST":"Secreted"
}
records = data.split(">")
proteins= []
correct_predictions = 0
total_predictions = 0

location_correct = {}
location_total = {}

for location in [
    "Nucleus", "Mitochondria" , "ER", "Cell Membrane", "Secreted"
]:
    location_correct[location] = 0
    location_total[location] = 0

for record in records[1:]:
    lines= record.split("\n")
    header= lines[0]
    sequence= ""
    for i in range(1, len(lines)):
        sequence= sequence+ lines[i].strip()
    proteins.append([header,sequence])

def has_signal_peptide(seq):
    first30 = seq[:30]
    hydrophobic = "AILMFWVPG"
    count = 0
    for aa in first30:
        if aa in hydrophobic:
            count += 1
    return count >= 15

def has_er_retention(seq):
    return (
        seq.endswith("KDEL") or
        seq.endswith("HDEL") or
        seq.endswith("RDEL")
    )

def has_nls(seq):
    basic = "KR"
    for i in range(len(seq)-4):
        window = seq[i:i+5]
        count = 0
        for aa in window:
            if aa in "KR":
                count += 1
        if count >= 4:
            return True
    return False
    
def has_mito_signal(seq):
    first30 =seq[:30]
    positive = "RKH"
    count =0
    for aa in first30:
        if aa in positive:
            count += 1
    return count >= 5

def has_transmembrane_helix(seq):
    hydrophobic = "AILMFWVPG"
    for i in range(len(seq)-20):
        window = seq[i:i+20]
        count = 0
        for aa in window:
            if aa in hydrophobic:
                count += 1
        if count >= 15:
            return True
    return False

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



for p in proteins:
    header = p[0]
    sequence = p[1]
    protein_name = header.split("|")[2].split()[0]
    true_location = true_locations[protein_name]
    location_total[true_location] += 1
    prediction = predict_location(sequence)
    correct = prediction == true_location
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
print("Accuracy: ", accuracy)
print()

for location in location_total:
    success_rate = location_correct[location] / location_total[location]
    print(location)
    print("Correct:", location_correct[location])
    print("Total:", location_total[location])
    print("Success rate:", success_rate)
    print()

