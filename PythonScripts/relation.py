
normalized_keywords = [
    "decrease anticoagulant activities",
    "increase anticoagulant activities",
    "serum concentration decreased",
    "serum concentration increased",
    "therapeutic efficacy decreased",
    "metabolism increased",
    "metabolism decreased",
    "adverse effects increased",
    "decrease cardiotoxic activities",
    "increase cardiotoxic activities",
    "increase qtcprolonging activities",
    "decrease excretion rate",
    "decrease antihypertensive activities",
    "increase neuroexcitatory activities",
    "decrease diuretic activities",
    "increase neutropenic activities",
    "increase immunosuppressive activities",
    "increase nephrotoxic activities",
    "decrease absorption",
    "increase hyperkalemic activities",
    "increase hypotensive activities",
    "increase serotonergic activities",
    "increase sedative activities",
    "increase hypoglycemic activities",
    "increase bradycardic activities",
    "increase central nervous system depressant activities",
    "increase antiplatelet activities",
    "decrease antineoplastic activities",
    "risk hypersensitivity reaction increased",
    "increase myopathic rhabdomyolysis activities",
    "increase hyponatremic activities",
    "increase anticholinergic activities",
    "increase neuromuscular blocking activities",
    "increase thrombogenic activities",
    "increase stimulatory activities",
    "increase hepatotoxic activities",
    "increase photosensitizing activities",
    "increase atrioventricular blocking activities",
    "increase adverse neuromuscular activities",
    "increase vasoconstricting activities",
    "serum concentration metabolites reduced",
    "increase antihypertensive activities",
    "increase absorption",
    "increase neurotoxic activities",
    "increase hypokalemic activities",
    "increase hypocalcemic activities",
    "increase hypertensive activities", 
    "increase constipating activities",
    "decrease the stimulatory activities",
    "bioavailability decreased", 
    "decrease sedative activities", 
    "increase respiratory depressant", 
    "decrease bronchodilatory activities", 
    "increase analgesic activities",
    "increase arrhythmogenic activities",
    "increase fluid retaining activities", 
    "increase tachycardic activities", 
    "decrease analgesic activities", 
    "decrease effectiveness",
    "decrease neuromuscular blocking activities", 
    "decrease effectiveness", 
    "increase vasopressor activities", 
    "increase vasodilatory activities",
    "decrease vasoconstricting activities",
    "increase ulcerogenic activities", 
    "increase hyperglycemic activities", 
    "may increase the excretion rate"
    
]

def inorder(relation, keywords):
    if len(keywords) == 0: return True
    idx = [relation.find(k) for k in keywords]
    if idx[0] == -1: return False
    for i in range(1, len(idx)):
        if idx[i] == -1 or idx[i] < idx[i-1]:
            return False
    return True

def normalizeRelation(relation):
    for n in normalized_keywords:
        if inorder(relation.relation, n.split()):
            return n
    return None


class Relation:
    def __init__(self):
        self.subject = None
        self.object = None
        self.orig_description = None
        self.relation = None
        self.normalizedRelation = None

    def __repr__(self):
        return "Subject:%s\nOriginal Description:%s\nDescription:%s\nNormalized:%s\nObject:%s" \
        % (self.subject, self.orig_description, self.relation, self.normalizedRelation, self.object)

def get_relation_class(drug_pair, relation): 
    parsed = parseEntry(drug_pair[0], drug_pair[1], relation)
    r = to_relation(parsed, relation)
    r.normalizedRelation = normalizeRelation(r)
    return r


def to_relation(entry, orig_description):
    assert(len(entry) == 4)
    relation = Relation()
    relation.orig_description = orig_description
    if entry[0] == '':
        relation.subject = entry[1]
        relation.relation = entry[2]
        relation.object = entry[3]
    else:
        relation.subject = entry[3]
        relation.relation = entry[0] + "__ " + entry[2]
        relation.object = entry[1]
    return relation

# given a drug pair, this function extracts the subject drug, the object drug and their relations
def parseEntry(e1, e2, description):
    e1 = e1.lower()
    e2 = e2.lower()
    description = description.lower()

    idx1 = description.find(e1)
    idx2 = description.find(e2)
    
    prefix = ""
    if min(idx1, idx2) != 0:
        prefix = description[:min(idx1, idx2)]
    
    if idx1 < idx2:
        relation = (prefix, e1, description[idx1 + len(e1): idx2].strip(), e2)
    else:
        relation = (prefix, e2, description[idx2 + len(e2): idx1].strip(), e1)
    return relation


def get_relations_for_drug(drug):
    relations = []
    for interacting_drug, description in drug.interactions.items():
        parsed_entry = parseEntry(drug.name.lower(), interacting_drug.lower(), description.lower())
        relation = to_relation(parsed_entry, description.lower())

        normalized = normalizeRelation(relation)
        relation.normalizedRelation = normalized

        relations.append(relation)
    return relations

def get_relations_for_drug_list(drug_list):
    relations = []
    for drug in drug_list:
        relations.extend(get_relations_for_drug(drug))
    return relations

def filter_unknown_relations(relations):
    return [relation for relation in relations if relation.normalizedRelation is not None]

if __name__ == "__main__":
    pass



