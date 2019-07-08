import rdkit as rd 
import pandas as pd 
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit.Chem import AllChem 
import numpy as np



# Constructs a molecule from smiles and returns the mol object 
# using the rdkit library
def get_profile(drug):
    #print(drug)
    drug1_mol = Chem.MolFromSmiles(drug)
    if drug1_mol is None: 
        return None
    return Chem.AddHs(drug1_mol)


# Computes a single number that corresponds to the similarity between
# two drugs 
def compute_similarity(profile_1, profile_2):
    fps_1 = AllChem.GetMorganFingerprint(profile_1, 2)
    fps_2 = AllChem.GetMorganFingerprint(profile_2, 2)

    if fps_1 is None or fps_2 is None: 
        return None
    
    score = DataStructs.DiceSimilarity(fps_1, fps_2)
    return score 


# Computes the tanimoto similarity vector by comparing a single drug to a
# all the drugs in the drug list. Can be used as a feature set


#given an input list of drugs, takes all possible combinations of drugs 
# in the list and computes a structure simiarity score for all pairs. 
def calculate_tanimoto_similarity(input_drug_list):
    structure_similarity_scores = {}

    for i, row in enumerate(drugbank_list.values):
        curr_drug = drugbank_list.index[i]
        structure_similarity_scores[curr_drug] = {}
        profile = get_profile(row[0])
        if profile is None:
            continue
        
        for j, row2 in enumerate(input_drug_list.values):
            curr_input_drug = input_drug_list.index[j]
            input_profile = get_profile(row2[0])
            if input_profile is None:
                continue
            
            score = compute_similarity(input_profile, profile)
            if score is not None: 
                structure_similarity_scores[curr_drug][curr_input_drug] = score  
    return structure_similarity_scores


# Featurize the SMILES into bit vectors using the molecular fingerprinting rules 
# Package: RDkit provides functionality
def featurize_smiles_mol_fingerprint_drug(drug):
    profile = get_profile(drug)
    arr = np.zeros((1,))
    if profile is not None:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(profile, 2)
        if fingerprint is not None:
            DataStructs.ConvertToNumpyArray(fingerprint,arr)
            return arr
        else:
            return None
    else:
        return None

#given a drug name, it finds the SMILES from the drug list 
def find_structure_in_list(drug, drug_structure_mapping):
    return drug_structure_mapping[drug] if drug in drug_structure_mapping else None

# creates a mapping between the drug name and drug structure
def get_structure_mapping(drug_list):
    mapping = {}
    for item in drug_list:
        mapping[item.name] = item.structure
    return mapping

# for every interaction label that is an integer value, there is a 
# interaction type. This function will genetate a label map going from 
# the real interaction type to the integer classes for modeling
def get_label_map(selected_labels):
    elems = sorted(set(selected_labels.values()))
    return {elems[i]:i for i in range(len(elems))}

# Important: given a drug_list, this function will featurize the 
# the structures and generate a feature matrix for ML models
def featurize_smiles_mol_fingerprint_all_pairs(labels, drug_list):
    structure_mapping = get_structure_mapping(drug_list)
    feature_dict = {}
    label_map = get_label_map(labels)
    X = []
    Y = []
    Z = []
    for drug_pair in labels:
        #print(drug_pair)
        drug_subject_structure = find_structure_in_list(drug_pair[0], structure_mapping)
        drug_object_structure = find_structure_in_list(drug_pair[1], structure_mapping)
        #print("Stucture_subject", drug_subject_structure)
        #print("Stucture_object", drug_object_structure)
        if (drug_subject_structure is None) or (drug_object_structure is None):
            continue

        feature_subject = None
        if drug_subject_structure in feature_dict:
            feature_subject = feature_dict[drug_subject_structure]
        else:
            feature_subject = featurize_smiles_mol_fingerprint_drug(drug_subject_structure)
            feature_dict[drug_subject_structure] = feature_subject


        #print(feature_subject)
        feature_object = None
        if drug_object_structure in feature_dict:
            feature_object = feature_dict[drug_object_structure]
        else:
            feature_object = featurize_smiles_mol_fingerprint_drug(drug_object_structure)
            feature_dict[drug_object_structure] = feature_object

        #print(feature_object)
        if (feature_subject is None) or (feature_object is None):
            continue
        features = np.concatenate((feature_subject, feature_object))

        interaction = labels[drug_pair]
        interaction_int = label_map[interaction]
        X.append(features)
        Y.append(interaction_int)
        Z.append(drug_pair)
    return X, Y, Z





