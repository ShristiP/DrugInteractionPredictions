from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np 
import pandas as pd 
from sklearn.externals import joblib
import pickle

from .model import website_model, fingerprints, reverse_label_map


class PredictionResult:
    def __init__(self, smiles=None, drug=None, interaction=None, probability=None):
        self.smiles = smiles
        self.drug = drug
        self.interaction = interaction
        self.probability = probability
        self.err = None

    def __repr__(self):
        return "interaction:%s, probability:%.4f" % (self.interaction, self.probability)


def get_fingerprint_for_smiles(smiles):
    if smiles is None:
        return None

    rdkit_object = Chem.MolFromSmiles(smiles)
    if rdkit_object is None:
        return None

    fp_bitvec = AllChem.GetMorganFingerprintAsBitVect(rdkit_object, 2)
    if fp_bitvec is None:
        return None

    fp_nparr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp_bitvec, fp_nparr)
    return fp_nparr


def predict_interaction(smiles, drug_name):
    prediction = PredictionResult()
    prediction.smiles = smiles
    prediction.drug = drug_name

    # obtain name to fingerprint dictionary
    fingerprints.index = fingerprints.index.str.lower()
    name_to_fingerprints = fingerprints[['Molecular_Fingerprint']]
    name_to_fingerprints_dictionary = name_to_fingerprints.to_dict()['Molecular_Fingerprint']


    #name_to_smiles = fingerprints[['SMILES']]
    
    #smiles_to_fingerprints = fingerprints[['SMILES', 'Molecular_Fingerprint']]

    #smiles_to_fingerprints = smiles_to_fingerprints.set_index(['SMILES'])
    #smiles_to_fingerprints_dictionary = smiles_to_fingerprints.to_dict()['Molecular_Fingerprint']

    if smiles == "New Target SMILES" or smiles == "":
        return [prediction]
    if drug_name == "please enter drug name here" or drug_name == "":
        return [prediction]

    # compute fingerprint from provided smiles
    print("compute fingerprint from provided smiles")
    fp_smiles = get_fingerprint_for_smiles(smiles)
    if fp_smiles is None:
        print("fp_smiles is None")
        prediction.err = "could not get fingerprint for smiles:%s" % smiles
        return [prediction]
    
    # compute fingerprint for the provided drug
    print("compute fingerprint for the provided drug")
    fp_drug = None
    if drug_name in name_to_fingerprints_dictionary:
        fp_drug = name_to_fingerprints_dictionary[drug_name]
    if fp_drug is None:
        print("fp_drug is None")
        prediction.err = "could not get fingerprint for drug:%s" % drug_name
        return [prediction]
    
    # predict on the combined fingerprint
    pair_fingerprint = np.concatenate((fp_drug, fp_smiles))
    result = website_model.predict([pair_fingerprint])[0]
    probas = website_model.predict_proba([pair_fingerprint])[0]

    top_pred_indices = np.flip(np.argsort(probas))

    #best_idx = np.argmax(probas)
    assert(top_pred_indices[0] == result and "mismatch between probas and prediction")

    top_preds = []
    for i in range(5):
        prediction = PredictionResult()
        prediction.smiles = smiles
        prediction.drug = drug_name
        prediction.interaction = reverse_label_map[top_pred_indices[i]]
        prediction.probability = probas[top_pred_indices[i]]
        top_preds.append(prediction)
    return top_preds

    '''


    

    # set prediction object
    prediction.interaction = reverse_label_map[result]
    prediction.probability = probas[best_idx]

    print("Predicted interaction is: ",prediction)
    return prediction
    '''


if __name__ == "__main__":
    res = predict_interaction("CCCCOC(=O)C1=CC=C(C=C1)N", "altretamine")
    print (res)
