import xml.etree.ElementTree as ET
import xmlschema
import pickle

class Drug:
    def __init__(self):
        self.name = None
        self.description = None
        self.indication = None
        self.structure = None
        self.water_solubility = None
        self.logP = None # logP is a measure of lipophilicity
        self.logS = None
        self.polarizibility = None
        self.refractivity = None
        self.physio_charge = None
        self.num_hbond_donors = None
        self.num_hbond_acceptors = None
        self.strongest_acid_pka = None
        self.strongest_basic_pka = None
        self.psa = None # Polar surface area
        self.interactions = None

    def __repr__(self):
        return "Name:%s Structure:%s" % (self.name, self.structure)



def read_data(xmlfile):
    my_tree = ET.parse(xmlfile)
    root = my_tree.getroot()
    count = 0
    drug_list = []
    for drug in root:
        d = Drug()
        #go through the drug descriptions and extract relevant features
        for i in range(0, len(drug)): 
            #print(drug[i].tag, drug[i].text)
            if (drug[i].tag == "{http://www.drugbank.ca}name"):
                d.name = (drug[i].text).lower()
            if (drug[i].tag == "{http://www.drugbank.ca}description"):
                d.description = drug[i].text
            if (drug[i].tag == "{http://www.drugbank.ca}indication"):
                d.indication = drug[i].text    
                
            # Now retrieve all the structural features of the drug as well as other properties that can be 
            # calculated using structural info
            if (drug[i].tag ==  "{http://www.drugbank.ca}calculated-properties"):
                calculated_properties = drug[i]
                properties = {}
                for calc_property in calculated_properties:
                    #print(calc_property[0].text, calc_property[1].text)
                    if calc_property[0].text == 'logS':
                        d.logS = calc_property[1].text

                    elif calc_property[0].text == 'SMILES':
                        d.structure = calc_property[1].text

                    elif calc_property[0].text == 'Water Solubility':
                        d.water_solubility = calc_property[1].text
                        
                    elif calc_property[0].text == 'logP':
                        d.logP = calc_property[1].text
                        
                    elif calc_property[0].text == 'Refractivity':
                        d.refractivity = calc_property[1].text
                        
                    elif calc_property[0].text == 'Polarizability':
                        d.polarizibility = calc_property[1].text
                        
                    elif calc_property[0].text == 'Physiological Charge':
                        d.physio_charge = calc_property[1].text
                    
                    elif calc_property[0].text == 'H Bond Donor Count':
                        d.num_hbond_donors = calc_property[1].text
                        
                    elif calc_property[0].text == 'pKa (strongest acidic)':
                        d.strongest_acid_pka = calc_property[1].text
                    
                    elif calc_property[0].text == 'pKa (strongest basic)':
                        d.strongest_basic_pka = calc_property[1].text
                        
                    elif calc_property[0].text == 'H Bond Acceptor Count':
                        d.num_hbond_acceptors = calc_property[1].text
                    
                    elif calc_property[0].text == 'Polar Surface Area (PSA)':
                        d.PSA = calc_property[1].text
                        
            # now take interactions which will be labels in machine learning model
            if (drug[i].tag ==  "{http://www.drugbank.ca}drug-interactions"):
                all_drug_interactions = {}
                drug_interactions = drug[i]
                for drug_interaction in drug_interactions:
                    interacting_drug = (drug_interaction[1].text).lower()
                    description = drug_interaction[2].text
                    all_drug_interactions[interacting_drug] = description
                d.interactions = all_drug_interactions
                #print(d.interactions)
    
        drug_list.append(d) 
    return drug_list

# cleans up those drugs from the list that either have no interaction or no structure information. 
def cleanup_data (drug_list):
    cleaned_druglist = []
    for drug in drug_list: 
        if drug.interactions == None or drug.structure == None:
            continue
        else:
            cleaned_druglist.append(drug)        
    return cleaned_druglist


# for particular drug in the druglist, checks if structure is present or not.
def drug_structure_present(drug, drug_list):
    for item in drug_list:
        if item.name == drug  and item.structure is not None: 
            return True
        else: 
            return False

# 
def cleanup_interactions(drug_list):
    cleaned_list = []
    for drug in drug_list:
        interactions = drug.interactions
        new_interactions = {}
        for key in interactions:
            #print(key)
            if drug_structure_present(key, drug_list) == True:
                value = interactions.get(key)
                new_interactions[key] =  value
        drug.interactions = new_interactions
        cleaned_list.append(drug)
    return cleaned_list



