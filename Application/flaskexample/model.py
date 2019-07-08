from joblib import dump, load
import pickle

website_model = load('/Users/shristi/Documents/InsightDataProject/Models/RFmodel01.joblib')

infile = open('/Users/shristi/Documents/InsightDataProject/Data/ProcessedData/fingerprints.pkl','rb')
fingerprints = pickle.load(infile)

label_map = {'adverse effects increased': 0, 'bioavailability decreased': 1, 'decrease absorption': 2, 'decrease anticoagulant activities': 3, 'decrease antihypertensive activities': 4, 'decrease bronchodilatory activities': 5, 'decrease cardiotoxic activities': 6, 'decrease diuretic activities': 7, 'decrease excretion rate': 8, 'decrease sedative activities': 9, 'decrease the stimulatory activities': 10, 'increase analgesic activities': 11, 'increase anticholinergic activities': 12, 'increase anticoagulant activities': 13, 'increase antihypertensive activities': 14, 'increase antiplatelet activities': 15, 'increase arrhythmogenic activities': 16, 'increase atrioventricular blocking activities': 17, 'increase bradycardic activities': 18, 'increase cardiotoxic activities': 19, 'increase central nervous system depressant activities': 20, 'increase fluid retaining activities': 21, 'increase hyperkalemic activities': 22, 'increase hypertensive activities': 23, 'increase hypoglycemic activities': 24, 'increase hypokalemic activities': 25, 'increase hypotensive activities': 26, 'increase immunosuppressive activities': 27, 'increase nephrotoxic activities': 28, 'increase neuroexcitatory activities': 29, 'increase neuromuscular blocking activities': 30, 'increase respiratory depressant': 31, 'increase sedative activities': 32, 'increase serotonergic activities': 33, 'metabolism decreased': 34, 'metabolism increased': 35, 'serum concentration decreased': 36, 'serum concentration increased': 37, 'serum concentration metabolites reduced': 38, 'therapeutic efficacy decreased': 39}

reverse_label_map = {v:k for k, v in label_map.items()}
