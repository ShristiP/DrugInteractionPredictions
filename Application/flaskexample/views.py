from flask import render_template
from flaskexample import app
from flask import request 
#from sqalchemy import create_engine 
import pandas as pd 
import sklearn
from sklearn.externals import joblib
from joblib import dump, load
from .predict import predict_interaction
#from sqlalchemy_utils import database_exists, create_database
#simport psycopg2

@app.route('/')
@app.route('/home')
def home():
   return render_template('home.html')


@app.route('/about')
def about():
   return render_template('about.html')
        

@app.route('/pred_output')
def pred_output():
    smiles = request.args.get('SMILES').strip()
    drugname = request.args.get('Drugname').strip().lower()
    print("predicting interaction for smiles:%s, drug:%s" % (smiles, drugname))

    # obtain prediction and display results
    predictions = predict_interaction(smiles, drugname)
    assert(len(predictions) > 0 and "at least one result should be returned")
    #prediction = prediction[0]
    #print("top pred:%s" % prediction)

    if len(predictions) == 1:
        prediction = predictions[0]
        err_msg = prediction.err if prediction.err is not None else "Unknown error"
        return render_template("error.html",
            smiles=prediction.smiles,
            drug=prediction.drug,
            err_msg=err_msg)
    else:
        #prediction.probability = 0.9
        print(predictions)
        return render_template("result.html", top_preds=predictions, smiles=smiles, input_drug=drugname)