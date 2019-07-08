from joblib import dump, load
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
from sklearn.metrics import f1_score, recall_score, precision_score, accuracy_score
from sklearn import metrics 
from sklearn.metrics import precision_recall_curve
from inspect import signature

# given a iteraction type that is an int, this will generate the 
# dictionary to map the integer outputs to real class labels 
def get_reverse_label_map(selected_interactions):
    elems = sorted(set(selected_interactions.values()))
    return {i:elems[i] for i in range(len(elems))}


# given a drugpair, this function will generate the associated 
# smiles pair 
def get_smiles(drug_pair, smiles_dict):
    #print(drug_pair) 
    drug1 = drug_pair[0]
    drug2 = drug_pair[1]
    if all(drug in list(smiles_dict.keys()) for drug in [drug1, drug2])== True:
    # print(drug1, drug2)
        smiles1 = smiles_dict[drug1]
        smiles2 = smiles_dict[drug2]
        smiles_pair =  (smiles1, smiles2)
    else:
        smiles_pair = None

    return smiles_pair

# This will generate a prediction label class for the drug pairs in 
# a validation set
def predict_val_set(model, drug_pair, fingerprints_dict):
    if drug_pair[0] not in fingerprints_dict.keys():
        print("fingerprints not available for %s" % drug_pair[0])
        return
    if drug_pair[1] not in fingerprints_dict.keys():
        print("fingerprints not available for %s" % drug_pair[1])
        return

    f1 = fingerprints_dict[drug_pair[0]]
    f2 = fingerprints_dict[drug_pair[1]]
    pair_fingerprint = np.concatenate((f1, f2))
    result = model.predict([pair_fingerprint])
    return result


# Plots confusion matrix for a multiclass classifier 
def plot_confusion_matrix(y_true, y_pred, classes,names,normalize=False,title=None, cmap=plt.cm.Blues):
    
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        #print("Normalized confusion matrix")
    #else:
        #print('Confusion matrix, without normalization')

    #print(cm)

    fig, ax = plt.subplots( figsize = (16, 16))
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=names, yticklabels=names,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    '''
    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    '''
    #fig.tight_layout()
    return ax

# reports accuracy, precision, recall and f1 score for a model
def generate_model_report(model, x_test, y_test):
    y_pred = model.predict(x_test)
    print("Accuracy: ", accuracy_score(y_test, y_pred))
    print("Precision: ", precision_score(y_test, y_pred, average = 'weighted'))
    print("Recall: ", recall_score(y_test, y_pred, average = 'weighted'))
    print("F1 Score: ", f1_score(y_test, y_pred, average = 'weighted'))


# for per class prevision recall computation, this function converts
# the multiclass predictios to a 2 class problem by taking class of interest as 1 
# and assuming all other classes to be 0. [For Probabilities]
def convert_to_2_class_proba(y_true, predict_proba, cls):
    new_ytrue = []
    new_ypred = []
    for i in range(len(y_true)):
        if y_true[i] == cls: 
            new_ytrue.append(1)
            new_ypred.append(predict_proba[i][cls])
        else:
            new_ytrue.append(0)
            new_ypred.append(predict_proba[i][cls])
    return new_ytrue, new_ypred

# for per class prevision recall computation, this function converts
# the multiclass predictios to a 2 class problem by taking class of interest as 1 
# and assuming all other classes to be 0. [For discrete predictions]
def convert_to_2_class(y_true, y_pred, cls):
    new_ytrue = []
    new_ypred = []
    for i in range(len(y_true)):
        if y_true[i] == cls and y_pred[i] == cls: 
            new_ytrue.append(1)
            new_ypred.append(1)
        elif y_true[i] == cls and y_pred[i] != cls:
            new_ytrue.append(1)
            new_ypred.append(0)
        elif y_true[i] != cls and y_pred[i] ==cls: 
            new_ytrue.append(0)
            new_ypred.append(1)
        elif y_true[i] != cls and y_pred[1] !=cls:
            new_ytrue.append(0)
            new_ypred.append(0)

    return new_ytrue, new_ypred

# Generates a per class report for the model
def generate_model_report_per_class(model, y_test, y_pred, classes):
    accuracy_per_class = {}
    precision_per_class = {}
    recall_per_class = {}
    f1score_per_class = {}
    for cls in classes:
        new_ytrue, new_ypred = convert_to_2_class(y_test, y_pred, cls)
        accuracy_per_class[cls] = accuracy_score(new_ytrue, new_ypred)
        precision_per_class[cls] = precision_score(new_ytrue, new_ypred)
        recall_per_class [cls] = recall_score(new_ytrue, new_ypred)
        f1score_per_class[cls] = f1_score(new_ytrue, new_ypred)
    return accuracy_per_class, precision_per_class, recall_per_class, f1score_per_class


# plots a pr curve for every class. 
def draw_pr_curves(y_true, predict_proba, classes):
    lines = []
    label_points =  classes #['increase antiplatelet activities', 'serum concentration metabolites reduced']
    plt.figure (figsize = (11,8))
    for i in range(len(classes)):
        new_ytrue, new_ypred = convert_to_2_class_proba(y_true, predict_proba, classes[i])
        #print(new_ytrue, new_)
        precision, recall, thresholds = precision_recall_curve(new_ytrue, new_ypred, 1)

        lines += plt.plot(recall, precision)
        #plt.step(recall, precision, color='b', alpha=0.2, where='post')
        #plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)
    plt.legend(classes, loc='lower center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision Recall curve by class')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
        #plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(average_precision))




