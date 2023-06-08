using CSV, DataFrames, PyCall, ScikitLearn, Statistics, Plots, Plots.PlotMeasures, LightXML, LinearAlgebra, ProgressBars, Base.Filesystem

train_test_split = pyimport("sklearn.model_selection").train_test_split
cross_val_score = pyimport("sklearn.model_selection").cross_val_score
skl = pyimport("sklearn.ensemble")
GridSearchCV =  pyimport("sklearn.model_selection").GridSearchCV

pd = pyimport("padelpy")

rdk = pyimport("rdkit.Chem")
AllChem = pyimport("rdkit.Chem.AllChem")
MACCSkeys = pyimport("rdkit.Chem.MACCSkeys")
Sheridan = pyimport("rdkit.Chem.AtomPairs.Sheridan")
Torsions = pyimport("rdkit.Chem.AtomPairs.Torsions")
Fingerprinter = pyimport("rdkit.Chem.EState.Fingerprinter")
FingerprintUtils = pyimport("rdkit.Chem.MolDb.FingerprintUtils")
rdReducedGraphs = pyimport("rdkit.Chem.rdReducedGraphs")
rdmolops = pyimport("rdkit.Chem.rdmolops")

jl = pyimport("joblib")