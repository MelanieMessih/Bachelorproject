### Bachelor project 
- Bachelor programme: Chemistry (joint degree UvA/VU)
- Duration: 3 months (April 2023 - June 2023)

# Molecular fingerprints

#### remove_features(X_data, headers=[])
    This function takes a Dataframe or Matrix as input. When a DataFrame is given, "headers" can be left empty. 
    When a Matrix is given, headers should be spicified to make sure the column headers are known when working 
    with the data. The function removes columns with missing, nan or inf values and returns the cleaned up 
    Dataframe or Matrix and its headers. 
    
    Parameters:
    - X_data
    - headers=[]

#### read_data(csv_file_name)
    This function takes a CSV file name as input and reads it as a DataFrame. The independent and dependent 
    variables are selected and feature columns containing missing, nan or inf values are removed. The 
    function returns the DataFrame, independent variables, feature names, and dependent variables.
    
    Parameters:
    - csv_file_name


#### create_directory(path::AbstractString)  
    This function takes a path. If the path does not exist, the specified path of directories is created.  
    
    Parameters:  
    - path: ...

#### create_all_directories()  
    This function creates all directories needed for this project. Specifically, the function creates directories for 
    summaries, figures, fingerprints and models, and each option has three subdirectories: PaDEL, RDKit and 
    combined_fingerprints.  
    
    Parameters:  
    -

#### RDKit_fingerprints(dataset, fp_name, index_col_nr, smiles_col_nr, inchikeys_col_nr; nBits=nothing, radius=nothing)    
    This function takes a DataFrame, a fingerprint name that corresponds to one in the dictionary 'dict' below, 
    and the column numbers of indices, SMILES and inchikeys in the DataFrame. nBits and radius are optional arguments 
    that only apply for certain fingerprints, specified in 'dict'. The function generates RDKit fingerprints for each 
    molecule in the data set using their SMILES. The function returns a DataFrame with indices, inchikeys and the 
    generated fingerprints, and the feature names. 
    
    Parameters:  
    - dataset
    - fp_name
    - index_col_nr
    - smiles_col_nr
    - inchikeys_col_nr
    - nBits=nothing
    - radius=nothing
    
#### PaDEL_fingerprints(dataset, fp_name, index_col_nr, smiles_col_nr, inchikeys_col_nr)
    This function takes a DataFrame, a fingerprint name that corresponds to one in the descriptors.XML file, 
    and the column numbers of indices, SMILES and inchikeys in the DataFrame. It generates PaDEL fingerprints for 
    each molecule in the data set using their SMILES. The function returns a DataFrame with indices, inchikeys and 
    the generated fingerprints, and the feature names. 
    
    Parameters:
    - dataset
    - fp_name
    - index_col_nr
    - smiles_col_nr
    - inchikeys_col_nr
    
#### change_xml(wanted_fingerprint)
    This function takes a fingerprint name that corresponds to one in the descriptors.xml file and modifies 
    the XML file by changing the selected fingerprint value to true and the remaining fingerprint values to 
    false. The modified XML file is saved back to the descriptors.xml file. 
    
    Parameters:
    - wanted_fingerprint
    
#### create_train_test_split(total_df, y_data; start_col_X_data=1, train_size=0.9, random_state=42)
    This function takes a DataFrame that contains independent variables, and a vector of dependent variables. Optional 
    arguments are the column number of the DataFrame where the independent variables start, a float between 0 and 1
    that represents the training set size, and an integer representing the random state. The function creates a train
    test split using the chosen parameters. The train and test set of the total DataFrame, the dependent variables and 
    the independent variables are returned as matrices. 
    
    Parameters:
    - total_df
    - y_data
    - start_col_X_data=1
    - train_size=0.9
    - random_state=42

#### train_model(X_train, y_train, X_test, y_test; n_estimators=600, min_samples_leaf=4, max_features="sqrt", random_state=42)
    This function takes train matrices of independent and dependent variables, and test matrices of independent and 
    dependent variables. Optional arguments are the number of estimators, a minimum number of samples in each leaf, a 
    maximum number of features to consider when looking for the best split, and an integer representing the random state.
    The function creates and predicts a Random Forest Regressor using the chosen parameters and returns the trained model,
    the train and test score, and the train and test predictions.
    
    Parameters:
    - X_train
    - y_train
    - X_test
    - y_test
    - n_estimators=600
    - min_samples_leaf=4
    - max_features="sqrt"
    - random_state=42
    
#### cross_validation(model, X_train, y_train, test_score; cv=3)
    This function takes a model, matrices of independent and dependent variables, and a test score. The 
    optional argument is an integer that represents the number of folds for the cross validation. The 
    function performs a cross validation and returns an array of the cross validation score and R, which 
    is the median of the test score and the cross validation scores. 
    
    Parameters:
    - model
    - X_train
    - y_train
    - test_score
    - cv=3
    
#### important_features(model, headers; variance_explained=85) 
    This function takes a model and an array of variable names. The optional argument is an integer between 
    0 and 100 that represents the variance percentage you want to have explained by the important features. 
    The function selects the most important features based on the chosen variance explained. The function 
    returns the sum of the importances, an array of the headers, an array of the importances and the number of 
    the selected important features. 
    
    Parameters:
    - model
    - headers 
    - variance_explained=85
    
#### grid_best_params(param_grid, X_train, y_train)
    This function takes a dictionary with "min_samples_leaf", "n_estimators" and "max_features" as keys and 
    corresponding vectors as values, an array with X data and a vector of y data. The function performs a 
    5-fold (default) 3D grid search and returns the best parameters for min_samples_leaf, n_estimators and 
    max_features. 
    
    Parameters:
    - param_grid
    - X_train
    - y_train
    
#### plot_training_testing(fp_name, software_name, train_predictions, y_train, test_predictions, y_test; x_label="Measured LC50 (log(mg/L))", y_label="Predicted LC50 (log(mg/L))")
    This function takes a fingerprint name, a software name, and arrays of train predictions, train data, test 
    predictions and test data. The optional arguments are labels for the x and y axes. The function creates a
    scatter plot of predicted values vs. measured values for the training and testing set, together with the line 
    y = x, and saves the figure. 
    
    Parameters:
    - fp_name
    - software_name 
    - train_predictions 
    - y_train 
    - test_predictions
    - y_test
    - x_label="Measured LC50 (log(mg/L))"
    - y_label="Predicted LC50 (log(mg/L))"
   
#### plot_residuals(fp_name, software_name, residuals_train, y_train, residuals_test, y_test; x_label="Measured LC50 (log(mg/L))", y_label="Residuals")
    This function takes a fingerprint name, a software name, and arrays of train residuals, train data, test 
    residuals and test data. The optional arguments are labels for the x and y axes. The function creates a
    scatter plot of residuals vs. measured values for the training and testing set, together with the horizontal 
    lines y1 = 1 and y2 = -1, and saves the figure. 
    
    Parameters:
    - fp_name
    - software_name
    - residuals_train
    - y_train
    - residuals_test
    - y_test
    - x_label="Measured LC50 (log(mg/L))"
    - y_label="Residuals"
    
#### plot_distribution_residuals(fp_name, software_name, residuals; x_label="Residuals", y_label="Frequency")
    This function takes a fingerprint name, a software name, and an array of residuals. The optional arguments are 
    labels for the x and y axes. The function creates a histogram that represents the distribution of the residuals, 
    together with the horizontal line y = 0, and saves the figure. 
    
    Parameters:
    - fp_name
    - software_name
    - residuals
    - x_label="Residuals"
    - y_label="Frequency"
    
#### plot_important_features(fp_name, software_name, important_features, importances; x_label="", y_label="Percentage importance")
    This function takes a fingerprint name, a software name, and arrays of important features and the 
    corresponding importances. The optional arguments are labels for the x and y axes. The function creates
    a bar plot of the important features and their importance percentages, and saves the figure. 
    
    Parameters:
    - fp_name
    - software_name
    - important_features
    - importances
    - x_label=""
    - y_label="Percentage importance"
    
#### function create_plots_of_summary(summary_opt_filename, software_name)
    This function takes a filename of a summary of optimized fingerprint models and the corresponding software 
    name. The function creates the following plots: 1) a scatter plot of predicted values vs. measured values 
    for the training and testing set, together with the line y = x; 2) a scatter plot of residuals vs. measured 
    values for the training and testing set, together with the horizontal lines y1 = 1 and y2 = -1; 3) a histogram 
    that represents the distribution of the residuals, together with the horizontal line y = 0; 4) a bar plot of 
    the important features and their importance percentages, and saves the figures. 
    
    Parameters:
    - summary_opt_filename
    - software_name
    
#### function remove_low_scores(summary_opt; limit_train_score=0.8, limit_test_score=0.5)
    This function takes a CSV.File object that contains a summary and two optional arguments which are floats 
    between 0 and 1 that represent the minimum train and test score for the selection of the fingerprints. The 
    function removes the fingerprints that do not meet the minimum train and test scores and returns the new 
    summary as a DataFrame.
    
    Parameters:
    - summary_opt: CSV.File object.
    - limit_train_score=0.8: Float between 0 and 1.
    - limit_test_score=0.5: Float between 0 and 1.
    
#### summary_best_scores(summary_opt_filename, software_name; limit_train_score=0.8, limit_test_score=0.5)
    This function takes a filename of a summary of optimized fingerprint models and the corresponding software 
    name. The optional arguments are two floats between 0 and 1, which represent the minimum train and test score 
    for the selection of the fingerprints. The function saves a subselection of the total summary containing
    scores, and creates a similar summary of scores with only the fingerprints that meet the chosen minimum train 
    and test scores. 
    
    Parameters:
    - summary_opt_filename
    - software_name
    - limit_train_score=0.8 
    - limit_test_score=0.5
    
#### create_combined_fp(filename_scores_PaDEL, filename_scores_RDK, filename_combined_fp)
    This function takes the existing filenames of selected PaDEL and RDKit summaries of scores, and a new filename 
    for the combined fingerprint. The function combines the most important features of the PaDEL and RDKit 
    fingerprints and stores them in a new CSV file.
    
    Parameters:
    - filename_scores_PaDEL
    - filename_scores_RDK
    - filename_combined_fp
    
#### train_optimize_fp(fp_filename, software_name, y_data, new_model_name, param_grid; variance_explained=85)
    This function takes a fingerprint filename, a software name, an array of y data, a new model name, and a dictionary 
    with "min_samples_leaf", "n_estimators" and "max_features" as keys and corresponding vectors as values. The optional 
    argument is an integer between 0 and 100 that represents the variance percentage you want to have explained by the 
    important features. The function creates a train test split, optimizes and trains a random forest regressor model, 
    performs 3-fold cross validation, collects the important features that explain the chosen variance. The function 
    saves a summary of the optimized fingerprint model and a summary of the scores as CSV files and the optimized model 
    as a JOBLIB file. 
    
    Parameters:
    - fp_filename
    - software_name
    - y_data
    - new_model_name
    - param_grid
    - variance_explained=85
    
#### str_to_floatVec(string)
    This function takes a string of substrings separated by semicolons and parses the substrings to separated parts of 
    type Float64. The function returns the resulting float vector. 
    
    Parameters:
    - string: String of substrings separated by semicolons. Substrings should be able to be transformed to floats. 
    
    Returns:
    - Vector with elements of type Float64.
    
#### str_to_strVec(string)
    This function takes a string of substrings separated by semicolons and splits the substrings to separated parts of
    type String. The function returns the resulting string vector. 
    
    Parameters:
    - string: String of substrings separated by semicolons.
    
    Returns:
    - Vector with elements of type String.

#### compute_residuals(y, predictions)
    This function takes a vector of y data and a vector of predictions and computes the residuals. The function returns 
    the resulting vector of residuals. 
    
    Parameters:
    - y: Vector of type Float or Int.
    - predictions: Vector of type Float or Int.
    
    Returns:
    - Vector of type Float or Int.

#### write_csv(csv_filename, data) = CSV.write(csv_filename, data)
    This function takes a CSV filename and a collection of type Dictionary, DataFrame or Table. The function saves the 
    given data as a CSV file with the chosen csv_filename. 
    
    Parameters:
    - csv_filename: String (ending with ".csv").
    - data: Dictionary, DataFrame or Table.

#### save_model(model, model_name) = jl.dump(model, model_name)

    This function takes 

#### read_csv(csv_filename) = CSV.File(csv_filename)

    Parameters:
    - csv_filename: String (ending with ".csv").

    Returns:
    - CSV.File object.





Het inroosteren van lessen is een ingewikkeld probleem. In deze case moet een weekrooster gemaakt worden voor een vakkenlijst op Science Park. 

Hiervoor moeten 131 activiteiten ingepland worden. Dit kunnen hoorcolleges, werkcolleges en practica zijn.
- Een activiteit duurt 2 uur (= tijdslot)
- Maximale groepsgrootte bij werkcolleges en practica

Verder zijn er 7 zalen waarin de activiteiten kunnen plaatsvinden.
- Alle zalen zijn voor alle soorten activiteiten geschikt
- Capaciteit verschilt per zaal

Elk van de vakken kan worden ingedeeld in een van de 145 tijdsloten. Dit zijn periodes van 2 uur.
- Elke zaal heeft vier tijdsloten overdag (9-11u, 11-13u, 13-15u, 15-17u)
- Grootste zaal heeft ook een avondslot (17-19u)

We hebben te maken met 609 Studenten.
- Elke student volgt maximaal 5 vakken


### Constraints

De hard constraints van onze case zijn als volgt:
- Alle activiteiten moeten worden ingeroosterd
- Maximaal één activiteit per tijdslot per zaal inplannen
- Student mag maximaal twee tussenuren hebben
- Houden aan de maximumgrootte van werkcolleges en practica
- Zo min mogelijk werkcollege- en practicumgroepen

Naast het genereren van een geldige oplossing wordt er gekeken naar de kwaliteit van het rooster. Er wordt een aantal maluspunten toegekend bij het overtreden van de volgende soft constraints:
- Studenten met tussenuren (α = # keren een tussenuur per dag per student)
- Studenten met 2 tussenuren (β = # keren twee tussenuren per dag per student)
- Studenten met twee activiteiten in hetzelfde tijdslot (γ = # lessen die overlappen per student)
- Gebruik van avondslot (δ = # gebruikte avondsloten per lokaal)
- Studenten die niet in het lokaal passen (ε = # studenten die niet in lokaal passen)

### Goal

De kwaliteit van het rooster wordt gemeten aan de hand van de volgende objective function, die geminimaliseerd moet worden:

- f(α, β, γ, δ, ε) = α + 3⋅β + γ + 5⋅δ + ε

### Requirements

This codebase is fully written in Julia [3.9.13]. The packages needed to sucessfully run the code are provided below:

```
using Pkg, CSV, DataFrames, PyCall, Conda, ScikitLearn, Statistics, Plots, Tables, Plots.PlotMeasures, LightXML, LinearAlgebra, ProgressBars, OrderedCollections, Base.Filesystem
```

Other packages are installed using pyimport and are given in the "import.jl" file.

### Use

An example of how to run the function is provided below:

```
include("final_function.jl")
create_best_fingerprint(2, index_col_nr=1, inchikeys_col_nr=4)
```

### Files needed

An overview of the files needed in the main directory is provided below:

```
descriptors.xml
toxicity_data_fish_desc.csv
final_functions.jl
final_imports.jl
final_function.jl
```

Het bestand geeft aan hoe verschillende functies en algoritmes gebruikt kunnen worden.

Indien er met behulp van de instructies in main.py voor gekozen is om een yaml file van het rooster te maken, is het ook mogelijk om een visualisatie van het rooster per lokaal op te vragen. Dit kan door het aanroepen van:

```
pdfschedule --font Courier --color data/room{room_name}.yaml figures/room{room_name}.pdf
```

Hierbij kan voor {room_name} een van de volgende lokalen ingevuld worden:
- A1.04, A1.06, A1.08, A1.10, B0.201, C0.110, C1.112

Een plot waarin het verloop van een van de Hillclimbers wordt getoond kan worden aangeroepen met de make_plot() functie.
Een histogram van de score van een aantal Random of Greedy roosters kan worden aangeroepen met de make_histogram() functie.
Allebei de functies kunnen worden gerund door het aanroepen van:
```
python experiments.py
```

Onderaan het bovengenoemde bestand kunnen de verschillende optionele arguments handmatig aangepast worden.

### Structure

The following list describes the directories and files from this project, and where these can be found:

- **/code**: contains the code of this project
  - **/code/algorithms**: contains the code for the algorithms
  - **/code/classes**: contains the classes needed for this case
  - **/code/experiments.py**: contains the code for performing the experiments
- **/data**: contains the different files needed to fill and visualize the timetables
- **/figures**: contains the different results of the experiments and the visualisation of the best timetable for each classroom
- **/presentation**: contains the final presentation of our project

## Authors
- Melanie Messih

## Supervisors
- dr. Saer Samanipour
- Viktoriia Turkina MSc

## Research group, research institute and university
Analytical Chemistry Group, Van 't Hoff Institute for Molecular Sciences (HIMS), University of Amsterdam (UvA)
