### Bachelor project 
- Bachelor programme: Chemistry (joint degree UvA/VU)
- Duration: 3 months (April 2023 - June 2023)

# Molecular fingerprints

SMALL DESCRIPTION OF THE PROJECT (GOAL AND WHAT THE FUNCTION DOES)

### Requirements

This codebase is fully written in Julia 1.6.7. The packages needed to sucessfully run the code are provided below:

```
using Pkg, CSV, DataFrames, PyCall, Conda, ScikitLearn, Statistics, Plots, Tables, Plots.PlotMeasures, LightXML, LinearAlgebra, ProgressBars, OrderedCollections, Base.Filesystem
```

To calculate molecular fingerprints, Java needs to be installed (https://www.java.com/en/download/manual.jsp). Other packages needed are installed using Conda and pyimport from PyCall.

The following list describes the folder structure that can be found on this page.
- **Bachelorproject**: all of the files in this folder should be included in the main directory when running the code.
  - **/final_function.jl**: Julia file needed to generate the best fingerprint based on the given dataset. 
  - **/final_imports.jl**: Julia file with all imports needed.
  - **/final_functions.jl**: Julia file with all functions needed.
  - **/toxicity_data_fish.csv**: CSV fish toxicity data set used for this project.
  - **/descriptors.xml**: XML file needed to generate PaDEL fingerprints from SMILES.

### Use

An example of how to run the function is provided below:

```
include("final_function.jl")
create_best_fingerprint(fish_toxicity_data, y_data, 2, index_col_nr=1, inchikeys_col_nr=4)
```

### Overview functions

- Included in final_function.jl:
#### create_best_fingerprint(dataset, y_data, smiles_col_nr; index_col_nr=nothing, inchikeys_col_nr=nothing, limit_train_score=0.8, limit_test_score=0.5, variance_explained=85)
    This function takes a DataFrame of the dataset that contains columns with X data, SMILES and optionally indices 
    and inchikeys, a vector of y data, and the column number of the DataFrame that contains SMILES. Optional
    parameters are the column numbers that contain indices and inchikeys, two floats between 0 and 1 that represent 
    the minimum train and test score for the selection of the fingerprints, and an integer between 0 and 100 that 
    represents the variance percentage you want to have explained by the important features. 
    
    The function generates PaDEL and RDKit fingerprints, creates and optimizes random forest regressor models with 
    them, selects the fingerprints of which the train and test score exceed the given limit_train_score and 
    limit_test_score, and combines the selected fingerprints into one fingerprint. This final fingerprint is used to
    create the final random forest model, which is trained and optimized. 
    
    During this process, summaries of the scores and results, the generated fingerprints, and the generated figures 
    are saved as CSV files and the created models are saved as JOBLIB files. 
    
    Parameters:
    - dataset: DataFrame that at least contains columns with X data, y data and SMILES of chemicals. 
    - y_data: Vector with elements of type Int or Float.
    - smiles_col_nr=nothing: Int.
    - index_col_nr=nothing: Int.
    - limit_train_score=0.8: Float between 0 and 1.
    - limit_test_score=0.5: Float between 0 and 1. 
    - variance_explained=85: Int between 0 and 100. 

- Included in final_functions.jl:
#### remove_features(X_data, headers=[])
    This function takes a DataFrame or Matrix as input. When a DataFrame is given, "headers" can be left empty. When 
    a Matrix is given, headers should be spicified to make sure the column headers are known when working with the 
    data. The function removes columns with missing, nan or inf values and returns the cleaned up Dataframe or Matrix 
    and its headers. 
    
    Parameters:
    - X_data: DataFrame or Matrix.
    - headers=[]: if X_data is a Matrix, headers should be specified as a Vector with elements of type Str.
    
    Returns:
    - X_data: DataFrame or Matrix.
    - headers: Vector with elements of type Str. 
    

#### read_data(csv_file_name)
    This function takes a CSV file name as input and reads it as a DataFrame. The independent and dependent variables 
    are selected and feature columns containing missing, nan or inf values are removed. The function returns the 
    DataFrame, independent variables, feature names, and dependent variables.
    
    Parameters:
    - csv_file_name: Str (ending with ".csv"); name of an existing CSV file.
    
    Returns:
    - data_name: DataFrame; contains dataset.
    - X_data: DataFrame with elements of type Float64; contains X data.
    - headers: Vector with elements of type Str. 
    - y_data: Vector with elements of type Float64; contains y data.


#### create_directory(path::AbstractString)  
    This function takes a path. If the path does not exist, the specified path of directories is created.  
    
    Parameters:  
    - path: Str.

#### create_all_directories()  
    This function creates all directories needed for this project. Specifically, the function creates directories for 
    summaries, figures, fingerprints and models, and each option has three subdirectories: PaDEL, RDKit and 
    combined_fingerprints.

#### RDKit_fingerprints(dataset, fp_name, index_col_nr, smiles_col_nr, inchikeys_col_nr; nBits=nothing, radius=nothing)    
    This function takes a DataFrame, a fingerprint name that corresponds to one in the dictionary 'dict', and the 
    column numbers of indices, SMILES and inchikeys in the DataFrame. nBits and radius are optional arguments that 
    only apply for certain fingerprints, specified in 'dict'. The function generates RDKit fingerprints for each 
    molecule in the data set using their SMILES. The function returns a DataFrame with indices, inchikeys and the 
    generated fingerprints, and the feature names. 
    
    Parameters:  
    - dataset: DataFrame; contains dataset. 
    - fp_name: Str; should match one of the RDKit fingerprint names.
    - index_col_nr: Int.
    - smiles_col_nr: Int.
    - inchikeys_col_nr: Int.
    - nBits=nothing: Int. 
    - radius=nothing: Int; 
    
    Returns:
    - total_df: DataFrame; contains optionally a column of indices of type Int64 and a column of inchikeys of type 
      Str31, followed by X data elements of type Float64. 
    - headers_RDK: Vector with elements of type Str.
    
#### PaDEL_fingerprints(dataset, fp_name, index_col_nr, smiles_col_nr, inchikeys_col_nr)
    This function takes a DataFrame, a fingerprint name that corresponds to one in the descriptors.XML file, 
    and the column numbers of indices, SMILES and inchikeys in the DataFrame. It generates PaDEL fingerprints for 
    each molecule in the data set using their SMILES. The function returns a DataFrame with indices, inchikeys and 
    the generated fingerprints, and the feature names. 
    
    Parameters:
    - dataset: DataFrame; contains dataset.
    - fp_name: Str; should match one of the PaDEL fingerprint names.
    - index_col_nr: Int.
    - smiles_col_nr: Int.
    - inchikeys_col_nr: Int.
    
    Returns:
    - total_df: DataFrame; contains optionally a column of indices of type Int64 and a column of inchikeys of type 
      Str31, followed by X data elements of type Float64. 
    - headers_PaDEL: Vector with elements of type Str.
    
    
#### change_xml(wanted_fingerprint)
    This function takes a fingerprint name that corresponds to one in the descriptors.xml file and modifies the XML 
    file by changing the selected fingerprint value to true and the remaining fingerprint values to false. The 
    modified XML file is saved back to the descriptors.xml file. 
    
    Parameters:
    - wanted_fingerprint: Str; should match one of the PaDEL fingerprint names.
    
#### create_train_test_split(total_df, y_data; start_col_X_data=1, train_size=0.9, random_state=42)
    This function takes a DataFrame that contains independent variables, and a vector of dependent variables. Optional 
    arguments are the column number of the DataFrame where the independent variables start, a float between 0 and 1
    that represents the training set size, and an integer representing the random state. The function creates a train
    test split using the chosen parameters. The train and test set of the total DataFrame, the dependent variables and 
    the independent variables are returned as matrices. 
    
    Parameters:
    - total_df: DataFrame, contains optionally a column of indices of type Int and a column of inchikeys of type Str, 
      followed by X data elements of type Float64. 
    - y_data: Vector with elements of type Float64. 
    - start_col_X_data=1: Int.
    - train_size=0.9: Float between 0 and 1. 
    - random_state=42: Int.
    
    Returns:
    - total_train: Matrix; contains amongst others X train data elements of type Float64. 
    - total_test: Matrix; contains amongst others X test data elements of type Float64. 
    - y_train: Vector with elements of type Float64. 
    - y_test: Vector with elements of type Float64. 
    - X_train: Matrix with elements of type Float64. 
    - X_test: Matrix with elements of type Float64. 

#### train_model(X_train, y_train, X_test, y_test; n_estimators=600, min_samples_leaf=4, max_features="sqrt", random_state=42)
    This function takes train arrays of X and y data, and test matrices of X and y data. Optional arguments are the number 
    of estimators, a minimum number of samples in each leaf, a maximum number of features to consider when looking for the 
    best split, and an integer representing the random state. The function creates and predicts a Random Forest Regressor 
    using the chosen parameters and returns the trained model, the train and test score, and the train and test predictions.
    
    Parameters:
    - X_train: Matrix with elements of type Float64. 
    - y_train: Vector with elements of type Float64. 
    - X_test: Matrix with elements of type Float64. 
    - y_test: Vector with elements of type Float64. 
    - n_estimators=600: Int.
    - min_samples_leaf=4: Int.
    - max_features="sqrt": Int or Str; 1.0, "sqrt" or "log2".
    - random_state=42: Int.
    
    Returns: 
    - model: Sklearn RandomForestRegressor. 
    - train_score: Float.
    - test_score: Float.
    - train_predictions: Vector with elements of type Float.
    - test_predictions: Vector with elements of type Float.
    
#### cross_validation(model, X_train, y_train, test_score; cv=3)
    This function takes a model, arrays of independent and dependent variables, and a test score. The optional argument is an 
    integer that represents the number of folds for the cross validation. The function performs a cross validation and returns 
    an array of the cross validation score and R, which is the median of the test score and the cross validation scores. 
    
    Parameters:
    - model: Any object.
    - X_train: Matrix with elements of type Float64. 
    - y_train: Vector with elements of type Float64. 
    - test_score: Float.
    - cv=3: Int.
    
    Returns:
    - cross_score: Vector with elements of type Float.
    - R: Float; median of test_score and the elements of cross_score.
    
#### important_features(model, headers; variance_explained=85) 
    This function takes a model and an array of variable names. The optional argument is an integer between 0 and 100 that 
    represents the variance percentage you want to have explained by the important features. The function selects the most 
    important features based on the chosen variance explained. The function returns the sum of the importances, an array of 
    the headers, an array of the importances and the number of the selected important features. 
    
    Parameters:
    - model: Any object.
    - headers: Vector with elements of type Str.
    - variance_explained=85: Int between 0 and 100.
    
    Returns:
    - sum_importance: Float.
    - important_headers: Vector with elements of type Str.
    - importances: Vector with elements of type Float. 
    - n_important_features: Int.
    
#### grid_best_params(param_grid, X_train, y_train)
    This function takes a dictionary with "min_samples_leaf", "n_estimators" and "max_features" as keys and corresponding 
    vectors as values, an array with X data and a vector of y data. The function performs a 5-fold (default) 3D grid search 
    and returns the best parameters for min_samples_leaf, n_estimators and max_features. 
    
    Parameters:
    - param_grid: Dict with "min_samples_leaf" (Vector of Int), "n_estimators" (Vector of Int) and "max_features" (Vector 
      of 1.0, "sqrt", and/or "log2") as keys (values). 
    - X_train: Matrix with elements of type Float64.
    - y_train: Vector with elements of type Float64.
    
    Returns:
    - min_samples_leaf_opt: Int.
    - n_estimators_opt: Int.
    - max_features_opt: Int or Str; 1.0, "sqrt" or "log2".
    
#### plot_training_testing(fp_name, software_name, train_predictions, y_train, test_predictions, y_test; x_label="Measured LC50 (log(mg/L))", y_label="Predicted LC50 (log(mg/L))")
    This function takes a fingerprint name, a software name, and arrays of train predictions, train data, test predictions and 
    test data. The optional arguments are labels for the x and y axes. The function creates a scatter plot of predicted values 
    vs. measured values for the training and testing set, together with the line y = x, and saves the figure. 
    
    Parameters:
    - fp_name: Str; should match one of the PaDEL or RDKit fingerprint names.
    - software_name: Str; "PaDEL", "RDKit" or "combined_fingerprints".
    - train_predictions: Vector with elements of type Int or Float. 
    - y_train: Vector with elements of type Int or Float. 
    - test_predictions: Vector with elements of type Int or Float. 
    - y_test: Vector with elements of type Int or Float. 
    - x_label="Measured LC50 (log(mg/L))": Str.
    - y_label="Predicted LC50 (log(mg/L))": Str.
   
#### plot_residuals(fp_name, software_name, residuals_train, y_train, residuals_test, y_test; x_label="Measured LC50 (log(mg/L))", y_label="Residuals")
    This function takes a fingerprint name, a software name, and arrays of train residuals, train data, test residuals and test 
    data. The optional arguments are labels for the x and y axes. The function creates a scatter plot of residuals vs. measured 
    values for the training and testing set, together with the horizontal lines y1 = 1 and y2 = -1, and saves the figure. 
    
    Parameters:
    - fp_name: Str; should match one of the PaDEL or RDKit fingerprint names.
    - software_name: Str; "PaDEL", "RDKit" or "combined_fingerprints".
    - residuals_train: Vector with elements of type Int or Float. 
    - y_train: Vector with elements of type Int or Float. 
    - residuals_test: Vector with elements of type Int or Float. 
    - y_test: Vector with elements of type Int or Float. 
    - x_label="Measured LC50 (log(mg/L))": Str.
    - y_label="Residuals": Str.
    
#### plot_distribution_residuals(fp_name, software_name, residuals; x_label="Residuals", y_label="Frequency")
    This function takes a fingerprint name, a software name, and an array of residuals. The optional arguments are labels for the 
    x and y axes. The function creates a histogram that represents the distribution of the residuals, together with the horizontal 
    line y = 0, and saves the figure. 
    
    Parameters:
    - fp_name: Str; should match one of the PaDEL or RDKit fingerprint names.
    - software_name: Str; "PaDEL", "RDKit" or "combined_fingerprints".
    - residuals: Vector with elements of type Int or Float. 
    - x_label="Residuals": Str.
    - y_label="Frequency": Str.
    
#### plot_important_features(fp_name, software_name, important_features, importances; x_label="", y_label="Percentage importance")
    This function takes a fingerprint name, a software name, and arrays of important features and the corresponding importances. 
    The optional arguments are labels for the x and y axes. The function creates a bar plot of the important features and their 
    importance percentages, and saves the figure. 
    
    Parameters:
    - fp_name: Str; should match one of the PaDEL or RDKit fingerprint names.
    - software_name: Str; "PaDEL", "RDKit" or "combined_fingerprints".
    - important_features: Vector with elements of type Str.
    - importances: Vector with elements of type Int or Float.
    - x_label="": Str.
    - y_label="Percentage importance": Str.
    
#### function create_plots_of_summary(summary_opt_filename, software_name)
    This function takes a filename of a summary of optimized fingerprint models and the corresponding software name. The function 
    creates the following plots: 1) a scatter plot of predicted values vs. measured values for the training and testing set, 
    together with the line y = x; 2) a scatter plot of residuals vs. measured values for the training and testing set, together 
    with the horizontal lines y1 = 1 and y2 = -1; 3) a histogram that represents the distribution of the residuals, together with 
    the horizontal line y = 0; 4) a bar plot of the important features and their importance percentages, and saves the figures. 
    
    Parameters:
    - summary_opt_filename: Str (ending with ".csv"); CSV filename should not exist already.
    - software_name: Str; "PaDEL", "RDKit" or "combined_fingerprints".
    
#### function remove_low_scores(summary_opt; limit_train_score=0.8, limit_test_score=0.5)
    This function takes a CSV.File object that contains a summary and two optional arguments which are floats between 0 and 1 that 
    represent the minimum train and test score for the selection of the fingerprints. The function removes the fingerprints that do 
    not meet the minimum train and test scores and returns the new summary as a DataFrame.
    
    Parameters:
    - summary_opt: CSV.File object.
    - limit_train_score=0.8: Float between 0 and 1.
    - limit_test_score=0.5: Float between 0 and 1.
    
    Returns:
    - new_summary: DataFrame.
    
#### summary_best_scores(summary_opt_filename, software_name; limit_train_score=0.8, limit_test_score=0.5)
    This function takes a filename of a summary of optimized fingerprint models and the corresponding software name. The optional 
    arguments are two floats between 0 and 1, which represent the minimum train and test score for the selection of the fingerprints. 
    The function saves a subselection of the total summary containing scores, and creates a similar summary of scores with only the 
    fingerprints that meet the chosen minimum train and test scores. 
    
    Parameters:
    - summary_opt_filename: Str (ending with ".csv"); name of an existing CSV file.
    - software_name: Str; "PaDEL", "RDKit" or "combined_fingerprints".
    - limit_train_score=0.8: Float between 0 and 1.
    - limit_test_score=0.5: Float between 0 and 1.
    
#### create_combined_fp(filename_scores_PaDEL, filename_scores_RDK, filename_combined_fp)
    This function takes the existing filenames of selected PaDEL and RDKit summaries of scores, and a new filename for the combined 
    fingerprint. The function combines the most important features of the PaDEL and RDKit fingerprints and stores them in a new CSV 
    file.
    
    Parameters:
    - filename_scores_PaDEL: Str (ending with ".csv"); name of an existing CSV file.
    - filename_scores_RDK: Str (ending with ".csv"); name of an existing CSV file.
    - filename_combined_fp: Str (ending with ".csv"); CSV filename should not exist already.
    
#### train_optimize_fp(fp_filename, software_name, y_data, new_model_name, param_grid; variance_explained=85)
    This function takes a fingerprint filename, a software name, an array of y data, a new model name, and a dictionary with 
    "min_samples_leaf", "n_estimators" and "max_features" as keys and corresponding vectors as values. The optional argument is an 
    integer between 0 and 100 that represents the variance percentage you want to have explained by the important features. The 
    function creates a train test split, optimizes and trains a random forest regressor model, performs 3-fold cross validation, 
    collects the important features that explain the chosen variance. The function saves a summary of the optimized fingerprint 
    model and a summary of the scores as CSV files and the optimized model as a JOBLIB file. 
    
    Parameters:
    - fp_filename: Str (ending with ".csv"); name of an existing CSV file. 
    - software_name: Str; "PaDEL", "RDKit" or "combined_fingerprints".
    - y_data: Vector with elements of type Float or Int.
    - new_model_name: Str; model name should not exist already.
    - param_grid: Dict with "min_samples_leaf" (Vector of Int), "n_estimators" (Vector of Int) and "max_features" (Vector of 1.0, 
      "sqrt", and/or "log2") as keys (values).
    - variance_explained=85: Int between 0 and 100.
    
#### str_to_floatVec(string)
    This function takes a string of substrings separated by semicolons and parses the substrings to separated parts of type Float64. 
    The function returns the resulting Float vector. 
    
    Parameters:
    - string: Str of substrings separated by semicolons. Substrings should be able to be transformed to floats. 
    
    Returns:
    - Vector with elements of type Float64.
    
#### str_to_strVec(string)
    This function takes a string of substrings separated by semicolons and splits the substrings to separated parts of type Str. The 
    function returns the resulting Str vector. 
    
    Parameters:
    - string: Str of substrings separated by semicolons.
    
    Returns:
    - Vector with elements of type Str.

#### compute_residuals(y, predictions)
    This function takes a vector of y data and a vector of predictions and computes the residuals. The function returns the resulting 
    vector of residuals. 
    
    Parameters:
    - y: Vector with elements of type Float or Int.
    - predictions: Vector with elements of type Float or Int.
    
    Returns:
    - Vector of type Float or Int.

#### write_csv(csv_filename, data) = CSV.write(csv_filename, data)
    This function takes a CSV filename and a collection of type Dict, DataFrame or Table. The function saves the given data as a CSV 
    file with the chosen csv_filename. 
    
    Parameters:
    - csv_filename: Str (ending with ".csv"); CSV filename should not exist already.
    - data: Dict, DataFrame or Table.

#### save_model(model, model_name) = jl.dump(model, model_name)

    This function takes a model and a chosen model name and saves the model as a JOBLIB file with the chosen model_name.
    
    Parameters:
    - model: Any object. 
    - model_name: Str (ending with ".joblib"); JOBLIB filename should not exist already.

#### read_csv(csv_filename) = CSV.File(csv_filename)

    Parameters:
    - csv_filename: Str (ending with ".csv"); name of an existing CSV file.

    Returns:
    - CSV.File object.

### Structure

The following list describes the folder structure that will be created when the code is run:
- **/fingerprints**: contains all fingerprints that are generated
  - **/fingerprints/PaDEL**: all PaDEL fingerprints
  - **/fingerprints/RDKit**: all RDKit fingerprints
  - **/fingerprints/combined_fingerprints**: all combined fingerprints
  
- **/summaries**: contains all summaries with scores and results that are created
  - **/summaries/PaDEL**: all PaDEL summaries
  - **/summaries/RDKit**: all RDKit summaries
  - **/summaries/combined_fingerprints**: all summaries of the combined fingerprints

- **/models**: contains all models that are created
  - **/models/PaDEL**: all PaDEL models
  - **/models/RDKit**: all RDKit models
  - **/models/combined_fingerprints**: all models of the combined fingerprints

- **/figures**: contains all figures that are created
  - **/figures/PaDEL**: all PaDEL figures
  - **/figures/RDKit**: all RDKit figures
  - **/figures/combined_fingerprints**: all figures of the combined fingerprints

## Authors
- Melanie Messih

## Supervisors
- dr. Saer Samanipour
- Viktoriia Turkina MSc

## Research group, research institute and university
Analytical Chemistry Group, Van 't Hoff Institute for Molecular Sciences (HIMS), University of Amsterdam (UvA)
