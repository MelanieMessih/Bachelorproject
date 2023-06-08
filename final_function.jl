include("final_imports.jl")
include("final_functions.jl")

# Get data set, X data, their names, and y data
fish_toxicity_data, X_data, headers, y_data = read_data("toxicity_data_fish_desc.csv")

function create_best_fingerprint(dataset, y_data, smiles_col_nr; index_col_nr=nothing, inchikeys_col_nr=nothing, limit_train_score=0.8, limit_test_score=0.5, variance_explained=85)
    """
    This function takes a DataFrame of the dataset that contains columns with X data, y data, SMILES and optionally 
    indices and inchikeys, a vector of y data, and the column number of the DataFrame that contains SMILES. Optional
    parameters are the column numbers that contain indices and inchikeys, two floats between 0 and 1 that represent 
    the minimum train and test score for the selection of the fingerprints, and an integer between 0 and 100 that 
    represents the variance percentage you want to have explained by the important features. 
    
    The function generates PaDEL and RDKit fingerprints, creates and optimizes random forest regressor models with 
    them, selects the fingerprints of which the train and test score exceed the given limit_train_score and 
    limit_test_score, and combines the selected fingerprints into one fingerprint. This final fingerprint is used to
    create the final random forest model, which is trained and optimized. 
    
    During this process, summaries of the scores and results, the generated fingerprints, and the generated figures 
    are saved as CSV files and the created models are saved as JOBLIB files. 
    """
    # Create directories needed for this function
    create_all_directories()

    # Define the hyperparameter search space for the grid search
    n_estimators = collect(100:20:600)
    min_samples_leaf = collect(2:2:8)
    max_features = [1.0, "sqrt", "log2"]

    param_grid = Dict("n_estimators" => n_estimators, "min_samples_leaf" => min_samples_leaf, "max_features" => max_features)

    # Define where the X data starts, depending on how many additional columns there will be before them 
    if typeof(index_col_nr) === nothing
        # No additional columns before the X data
        if typeof(inchikeys_col_nr) === nothing
            start_col_X_data = 1
        # Only an inchikeys column before the X data
        else
            start_col_X_data = 2
        end
    # Only an indices column before the X data
    elseif typeof(inchikeys_col_nr) === nothing
        start_col_X_data = 2
    # Both an indices and an inchikeys column before the X data
    else
        start_col_X_data = 3
    end

    # Generate PaDEL fingerprints
    summary_PaDEL_opt = DataFrame()
    for name_fp in ["Fingerprinter", "ExtendedFingerprinter", "GraphOnlyFingerprinter", "PubchemFingerprinter", 
                    "SubstructureFingerprintCount", "KlekotaRothFingerprintCount", "AtomPairs2DFingerprintCount"]

        # Generate fingerprints for each molecule in the data set and save as CSV file
        total_df, headers_PaDEL = PaDEL_fingerprints(dataset, name_fp, index_col_nr, smiles_col_nr, inchikeys_col_nr)
        CSV.write("fingerprints/PaDEL/$name_fp.csv", total_df)

        # Create train test split of total DataFrame, y data and X data
        total_train_PaDEL, total_test_PaDEL, y_train_PaDEL, y_test_PaDEL, X_train_PaDEL, X_test_PaDEL = create_train_test_split(total_df, y_data, start_col_X_data=start_col_X_data)

        # Perform grid search and select best parameters
        min_samples_leaf_opt, n_estimators_opt, max_features_opt = grid_best_params(param_grid, X_train_PaDEL, y_train_PaDEL)

        # Train model with best parameters and perform 3-fold cross validation
        opt_model, opt_train_score, opt_test_score, opt_train_predictions, opt_test_predictions = train_model(X_train_PaDEL, y_train_PaDEL, X_test_PaDEL, y_test_PaDEL, n_estimators=n_estimators_opt, min_samples_leaf=min_samples_leaf_opt, max_features=max_features_opt)
        opt_cross_score, opt_R = cross_validation(opt_model, X_train_PaDEL, y_train_PaDEL, opt_test_score)
        
        # Get important feature information
        sum_importance, important_headers, importances, n_important_features = important_features(opt_model, headers_PaDEL)

        # Compute test and train residuals
        residuals_test = compute_residuals(y_test_PaDEL, opt_test_predictions)
        residuals_train = compute_residuals(y_train_PaDEL, opt_train_predictions)
            
        # Make sure the max feature column contains one type (String) only to avoid errors when combining them
        if typeof(max_features_opt) != String
            max_features_opt = string(max_features_opt)
        end

        # Create summary of the optimized fingerprint type
        summary_PaDEL_opt_fp = DataFrame(Name="$(name_fp)_opt_PaDEL", Length=length(headers_PaDEL), ModelName="opt_model_PaDEL_$(name_fp)", MinSamplesLeaf=min_samples_leaf_opt, 
                                    nEstimators=n_estimators_opt, MaxFeatures=max_features_opt, TrainScore=opt_train_score, TestScore=opt_test_score, 
                                    CrossValScore1=opt_cross_score[1], CrossValScore2=opt_cross_score[2], CrossValScore3=opt_cross_score[3], R=opt_R, yTrain=join(y_train_PaDEL, "; "), 
                                    yTest=join(y_test_PaDEL, "; "), TrainPredictions=join(opt_train_predictions, "; "), TestPredictions=join(opt_test_predictions, "; "), 
                                    Features=join(headers_PaDEL, "; "), ImportantFeatures=join(important_headers, "; "), nImportantFeatures=n_important_features, 
                                    SumImportantFeatures=sum_importance, Importances=join(importances, "; "), ResidualsTrain=join(residuals_train, "; "), 
                                    ResidualsTest=join(residuals_test, "; "))

        # Add the summary to the complete summary of all PaDEL optimized fingerprint types
        append!(summary_PaDEL_opt, summary_PaDEL_opt_fp)

        # Save summary and model of the optimized fingerprint type as CSV and JOBLIB files, respectively
        write_csv("summaries/PaDEL/summary_PaDEL_opt.csv", summary_PaDEL_opt)
        save_model(opt_model,"models/PaDEL/opt_model_PaDEL_$(name_fp).joblib")
    end

    # Generate RDKit fingerprints
    summary_RDK_opt = DataFrame()
    for name_fp in ["MACCS", "HashedMorgan_ECFP6", "ErG", "Avalon", "EState", "HashedTopologicalTorsion", "RDK", 
                    "RDKFP", "Layered", "Pattern"]
        # Generate fingerprints for each molecule in the data set and save as CSV file
        total_df, headers_RDK = RDKit_fingerprints(dataset, name_fp, index_col_nr, smiles_col_nr, inchikeys_col_nr)
        CSV.write("fingerprints/RDKit/$name_fp.csv", total_df)

        # Create train test split of total DataFrame, y data and X data
        total_train_RDK, total_test_RDK, y_train_RDK, y_test_RDK, X_train_RDK, X_test_RDK = create_train_test_split(total_df, y_data, start_col_X_data=start_col_X_data)

        # Perform grid search and select best parameters
        min_samples_leaf_opt, n_estimators_opt, max_features_opt = grid_best_params(param_grid, X_train_RDK, y_train_RDK)

        # Train model with best parameters and perform 3-fold cross validation
        opt_model, opt_train_score, opt_test_score, opt_train_predictions, opt_test_predictions = train_model(X_train_RDK, y_train_RDK, X_test_RDK, y_test_RDK, n_estimators=n_estimators_opt, min_samples_leaf=min_samples_leaf_opt, max_features=max_features_opt)
        opt_cross_score, opt_R = cross_validation(opt_model, X_train_RDK, y_train_RDK, opt_test_score)

        # Get important feature information
        sum_importance, important_headers, importances, n_important_features = important_features(opt_model, headers_RDK)

        # Compute test and train residuals
        residuals_test = compute_residuals(y_test_RDK, opt_test_predictions)
        residuals_train = compute_residuals(y_train_RDK, opt_train_predictions)

        # Make sure the max feature column contains one type (String) only to avoid errors when combining them
        if typeof(max_features_opt) != String
            max_features_opt = string(max_features_opt)
        end

        # Create summary of the optimized fingerprint type
        summary_RDK_opt_fp = DataFrame(Name="$(name_fp)_opt_RDK", Length=length(headers_RDK), ModelName="opt_model_RDK_$(name_fp)", MinSamplesLeaf=min_samples_leaf_opt, 
                                    nEstimators=n_estimators_opt, MaxFeatures=max_features_opt, TrainScore=opt_train_score, TestScore=opt_test_score, 
                                    CrossValScore1=opt_cross_score[1], CrossValScore2=opt_cross_score[2], CrossValScore3=opt_cross_score[3], R=opt_R, yTrain=join(y_train_RDK, "; "), 
                                    yTest=join(y_test_RDK, "; "), TrainPredictions=join(opt_train_predictions, "; "), TestPredictions=join(opt_test_predictions, "; "), 
                                    Features=join(headers_RDK, "; "), ImportantFeatures=join(important_headers, "; "), nImportantFeatures=n_important_features, 
                                    SumImportantFeatures=sum_importance, Importances=join(importances, "; "), ResidualsTrain=join(residuals_train, "; "), 
                                    ResidualsTest=join(residuals_test, "; "))

        # Add the summary to the complete summary of all RDKit optimized fingerprint types
        append!(summary_RDK_opt, summary_RDK_opt_fp)

        # Save summary and model of the optimized fingerprint type as CSV and JOBLIB files, respectively
        write_csv("summaries/RDKit/summary_RDK_opt.csv", summary_RDK_opt)
        save_model(opt_model,"models/RDKit/opt_model_RDK_$(name_fp).joblib")
    end

    # Create plots of PaDEL and RDKit summaries
    create_plots_of_summary("summary_PaDEL_opt.csv", "PaDEL")
    create_plots_of_summary("summary_RDK_opt.csv", "RDKit")

    # Create PaDEL and RDKit summaries of best scores
    summary_best_scores("summary_PaDEL_opt.csv", "PaDEL", limit_train_score=limit_train_score, limit_test_score=limit_test_score)
    summary_best_scores("summary_RDK_opt.csv", "RDKit", limit_train_score=limit_train_score, limit_test_score=limit_test_score)

    # Create one combined fingerprint of the most important features of the selected PaDEL and RDKit fingerprints 
    create_combined_fp("summary_PaDEL_selected_scores.csv", "summary_RDK_selected_scores.csv", "combined_fp.csv")

    # Train and optimize model with the combined fingerprint
    train_optimize_fp("combined_fp.csv", "combined_fingerprints", y_data, "combined_fp", param_grid, variance_explained=variance_explained)

    # Create plots of the combined fingerprints summary
    create_plots_of_summary("summary_opt_combined_fp.csv", "combined_fingerprints")
end

create_best_fingerprint(fish_toxicity_data, y_data, 2, index_col_nr=1, inchikeys_col_nr=4)