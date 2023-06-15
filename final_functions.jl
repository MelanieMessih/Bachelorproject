# Only for fish toxicity data set 
function remove_features(X_data, headers=[])
    """
    This function takes a Dataframe or Matrix as input. When a DataFrame is given, "headers" can be left empty. 
    When a Matrix is given, headers should be spicified to make sure the column headers are known when working 
    with the data. The function removes columns with missing, nan or inf values and returns the cleaned up 
    Dataframe or Matrix and its headers. 
    """
    # Headers can be defined by selecting the names of the DataFrame
    if isempty(headers)        
        headers = names(X_data)    
    end    

    # Remove columns with missing, nan or inf values
    headers = headers[vec(.!any(ismissing.(Matrix(X_data)), dims=1))]    
    X_data = X_data[:, vec(.!any(ismissing.(Matrix(X_data)), dims=1))]    
    headers = headers[vec(.!any(isnan.(Matrix(X_data)), dims=1))]    
    X_data = X_data[:, vec(.!any(isnan.(Matrix(X_data)), dims=1))]
    headers = headers[vec(.!any(Inf.==(Matrix(X_data)), dims=1))]    
    X_data = X_data[:, vec(.!any(Inf.==(Matrix(X_data)), dims=1))]    
    headers = headers[vec(.!any((-Inf).==(Matrix(X_data)), dims=1))]    
    X_data = X_data[:, vec(.!any((-Inf).==(Matrix(X_data)), dims=1))]    

    return X_data, headers
end

# Only for fish toxicity data set
function read_data(csv_file_name)
    """
    This function takes a CSV file name as input and reads it as a DataFrame. The independent and dependent 
    variables are selected and feature columns containing missing, nan or inf values are removed. The 
    function returns the DataFrame, independent variables, feature names, and dependent variables.
    """
    # Read data set
    data_name = CSV.read(csv_file_name, DataFrame)

    # Select independent variables, without text values
    X_data = data_name[:, 8:end]
    X_data, headers = remove_features(X_data)

    # Select dependent variables
    y_data = data_name[:, 6]

    return data_name, X_data, headers, y_data
end

function create_directory(path)
    """
    This function takes a path. If the path does not exist, the specified path of directories 
    is created.
    """
    if !isdir(path)
        mkpath(path)
        println("Directory created: $path")
    else
        println("Directory already exists: $path")
    end
end

function create_all_directories()
    """
    This function creates all directories needed for this project. Specifically, the function
    creates directories for summaries, figures, fingerprints and models, and each option has
    three subdirectories: PaDEL, RDKit and combined_fingerprints.
    """
    create_directory("summaries/PaDEL")
    create_directory("summaries/RDKit")
    create_directory("summaries/combined_fingerprints")

    create_directory("figures/PaDEL")
    create_directory("figures/RDKit")
    create_directory("figures/combined_fingerprints")

    create_directory("fingerprints/PaDEL")
    create_directory("fingerprints/RDKit")
    create_directory("fingerprints/combined_fingerprints")

    create_directory("models/PaDEL")
    create_directory("models/RDKit")
    create_directory("models/combined_fingerprints")
end

function RDKit_fingerprints(dataset, fp_name, index_col_nr, smiles_col_nr, inchikeys_col_nr; nBits=nothing, radius=nothing)    
    """
    This function takes a DataFrame, a fingerprint name that corresponds to one in the dictionary 'dict' below, 
    and the column numbers of indices, SMILES and inchikeys in the DataFrame. nBits and radius are optional arguments 
    that only apply for certain fingerprints, specified in 'dict'. The function generates RDKit fingerprints for each 
    molecule in the data set using their SMILES. The function returns a DataFrame with indices, inchikeys and the 
    generated fingerprints, and the feature names. 
    """
    # Dictionary containing information about various RDKit fingerprint types
    dict = Dict("MACCS" => Dict("call" => AllChem.GetMACCSKeysFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 167),
                "Morgan_ECFP2" => Dict("call" => AllChem.GetMorganFingerprint, "params" => Dict("nBits" => nothing, "radius"=> 2), "length_fp" => 2048),
                "Morgan_ECFP4" => Dict("call" => AllChem.GetMorganFingerprint, "params" => Dict("nBits" => nothing, "radius"=> 4), "length_fp" => 2048),
                "Morgan_ECFP6" => Dict("call" => AllChem.GetMorganFingerprint, "params" => Dict("nBits" => nothing, "radius"=> 6), "length_fp" => 2048),
                "Morgan_ECFP8" => Dict("call" => AllChem.GetMorganFingerprint, "params" => Dict("nBits" => nothing, "radius"=> 8), "length_fp" => 2048),
                "MorganAsBitVect_ECFP2" => Dict("call" => AllChem.GetMorganFingerprintAsBitVect, "params" => Dict("nBits" => 2048, "radius" => 2)),
                "MorganAsBitVect_ECFP4" => Dict("call" => AllChem.GetMorganFingerprintAsBitVect, "params" => Dict("nBits" => 2048, "radius" => 4)),
                "MorganAsBitVect_ECFP6" => Dict("call" => AllChem.GetMorganFingerprintAsBitVect, "params" => Dict("nBits" => 2048, "radius" => 6)),
                "MorganAsBitVect_ECFP8" => Dict("call" => AllChem.GetMorganFingerprintAsBitVect, "params" => Dict("nBits" => 2048, "radius" => 8)),
                "HashedMorgan_ECFP2" => Dict("call" => AllChem.GetHashedMorganFingerprint, "params" => Dict("nBits" => nothing, "radius" => 2), "length_fp" => 2048),
                "HashedMorgan_ECFP4" => Dict("call" => AllChem.GetHashedMorganFingerprint, "params" => Dict("nBits" => nothing, "radius" => 4), "length_fp" => 2048),
                "HashedMorgan_ECFP6" => Dict("call" => AllChem.GetHashedMorganFingerprint, "params" => Dict("nBits" => nothing, "radius" => 6), "length_fp" => 2048),
                "HashedMorgan_ECFP8" => Dict("call" => AllChem.GetHashedMorganFingerprint, "params" => Dict("nBits" => nothing, "radius" => 8), "length_fp" => 2048),
                "AtomPair" => Dict("call" => AllChem.GetAtomPairFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 1024),
                "HashedAtomPair" => Dict("call" => AllChem.GetHashedAtomPairFingerprint, "params" => Dict("nBits" => 1024, "radius" => nothing)),
                "HashedAtomPairAsBitVect" => Dict("call" => AllChem.GetHashedAtomPairFingerprintAsBitVect, "params" => Dict("nBits" => 1024, "radius" => nothing)),
                "ErG" => Dict("call" => rdReducedGraphs.GetErGFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 315),
                "Avalon" => Dict("call" => FingerprintUtils.BuildAvalonFP, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 512),
                "EState" => Dict("call"=> Fingerprinter.FingerprintMol, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 79),
                "TopologicalTorsion" => Dict("call" => AllChem.GetTopologicalTorsionFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 2048),
                "HashedTopologicalTorsion" => Dict("call" => AllChem.GetHashedTopologicalTorsionFingerprint, "params" => Dict("nBits" => 2048, "radius" => nothing)),
                "HashedTopologicalTorsionAsBitVect" => Dict("call" => AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect, "params" => Dict("nBits" => 2048, "radius" => nothing)),
                "RDK" => Dict("call" => AllChem.RDKFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 2048),
                "RDKFP" => Dict("call" => FingerprintUtils.BuildRDKitFP, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 2048),
                "UnfoldedRDKCountBased" => Dict("call" => rdmolops.UnfoldedRDKFingerprintCountBased, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 2048),
                "Layered" => Dict("call" => rdmolops.LayeredFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 2048),
                "Pattern" => Dict("call" => rdmolops.PatternFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 2048),
                "BP" => Dict("call" => Sheridan.GetBPFingerprint, "params" => Dict("nBits" => nothing, "radius" => nothing), "length_fp" => 2048))

    # Determine the additional parameters and length of the fingerprint based on the specified fingerprint type
    # mol
    if isnothing(dict[fp_name]["params"]["nBits"]) && isnothing(dict[fp_name]["params"]["radius"])
        add_params = nothing
        length_fp = dict[fp_name]["length_fp"]
 
    # mol, nBits (given)
    elseif isnothing(dict[fp_name]["params"]["radius"]) && !isnothing(nBits)
        add_params = nBits
        length_fp = nBits
    
    # mol, nBits (default)
    elseif isnothing(dict[fp_name]["params"]["radius"]) && isnothing(nBits)
        add_params = dict[fp_name]["params"]["nBits"]
        length_fp = dict[fp_name]["params"]["nBits"]
    
    # mol, radius (given)
    elseif isnothing(dict[fp_name]["params"]["nBits"]) && !isnothing(radius)
        add_params = radius
        length_fp = dict[fp_name]["length_fp"]
    
    # mol, radius (default)
    elseif isnothing(dict[fp_name]["params"]["nBits"]) && isnothing(radius)
        add_params = dict[fp_name]["params"]["radius"]
        length_fp = dict[fp_name]["length_fp"]
    
    else
        # mol, nBits (default), radius (default)
        if isnothing(nBits) && isnothing(radius)
            add_params = (dict[fp_name]["params"]["radius"], dict[fp_name]["params"]["nBits"])
            length_fp = dict[fp_name]["params"]["nBits"]
            # println(add_params)
    
        # mol, nBits (given), radius (given)
        elseif isnothing(nBits) && !isnothing(radius)
            add_params = (radius, nBits)
            length_fp = nBits
    
        # mol, nBits (given), radius (default)
        elseif isnothing(radius)
            add_params = (dict[fp_name]["params"]["radius"], nBits)
            length_fp = nBits
        
        # mol, nBits (default), radius (given)
        else
            add_params = (radius, dict[fp_name]["params"]["nBits"])
            length_fp = dict[fp_name]["params"]["nBits"]
        end
    end 

    # Create feature names and create a DataFrame to store the generated fingerprints
    nr_rows = size(dataset, 1)
    X_data_RDK = DataFrame(["$(fp_name) $i"=>Vector{Float64}(undef, nr_rows) for i in 1:length_fp])
    headers_RDK = names(X_data_RDK)
    
    # Generate fingerprints for each molecule in the dataset
    for i in 1:nr_rows
        # Select SMILES from the correct column in the data set and convert it to a Mol file 
        smiles = dataset[i,smiles_col_nr]
        mol = rdk.MolFromSmiles(smiles)

        # Generate the correct fingerprint type 
        if fp_name == "EState"
            one_fp = dict[fp_name]["call"](mol)[1] 
        elseif isnothing(add_params)
            one_fp = dict[fp_name]["call"](mol)
        elseif length(add_params) == 2
            one_fp = dict[fp_name]["call"](mol, add_params[1], add_params[2])
        else
            one_fp = dict[fp_name]["call"](mol, add_params)
        end
        
        # Convert generated fingerprints to the corresponding row of the DataFrame
        one_fp = [one_fp[i] for i=1:length_fp]
        X_data_RDK[i, :] = one_fp
    end

    # Define whether there is an indices and/or inchikeys column, and combine them with the generated fingerprints
    if typeof(index_col_nr) === nothing
        # No indices or inchikeys columns
        if typeof(inchikeys_col_nr) === nothing
            total_df = X_data_RDK
        # Only an inchikeys column
        else
            inchikeys = select(dataset, inchikeys_col_nr)
            total_df = hcat(inchikeys, X_data_RDK)
        end
    # Only an indices column
    elseif typeof(inchikeys_col_nr) === nothing
        indices = select(dataset, index_col_nr)
        total_df = hcat(indices, X_data_RDK)
    # Both an indices and an inchikeys column
    else
        indices_inchikeys = select(dataset, index_col_nr, inchikeys_col_nr)
        total_df = hcat(indices_inchikeys, X_data_RDK)
    end

    return total_df, headers_RDK
end

function PaDEL_fingerprints(dataset, fp_name, index_col_nr, smiles_col_nr, inchikeys_col_nr)
    """
    This function takes a DataFrame, a fingerprint name that corresponds to one in the descriptors.XML file, 
    and the column numbers of indices, SMILES and inchikeys in the DataFrame. It generates PaDEL fingerprints for 
    each molecule in the data set using their SMILES. The function returns a DataFrame with indices, inchikeys and 
    the generated fingerprints, and the feature names. 
    """
    # Select feature names and create a DataFrame to store the generated fingerprints
    nr_rows = size(dataset, 1)
    X_data_PaDEL = DataFrame()

    # Select the correct fingerprint type in the XML file
    change_xml(fp_name)

    # Generate fingerprints for each molecule in the dataset
    for i in 1:nr_rows
        # Select SMILES from the correct column in the data set
        smiles = dataset[i,smiles_col_nr]
        # Generate the specific fingerprint type and store it 
        one_fp = DataFrame(pd.from_smiles(smiles, fingerprints=true, descriptors=false))
        append!(X_data_PaDEL, one_fp)
    end

    headers_PaDEL = names(X_data_PaDEL)

    # Define whether there is an indices and/or inchikeys column, and combine them with the generated fingerprints
    if typeof(index_col_nr) === nothing
        # No indices or inchikeys columns
        if typeof(inchikeys_col_nr) === nothing
            total_df = X_data_PaDEL
        # Only an inchikeys column
        else
            inchikeys = select(dataset, inchikeys_col_nr)
            total_df = hcat(inchikeys, X_data_PaDEL)
        end
    # Only an indices column
    elseif typeof(inchikeys_col_nr) === nothing
        indices = select(dataset, index_col_nr)
        total_df = hcat(indices, X_data_PaDEL)
    # Both an indices and an inchikeys column
    else
        indices_inchikeys = select(dataset, index_col_nr, inchikeys_col_nr)
        total_df = hcat(indices_inchikeys, X_data_PaDEL)
    end

    return total_df, headers_PaDEL
end

function change_xml(wanted_fingerprint)
    """
    This function takes a fingerprint name that corresponds to one in the descriptors.xml file and modifies 
    the XML file by changing the selected fingerprint value to true and the remaining fingerprint values to 
    false. The modified XML file is saved back to the descriptors.xml file. 
    """
    data = parse_file("descriptors.xml")
    xroot = root(data)

    # Navigate to fingerprints in XML file
    for c in child_elements(xroot)
        if has_children(c)
            if attribute(c, "name") == "Fingerprint"
                # Change selected fingerprint value to true and the remaining fingerprint values to false
                for fingerprint in child_elements(c)
                    if attribute(fingerprint, "name") == wanted_fingerprint
                        set_attributes(fingerprint, name=attribute(fingerprint, "name"), value="true")
                    else
                        set_attributes(fingerprint, name=attribute(fingerprint, "name"), value="false")
                    end
                end
            end
        end
        # println(c)
    end

    # Save modified XML file 
    save_file(data, "descriptors.xml")
end

function create_train_test_split(total_df, y_data; start_col_X_data=1, train_size=0.9, random_state=42)
    """
    This function takes a DataFrame that contains independent variables, and a vector of dependent variables. Optional 
    arguments are the column number of the DataFrame where the independent variables start, a float between 0 and 1
    that represents the training set size, and an integer representing the random state. The function creates a train
    test split using the chosen parameters. The train and test set of the total DataFrame, the dependent variables and 
    the independent variables are returned as matrices. 
    """
    # Create train test split of total DataFrame and dependent variables using the chosen parameters
    total_train, total_test, y_train, y_test = train_test_split(Matrix(total_df), y_data, train_size=train_size, random_state=random_state)

    # Select train and test set of independent variables 
    X_train = total_train[:, start_col_X_data:end]
    X_test = total_test[:, start_col_X_data:end]
    
    return total_train, total_test, y_train, y_test, X_train, X_test
end

function train_model(X_train, y_train, X_test, y_test; n_estimators=600, min_samples_leaf=4, max_features="sqrt", random_state=42)
    """
    This function takes train matrices of independent and dependent variables, and test matrices of independent and 
    dependent variables. Optional arguments are the number of estimators, a minimum number of samples in each leaf, a 
    maximum number of features to consider when looking for the best split, and an integer representing the random state.
    The function creates and predicts a Random Forest Regressor using the chosen parameters and returns the trained model,
    the train and test score, and the train and test predictions.
    """
    # Create and predict a Random Forest Regressor
    rfr = skl.RandomForestRegressor
    model = rfr(n_estimators=n_estimators, min_samples_leaf=min_samples_leaf, random_state=random_state, max_features=max_features).fit(X_train, y_train)

    train_score = model.score(X_train, y_train)
    test_score = model.score(X_test, y_test)

    train_predictions = model.predict(X_train)
    test_predictions = model.predict(X_test)

    return model, train_score, test_score, train_predictions, test_predictions
end

function cross_validation(model, X_train, y_train, test_score; cv=3)
    """
    This function takes a model, matrices of independent and dependent variables, and a test score. The 
    optional argument is an integer that represents the number of folds for the cross validation. The 
    function performs a cross validation and returns an array of the cross validation score and R, which 
    is the median of the test score and the cross validation scores. 
    """
    # Cross validation of model
    cross_score = cross_val_score(model, X_train, y_train, cv=cv)

    # Compute median of cross-validation scores and test score
    to_median = push!(cross_score, test_score)
    R = median(to_median)

    return cross_score, R
end

function important_features(model, headers; variance_explained=85) 
    """
    This function takes a model and an array of variable names. The optional argument is an integer between 
    0 and 100 that represents the variance percentage you want to have explained by the important features. 
    The function selects the most important features based on the chosen variance explained. The function 
    returns the sum of the importances, an array of the headers, an array of the importances and the number of 
    the selected important features. 
    """
    # Get feature importance percentages of all features, and sort them from large to small
    importance = model.feature_importances_ .* 100
    importance_rev = sortperm(importance, rev=true)
    
    i = 0
    sum_importances = 0

    # Select the most important features based on the chosen variance explained
    while sum_importances < variance_explained
        i += 1
        sum_importances += importance[importance_rev[i]]
    end
    println(i, sum_importances, importance_rev[1:i])

    important_headers = headers[importance_rev[1:i]]
    importances = importance[importance_rev[1:i]]
    n_important_features = length(importances)

    return sum_importances, important_headers, importances, n_important_features
end

function grid_best_params(param_grid, X_train, y_train)
    """
    This function takes a dictionary with "min_samples_leaf", "n_estimators" and "max_features" as keys and 
    corresponding vectors as values, an array with X data and a vector of y data. The function performs a 
    5-fold (default) 3D grid search and returns the best parameters for min_samples_leaf, n_estimators and 
    max_features. 
    """
    # Perform grid search
    rf = skl.RandomForestRegressor() # YOU MAY HAVE TO MOVE THIS TO OUTSIDE THE FUNCTION
    rf_grid = GridSearchCV(estimator=rf,param_grid=param_grid, verbose=4)
    rf_grid.fit(X_train, y_train)

    # Select best parameters
    min_samples_leaf_opt = rf_grid.best_params_["min_samples_leaf"]
    n_estimators_opt = rf_grid.best_params_["n_estimators"]
    max_features_opt = rf_grid.best_params_["max_features"]

    return min_samples_leaf_opt, n_estimators_opt, max_features_opt
end

function plot_training_testing(fp_name, software_name, train_predictions, y_train, test_predictions, y_test; x_label="Measured LC50 (log(mg/L))", y_label="Predicted LC50 (log(mg/L))")
    """
    This function takes a fingerprint name, a software name, and arrays of train predictions, train data, test 
    predictions and test data. The optional arguments are labels for the x and y axes. The function creates a
    scatter plot of predicted values vs. measured values for the training and testing set, together with the line 
    y = x, and saves the figure. 
    """
    Plots.gr(dpi=300)

    # Add line y = x
    plot(y_train, y_train, linecolor="blue", label=nothing)

    # Plot training and testing predictions vs. measurements
    scatter!(y_train, train_predictions, markercolor="cornflowerblue", label="training set", grid=true)
    scatter!(y_test, test_predictions, markercolor="orange", label="testing set")
    xlabel!(x_label)
    ylabel!(y_label)

    savefig("figures/$(software_name)/$(fp_name)_train_test.png")
end

function plot_residuals(fp_name, software_name, residuals_train, y_train, residuals_test, y_test; x_label="Measured LC50 (log(mg/L))", y_label="Residuals")
    """
    This function takes a fingerprint name, a software name, and arrays of train residuals, train data, test 
    residuals and test data. The optional arguments are labels for the x and y axes. The function creates a
    scatter plot of residuals vs. measured values for the training and testing set, together with the horizontal 
    lines y1 = 1 and y2 = -1, and saves the figure. 
    """
    Plots.gr(dpi=300)

    # Add horizontal lines at y=1 and y = -1 that span the width of the plot
    hline([1], linecolor="blue", label=nothing)
    hline!([-1], linecolor="blue", label=nothing)

    # Plot training and testing residuals vs. measurements
    scatter!(y_train, residuals_train, markercolor="cornflowerblue", label="training set", grid=true)
    scatter!(y_test, residuals_test, markercolor="orange", label="testing set")
    xlabel!(x_label)
    ylabel!(y_label)

    savefig("figures/$(software_name)/$(fp_name)_residuals.png")
end

function plot_distribution_residuals(fp_name, software_name, residuals; x_label="Residuals", y_label="Frequency")
    """
    This function takes a fingerprint name, a software name, and an array of residuals. The optional arguments are 
    labels for the x and y axes. The function creates a histogram that represents the distribution of the residuals, 
    together with the horizontal line y = 0, and saves the figure. 
    """
    Plots.gr(dpi=300)

    # Add horizontal line at y = 0 that spans the width of the plot
    hline([0], linecolor="blue", label=nothing)

    # Plot distribution of residuals for given residuals 
    histogram!(residuals, bins=50, label=nothing, color="cornflowerblue", xlabel=x_label, ylabel=y_label)

    savefig("figures/$(software_name)/$(fp_name)_distribition.png")
end

function plot_important_features(fp_name, software_name, important_features, importances; x_label="", y_label="Percentage importance")
    """
    This function takes a fingerprint name, a software name, and arrays of important features and the 
    corresponding importances. The optional arguments are labels for the x and y axes. The function creates
    a bar plot of the important features and their importance percentages, and saves the figure. 
    """
    Plots.gr(dpi=300)

    # Plot important features and their importance percentages
    bar_plot = bar(important_features, importances, xlabel=x_label, ylabel=y_label, legend=false, xguidefontsize=8, xrotation=80)
    plot!(bar_plot, bottom_margin=20mm)

    savefig("figures/$(software_name)/$(fp_name)_importances.png")
end

function create_plots_of_summary(summary_opt_filename, software_name)
    """
    This function takes a filename of a summary of optimized fingerprint models and the corresponding software 
    name. The function creates the following plots: 1) a scatter plot of predicted values vs. measured values 
    for the training and testing set, together with the line y = x; 2) a scatter plot of residuals vs. measured 
    values for the training and testing set, together with the horizontal lines y1 = 1 and y2 = -1; 3) a histogram 
    that represents the distribution of the residuals, together with the horizontal line y = 0; 4) a bar plot of 
    the important features and their importance percentages, and saves the figures. 
    """
    summary_opt = read_csv("summaries/$(software_name)/$(summary_opt_filename)")

    # Create plots for each fingerprint documented in summary_opt
    for row in summary_opt
        name_opt = row.Name
    
        y_train_opt = str_to_floatVec(row.yTrain)
        y_test_opt = str_to_floatVec(row.yTest)
        train_predictions_opt = str_to_floatVec(row.TrainPredictions)
        test_predictions_opt = str_to_floatVec(row.TestPredictions)
        # Plot the measured vs the predicted 96h LC50 values (training and testing set)
        plot_training_testing(name_opt, software_name, train_predictions_opt, y_train_opt, test_predictions_opt, y_test_opt)
    
        residuals_train_opt = str_to_floatVec(row.ResidualsTrain)
        residuals_test_opt = str_to_floatVec(row.ResidualsTest)
        # Plot residuals vs the measured 96h LC50 values (training and testing set)
        plot_residuals(name_opt, software_name, residuals_train_opt, y_train_opt, residuals_test_opt, y_test_opt)
    
        # Combine all residuals (training and testing set) and plot the distribution of the residuals
        residuals_total_opt = vcat(residuals_train_opt, residuals_test_opt)
        plot_distribution_residuals(name_opt, software_name, residuals_total_opt)
    
        important_features_opt = str_to_strVec(row.ImportantFeatures)
        importances_opt = str_to_floatVec(row.Importances)
        # Plot most important features with their importances
        plot_important_features(name_opt, software_name, important_features_opt, importances_opt)
    end
end

function remove_low_scores(summary_opt; limit_train_score=0.8, limit_test_score=0.5)
    """
    This function takes a CSV.File object that contains a summary and two optional arguments which are floats 
    between 0 and 1 that represent the minimum train and test score for the selection of the fingerprints. The 
    function removes the fingerprints that do not meet the minimum train and test scores and returns the new 
    summary as a DataFrame.
    """
    summary_opt_df = DataFrame(summary_opt)
    list_indices = []

    # Collect indices of fingerprints that do not meet the minimum train and test scores
    for i in 1:size(summary_opt_df, 1)
        train_score = summary_opt_df[i, 7]
        test_score = summary_opt_df[i, 8]
        if train_score .< limit_train_score
            append!(list_indices, i)
        elseif test_score .< limit_test_score
            append!(list_indices, i)
        end
    end
    
    # Delete fingerprints that do not meet the minimum train and test scores
    new_summary = delete!(summary_opt_df, list_indices)
    return new_summary
end 

function summary_best_scores(summary_opt_filename, software_name; limit_train_score=0.8, limit_test_score=0.5)
    """
    This function takes a filename of a summary of optimized fingerprint models and the corresponding software 
    name. The optional arguments are two floats between 0 and 1, which represent the minimum train and test score 
    for the selection of the fingerprints. The function saves a subselection of the total summary containing
    scores, and creates a similar summary of scores with only the fingerprints that meet the chosen minimum train 
    and test scores. 
    """
    # Define names needed to save files to the correct location with the correct filenames
    if software_name == "PaDEL"
        short_name = "PaDEL"
    elseif software_name == "RDKit"
        short_name = "RDK"
    else
        println("The software name should either be PaDEL or RDKit (String)")
    end

    # Read summary of optimized fingerprint models and create a summary of scores (subselection)
    summary_opt_df = CSV.read("summaries/$(software_name)/$(summary_opt_filename)", DataFrame)
    summary_opt_scores = select(summary_opt_df, [:Name, :Length, :ModelName, :MinSamplesLeaf, :nEstimators, :MaxFeatures, 
                                                        :TrainScore, :TestScore, :CrossValScore1, :CrossValScore2, :CrossValScore3, :R,
                                                        :ImportantFeatures, :nImportantFeatures, :SumImportantFeatures, :Importances])
    write_csv("summaries/$(software_name)/summary_$(short_name)_opt_scores.csv", summary_opt_scores)

    # Copy the summary of scores and remove fingerprints that do not meet the minimum train and test scores
    cp("summaries/$(software_name)/summary_$(short_name)_opt_scores.csv", "summaries/$(software_name)/summary_$(short_name)_selected_scores.csv", force=true)
    summary_selected_scores = read_csv("summaries/$(software_name)/summary_$(short_name)_selected_scores.csv")
    summary_selected_scores = remove_low_scores(summary_selected_scores, limit_train_score=limit_train_score, limit_test_score=limit_test_score)
    write_csv("summaries/$(software_name)/summary_$(short_name)_selected_scores.csv", summary_selected_scores)
end

function create_combined_fp(filename_scores_PaDEL, filename_scores_RDK, filename_combined_fp)
    """
    This function takes the existing filenames of selected PaDEL and RDKit summaries of scores, and a new filename 
    for the combined fingerprint. The function combines the most important features of the PaDEL and RDKit 
    fingerprints and stores them in a new CSV file.
    """
    # Open PaDEL and RDKit summaries of selected scores
    summary_PaDEL_selected_scores = read_csv("summaries/PaDEL/$(filename_scores_PaDEL)")
    summary_RDK_selected_scores = read_csv("summaries/RDKit/$(filename_scores_RDK)")

    # DataFrame to store all selected PaDEL and RDKit fingerprints
    combined_fingerprints_df = DataFrame()

    for row in 1:length(summary_PaDEL_selected_scores)
        # Open the PaDEL fingerprint file for the selected fingerprint
        name_fp = summary_PaDEL_selected_scores[row].Name[1:end-10]
        fp = read_csv("fingerprints/PaDEL/$(name_fp).csv")

        # Split string of important features 
        str_important_features = summary_PaDEL_selected_scores[row].ImportantFeatures
        split_important_features = split(str_important_features, "; ")

        # Include each important feature in the combined fingerprints DataFrame
        for important_feature in split_important_features
            important_feature = String(important_feature)
            combined_fingerprints_df[!, important_feature] = fp[important_feature]
        end
    end

    for row in 1:length(summary_RDK_selected_scores)
        # Open the RDKit fingerprint file for the selected fingerprint
        name_fp = summary_RDK_selected_scores[row].Name[1:end-8]
        fp = read_csv("fingerprints/RDKit/$(name_fp).csv")

        # Split string of important features 
        str_important_features = summary_RDK_selected_scores[row].ImportantFeatures
        split_important_features = split(str_important_features, "; ")

        # Include each important feature in the combined fingerprints DataFrame
        for important_feature in split_important_features
            important_feature = String(important_feature)
            combined_fingerprints_df[!, important_feature] = fp[important_feature]
        end
    end

    # Save combined fingerprint consisting of the most important features of the PaDEL and RDKit fingerprints 
    CSV.write("fingerprints/combined_fingerprints/$(filename_combined_fp)", combined_fingerprints_df)
end

function train_optimize_fp(fp_filename, software_name, y_data, new_model_name, param_grid; variance_explained=85)
    """
    This function takes a fingerprint filename, a software name, an array of y data, a new model name, and a dictionary 
    with "min_samples_leaf", "n_estimators" and "max_features" as keys and corresponding vectors as values. The optional 
    argument is an integer between 0 and 100 that represents the variance percentage you want to have explained by the 
    important features. The function creates a train test split, optimizes and trains a random forest regressor model, 
    performs 3-fold cross validation, collects the important features that explain the chosen variance. The function 
    saves a summary of the optimized fingerprint model and a summary of the scores as CSV files and the optimized model 
    as a JOBLIB file. 
    """
    combined_fingerprints_df = CSV.read("fingerprints/$(software_name)/$(fp_filename)", DataFrame)

    # Create train test split of total DataFrame, dependent variables and independent variables
    total_train, total_test, y_train, y_test, X_train, X_test = create_train_test_split(combined_fingerprints_df, y_data)

    # Perform grid search and select best parameters
    min_samples_leaf_opt, n_estimators_opt, max_features_opt = grid_best_params(param_grid, X_train, y_train)

    # Train model with best parameters and perform 3-fold cross validation
    opt_model, opt_train_score, opt_test_score, opt_train_predictions, opt_test_predictions = train_model(X_train, y_train, X_test, y_test, n_estimators=n_estimators_opt, min_samples_leaf=min_samples_leaf_opt, max_features=max_features_opt)
    opt_cross_score, opt_R = cross_validation(opt_model, X_train, y_train, opt_test_score)

    # Define feature names of the combined fingerprints 
    headers = names(combined_fingerprints_df)

    # Get important feature information
    sum_importance, important_headers, importances, n_important_features = important_features(opt_model, headers, variance_explained=variance_explained)
    
    # Compute test and train residuals
    residuals_test = compute_residuals(y_test, opt_test_predictions)
    residuals_train = compute_residuals(y_train, opt_train_predictions)

    # Create summary of the optimized combined fingerprint
    summary_opt = DataFrame(Name="$(new_model_name)_opt", Length=length(headers), ModelName="opt_model_$(new_model_name)", MinSamplesLeaf=min_samples_leaf_opt, 
                                nEstimators=n_estimators_opt, MaxFeatures=max_features_opt, TrainScore=opt_train_score, TestScore=opt_test_score, 
                                CrossValScore1=opt_cross_score[1], CrossValScore2=opt_cross_score[2], CrossValScore3=opt_cross_score[3], R=opt_R, yTrain=join(y_train, "; "), 
                                yTest=join(y_test, "; "), TrainPredictions=join(opt_train_predictions, "; "), TestPredictions=join(opt_test_predictions, "; "), 
                                Features=join(headers, "; "), ImportantFeatures=join(important_headers, "; "), nImportantFeatures=n_important_features, 
                                SumImportantFeatures=sum_importance, Importances=join(importances, "; "), ResidualsTrain=join(residuals_train, "; "), 
                                ResidualsTest=join(residuals_test, "; "))

    # Create a summary of scores of the optimized combined fingerprint (subselection)
    summary_opt_scores = select(summary_opt, [:Name, :Length, :ModelName, :MinSamplesLeaf, :nEstimators, :MaxFeatures, 
                                                            :TrainScore, :TestScore, :CrossValScore1, :CrossValScore2, :CrossValScore3, :R,
                                                            :ImportantFeatures, :nImportantFeatures, :SumImportantFeatures, :Importances])

    # Save summaries and model of the optimized combined fingerprint as CSV and JOBLIB files, respectively
    save_model(opt_model,"models/$(software_name)/opt_model_$(new_model_name).joblib")
    write_csv("summaries/$(software_name)/summary_opt_$(new_model_name).csv", summary_opt)
    write_csv("summaries/$(software_name)/summary_opt_scores_$(new_model_name).csv", summary_opt_scores)
end

str_to_floatVec(string) = parse.(Float64, split(string, "; "))
str_to_strVec(string) = split(string, "; ")
compute_residuals(y, predictions) = y - predictions

write_csv(csv_filename, data) = CSV.write(csv_filename, data)
save_model(model, model_name) = jl.dump(model, model_name)
read_csv(csv_filename) = CSV.File(csv_filename)
