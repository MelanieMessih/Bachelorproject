### Bachelor project 
- Bachelor programme: Chemistry (joint degree UvA/VU)
- Duration: 3 months (April 2023 - June 2023)

# Molecular fingerprints

create_directory(path::AbstractString)
  This function takes a path. If the path does not exist, the specified path of directories 
    is created.
Parameters:
- path: ...

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
