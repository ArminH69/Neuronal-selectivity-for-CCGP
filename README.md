\# Neuronal Selectivity and Geometric Alignment in Hippocampal Abstract Generalization



This repository contains the analysis code associated with the manuscript:



\*\*“Neuronal selectivity and geometric alignment in the human hippocampus support abstract generalization.”\*\*



The study investigates how neuronal selectivity and population-level representational organization contribute to abstract generalization. The analysis consists of two complementary parts: a computational toy model and a re-analysis of previously published human hippocampal single-unit recordings.



\## Repository Structure



```text

.

├── toy\_model/

│   └── toy\_model\_analysis.py

│

├── hippocampal\_analysis/

│   ├── Category\_neurons\_main\_code.m

│   ├── Identity\_neurons\_main\_code.m

│   ├── build\_pseudopopulation\_independent.m

│   ├── compute\_category\_axis\_alignment.m

│   ├── compute\_category\_CCGP.m

│   ├── compute\_condition\_means.m

│   ├── compute\_context\_CCGP.m

│   ├── cosine\_similarity.m

│   ├── create\_independent\_trial\_split.m

│   ├── find\_min\_K\_independent.m

│   ├── getPredictorP.m

│   ├── sample\_neurons.m

│   ├── scenario\_encoding\_summary.m

│   └── train\_test\_decoder.m

│

└── README.md

```



## Toy model analysis

The `toy_model` directory contains the Python implementation of the computational
toy model used to investigate how neuronal population composition influences
abstract representation and cross-condition generalization.

The model generates artificial neural populations with different selectivity
profiles, including stimulus-, response-, context-, noise-, and
category-selective neurons. Population composition is systematically manipulated
to test how different neuronal subpopulations influence task-variable encoding
and AC-vs-BD cross-condition generalization performance (CCGP).

The analysis includes:

- generation of heterogeneous mixed-selectivity neural populations;
- comparison of different population composition scenarios;
- linear encoding analysis of task variables;
- cross-condition decoding analysis (CCGP);
- category-neuron population sweep, where the number of category-selective
  neurons is varied while other neuronal populations are kept constant.

The main script:

`Toy_model_analysis.py`

generates all toy-model simulations, saves numerical results, and produces the
figures associated with the category-neuron sweep.

\### Python requirements

The toy-model analysis requires:


\* Python 3

\* NumPy

\* pandas

\* Matplotlib

\* scikit-learn


\## Human Hippocampal Analysis


The `hippocampal_analysis` directory contains the MATLAB code used to test whether the principles identified in the computational toy model also apply to biological neural populations. Using previously published human hippocampal single-neuron recordings, the analysis examines whether abstract cross-context generalization depends simply on increasing the amount of stimulus-related information in the population, or instead on the selective contribution of neurons whose response profiles reflect the latent AC-vs-BD category structure. Hippocampal neurons are first characterized according to their task-variable selectivity, and stimulus-related neurons are further separated into **category-like** and **identity-like** populations using an independent subset of trials. Controlled pseudopopulations are then constructed by progressively increasing either category-like or identity-like neurons while keeping the other neuronal classes fixed. This design allows direct comparison of how the two neuronal populations influence stimulus encoding, AC-vs-BD CCGP, and population representational geometry. Additional analyses in the category-like manipulation examine whether improved generalization is associated primarily with stronger category separation or with more consistent alignment of the category axis across contexts, and whether changes in geometric alignment statistically mediate the relationship between category-like neuronal composition and CCGP.

The analysis is organized around two main scripts:


The `hippocampal_analysis` directory contains the MATLAB code used to analyze
the human hippocampal single-neuron data.

The analysis is organized around two main scripts:

### 1. `Category_neurons_main_code.m`

This is the main analysis for examining the effect of progressively adding
**category-like neurons** to the hippocampal population while keeping the
identity-like population fixed.

The script performs the complete category-neuron analysis pipeline, including:

1. Single-neuron encoding analysis for stimulus, context, response, and reward.
2. Functional classification of hippocampal neurons.
3. Merging single-variable and dominant-selectivity neuronal classes.
4. Independent splitting of trials into classification and population-analysis sets.
5. Classification of stimulus-selective neurons as category-like or identity-like.
6. Definition of the neuronal pools used for pseudopopulation construction.
7. Construction and repeated sampling of real-data population scenarios with
   progressively increasing numbers of category-like neurons.
8. Summary of encoding and CCGP across population scenarios.
9. AC-vs-BD CCGP analysis and visualization.
10. Comparison of task-variable encoding across population scenarios.
11. PCA-based representational geometry analysis and quantification of AC-vs-BD
    category separation.
12. Quantification of category-axis alignment and category-axis strength across
    independent pseudopopulation realizations.
13. Summary of geometric measures across neuronal-composition scenarios.
14. Correlation of category-axis alignment and category-axis strength with CCGP.
15. Partial-correlation analysis controlling for neuronal composition and the
    alternative geometric measure, including FDR correction.
16. Visualization of the relationships between geometric measures and CCGP.
17. Regression analyses relating neuronal composition, population geometry,
    and CCGP.
18. Mediation analysis testing whether category-axis alignment mediates the
    relationship between category-like neuron number and CCGP.
19. Bootstrap estimation of the mediation effect using 5,000 resamples.

The population manipulation increases the number of category-like neurons
(0, 5, 10, and 13) while keeping the identity-like population fixed at
13 neurons. Other neuronal populations are held constant across scenarios.

### 2. `Identity_neurons_main_code.m`

This script performs the complementary analysis for **identity-like neurons**.

The category-like population is kept fixed at 13 neurons while the number of
identity-like neurons is progressively increased (0, 5, 10, and 13).

The script includes:

1. Single-neuron encoding analysis.
2. Functional classification of hippocampal neurons.
3. Independent trial splitting.
4. Classification of stimulus-selective neurons into category-like and
   identity-like subtypes.
5. Construction of controlled pseudopopulation scenarios with increasing
   numbers of identity-like neurons.
6. Encoding and AC-vs-BD CCGP analyses across population scenarios.
7. Summary and visualization of encoding and CCGP changes.
8. PCA-based representational geometry analysis and quantification of
   AC-vs-BD category separation.


### Supporting functions

The remaining MATLAB files are supporting functions called by the two main
analysis scripts:

- `create_independent_trial_split.m` – creates independent trial sets for
  neuronal subtype classification and subsequent population analyses.
- `find_min_K_independent.m` – determines the maximum common number of
  held-out trials that can be sampled across neurons and task conditions.
- `build_pseudopopulation_independent.m` – constructs pseudopopulation activity
  matrices using only held-out analysis trials.
- `sample_neurons.m` – randomly samples neurons from the specified neuronal pool.
- `scenario_encoding_summary.m` – summarizes encoding strength for the neurons
  included in each population realization.
- `compute_category_CCGP.m` – computes AC-vs-BD cross-context generalization.
- `compute_context_CCGP.m` – computes cross-condition generalization of context.
- `train_test_decoder.m` – trains and evaluates the population decoder.
- `compute_condition_means.m` – calculates population responses for the eight
  stimulus × context conditions.
- `compute_category_axis_alignment.m` – calculates category-axis alignment and
  strength across contexts.
- `cosine_similarity.m` – computes cosine similarity between population axes.
- `getPredictorP.m` – evaluates the significance of predictors in the
  single-neuron encoding models.

\## Data


The dataset will be referenced through the original publication after manuscript submission.


\## Manuscript Status



The manuscript associated with this repository is currently in preparation/submission.



Citation information and a link to the manuscript will be added when available.



\## Contact



For questions regarding the code or analyses, please open an issue in this repository.



