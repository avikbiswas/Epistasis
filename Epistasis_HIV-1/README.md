#Potts model scan of double mutant cycles

Code and data for the Potts model analyses.

Models: The inferred Potts models for PR, RT, and IN from multiple sequence alignments derived from the Stanford HIV drug resistance database are provided in the data/models directory. For more details, please refere to the Materials and Methods section.

Data:
All relevant data are provided within the data directory.
alphabet_maps: contain the mapping from the full 21-letter amino acid alphabet to a reduced 4-letter alphabet for each protein. For more details, please refere to the Materials and Methods section.
ddEdata: contain the the Potts model double mutant cycle data for all 3 proteins in the reference nl4-3 sequence.
{PR, RT, IN}: relevant reference sequences in the full 21 letter amino acid alphabet, in the reduced 4-letter alphabet, and other relevant files used in the notebooks.

Code:
Python code to calculate the all possible double mutant cycles (in reduced 4 letter alphabet; refer to Methods) in a given reference sequence is provided in the directory Potts_double_mutant_cycle calculation, separately for PR, RT, and IN.

Additionally, Jupyter notebooks and related utility functions are provided. The notebooks can be used for data analysis as well as to generate figures presented in Figure 1. Run the notebooks from inside `notebooks/` .

```
Epistasis_HIV-1/
|-- notebooks/
|   |-- histogram_main.ipynb   # middle panels: ddE histograms for PR, RT, IN
|   |-- Circos_main.ipynb      # side panels: synergistic / antagonistic circos plots
|-- Potts_double_mutant_cycle_calculation/
|   |-- get_doublemutddE_all_{PR,RT,IN}.py   # compute ddE for all double mutant cycles
|   |-- README.md
|-- data/
|   |-- {IN,PR,RT}/
|   |   |-- {in,pr,rt}.nl43.seq                  # NL4-3 sequence, full alphabet
|   |   |-- {in,pr,rt}.nl43.reduce4.seq          # NL4-3 sequence, reduced alphabet
|   |   |-- {in,pr,rt}.consensus.seq             # consensus sequence, full alphabet
|   |   |-- {in,pr,rt}.consensus.reduce4.seq     # consensus sequence, reduced alphabet
|   |   |-- top_5_DRM.txt / bottom_5_DRM.txt     # dotted lines on histograms
|   |   |-- top_20_DRM.txt / bottom_20_DRM.txt   # highlighted chords in circos (Table S1)
|   |-- alphabet_maps/
|   |   |-- {IN,PR,RT}.alphamap.4letter.txt      # 21-letter -> 4-letter alphabet maps (reduced-alphabet dictionary)
|   |-- ddE_data/
|   |   |-- {in,pr,rt}.nl43.ddE.all.txt          # Potts ddE for all double mutants, NL4-3 reference
|   |   |-- {in,pr,rt}.consensus.ddE.all.txt     # Potts ddE for all double mutants, consensus reference
|   |-- models/
|   |   |-- J_PR.npy   # Potts couplings J for PR, shape (4851, 16)
|   |   |-- J_RT.npy   # Potts couplings J for RT, shape (17578, 16)
|   |   |-- J_IN.npy   # Potts couplings J for IN, shape (34453, 16)
|-- utilities/
    |-- functions.py
    |-- circ.py
```
