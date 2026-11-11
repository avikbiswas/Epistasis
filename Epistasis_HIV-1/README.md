# Figure 1: Potts model scan of double mutant cycles

Code and data for Figure 1 (panels B-D) in `Epistasis Manuscript`.
Panel A (Potts model schematic) is not included here.

Run the notebooks from inside `notebooks/` .

```
epistasis_HIV-1_figure_1/
|-- notebooks/
|   |-- histogram_main.ipynb   # middle panels: ddE histograms for PR, RT, IN
|   |-- Circos_main.ipynb      # side panels: synergistic / antagonistic circos plots
|-- data/
|   |-- {IN,PR,RT}/
|   |   |-- kn.her2.all                # Potts ddE for all double mutants
|   |   |-- {in,pr,rt}.reduce4.redux   # reduced-alphabet dictionary
|   |   |-- top_5_DRM.txt / bottom_5_DRM.txt    # dotted lines on histograms
|   |   |-- top_20_DRM.txt / bottom_20_DRM.txt  # highlighted chords in circos (Table S1)
|   |-- all_data/
|   |   |-- PR_ddE.txt   # Potts ddE for all PR double mutants
|   |   |-- RT_ddE.txt   # Potts ddE for all RT double mutants
|   |   |-- IN_ddE.txt   # Potts ddE for all IN double mutants
|   |   |-- {pr,rt,in}.reduce4.redux   # reduced-alphabet maps
|   |-- models/
|   |   |-- J_PR.npy   # Potts couplings J for PR, shape (4851, 16)
|   |   |-- J_RT.npy   # Potts couplings J for RT, shape (17578, 16)
|   |   |-- J_IN.npy   # Potts couplings J for IN, shape (34453, 16)
|   |-- {pr,rt,in}.nl43.seq           # NL4-3 sequence, full alphabet
|   |-- {pr,rt,in}.nl43.reduce4.seq   # NL4-3 sequence, reduced alphabet
|-- output/                    # {IN,PR,RT}_histogram.png, {IN,PR,RT}_{syn,antag}.png
|-- utilities/
    |-- functions.py
    |-- circ.py
```

