# collapsar
The goal is to model the evolution of the accretion disk around a black hole (BH) during the core-collapse of a rapidly rotating massive star. The input of the model are the stellar density and rotational angular frequency profiles and the output is the accretion history of the black hole.

The user should first download all three pieces of ".py" codes: "collapsar.py" does the main computation, "run_collapsar.py" is a wrapper that handles the input and output, and "mplchange.py" is only used if one wants to plot the evolutionary history. After downloading the codes, the user needs to change the variable "fdir" in both "collapsar.py" and "run_collapsar.py" into the name of the directory which contains the MESA stellar profile file. Possible file output will be saved in that directory as well.

There are two ways of running the model.

(1) Type "Python run_collapsar.py" in the command line. The user needs to supply the list of profile file names as a list variable "fname_list". Then the code will go over each of the listed files and print the following numbers: final BH mass "Mbh", final BH spin "abh", total energy carried away by the disk wind "Ewind", total accretion energy released near the ISCO "Eacc" (including both neutrino emission and wind kinetic energy), accretion energy released near the ISCO when the innermost disk regions are in the ADAF regime "Eacc_ADAF" (only including wind kinetic energy from near the ISCO).

(2) Type "Python collapsar.py" in the command line. The user needs to supply a single profile file name as "fname". Then, the code will plot the evolutionary history of the accretion disk for this model and save it as a ".png" figure.

The model contains three free parameters: "s_PL" (between 0 and 1, default 0.5) is the power-law index of the radial accretion rate profile in the ADAF regime, "alp_ss" (should be between 0.01 and 0.1, default 0.03) is the viscosity parameter, "themin" (between 0 and 90 degrees, default 30 degrees) is the minimum polar angle (in radian) below which fallback is impeded by the mechanical feedback from the disk. There is another parameter "dt_over_tvis" (<< 1, default value 0.01), which is the ratio between the time step and the current viscous time.

The detailed description of the model can be found in the document "model_description.pdf". 
