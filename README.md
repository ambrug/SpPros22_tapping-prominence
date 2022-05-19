# SpPros22_tapping-prominence

Analysis of tapping data and RPT judgments October 2021.
These data are from the 2nd experiment reported in Speech Prosody paper: Beware of the individual: Evaluating prominence perception in spontaneous speech.

* Study design: Petra Wagner, Leonie Schade, Anna Bruggeman
* Data collection & preparation: Leonie Schade, Anna Bruggeman
* Statistical Analysis: Anna Bruggeman

# Script files:
* 1_tap_vs_RPT.R Compares the tapping data with the RPT binary judgments. 
* 2_RFanalysis.R Computes random forests and visualises variable importances. Starts with preprocessing of acoustic data.

# Data files:
* extracted_acoustic.csv Contains acoustic data (Praat extracted) for each syllable in each stimulus
* pveloc.csv Contains tap duration and force data, binary RPT prom ratings, and pitch accentuation
* tapcorrect.txt List of stimulus/listener combinations that had correct number of taps
* tapwrong.txt List of stimulus/listener combinations that had correct number of taps
* wordinfo_manual.csv Lexical information about words and syllables 
* GECO_stimuli_overview.csv  Contextual info and info about presence of hesitations for each stimulus
