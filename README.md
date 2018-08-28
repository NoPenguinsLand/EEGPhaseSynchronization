# EEGPhaseSynchronization
EEG Phase Synchronization: Rotation Spring 2018
This code computes phase-locking value (PLV), phase-lag index (PLI) in both time (Stam 2007) and frequency domain (Vinck 2011), coherency (C), and weighted phase-lag index (Vinck 2011). 
Since these codes perform EEG data analysis, raw data is needed. Please ask me for EEG data.

Contact info: mhn at bu dot edu

ELI5: EEG (electroencephalogram) is a promising recording technology that will help researchers and medical professionals better understand the activity of the brain in vivo. However, EEG recording data is only as useful as far as its interpretation goes. One particular feature we're looking for is how similar the electrical signals from different parts of the brain are, (which could mean many things but nonetheless still warrants further investigation). However, the electrical signals propagate throughout the brain and can be picked up by sensors far from the sources of the electrical signals. This is known as volume-conduction. (The sensors are located on the scalp of the brain but the sources are deep within the brain). Thus this provides serious limitation in understanding EEG data. Several different scientists proposed different metrics for measuring synchronization between electrical signals from different sensors while rejecting signals due to volume conduction and noise. This code computes and visualizes different metrics, which are PLV, PLI, C, and WPLI.
