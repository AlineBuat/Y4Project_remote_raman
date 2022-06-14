# Y4Project_public
Final Year MEng Biomedical Engineering Individual Project code

Project title "Detecting Bioweapons with Standoff Raman Spectroscopy"
These are code and files used during the project and used to generate the figures in the Final Report.

This repository contains the following folders:
- the code used to analyse data
- Arduino Control: arduino programs used to send positional data to the DACs, used to control Galvo Mirrors. Arduino programs to read data from a Quadrant photodiode, and translate it into positional information. program to Send and receive packets using NRF24L01 transceiver modules on Arduino Nano (transmitter) and Arduino Due (Receiver).
- EagleShield: Arduino Shield layout (schematic, board, project files) to control 2 galvo mirrors
- Raman: spectral data (.csv)  in folders. MATLAB programs to analyse the spectral data (plotRaman_20220601.m). requires the "BEADS" MATLAB package to run.
- BeamQuality: analyse images containing laser spots as imaged by a camera evaluated along a distance around the focal point, to find the spot size. Finally, plot spot size vs distance around the focal point to find M squared quality of the beam.
