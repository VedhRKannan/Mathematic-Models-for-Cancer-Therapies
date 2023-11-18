Readme for:
Integrating Multi-Omics Data with Mechanistic Pan-Cancer Pathway Models to Predict Stochastic Cell Fate

Created by Mehdi Bouhaddou (October 2017)

To run the model you must install Sundials CVODES interface for MATLAB (sundialsTB v.2.4.0).

1. SIMULATIONS: To reproduce simulations used for figures run "Run_Figures_Pub.m". This may take several days depending on the type of simulation. This script generates matfiles that are used to create the figures. We recommend running each section independently (make sure you run the “UBIQUITOUS CODE” section before running). Example mayflies have been provided for making figures 2C and S3B.

2. FIGURES: To reproduce figures once matfiles (simulation results) are generated, run "Plot_Figures_Pub.m". This may take several minutes depending on figure. Again, we recommend running each section independently as above.

3. TO RUN THE MODEL: To run a single run of the model, run "RunModel.m".
For example, to run model for a 24 hour simulation of a single serum starved cell, run:
flagD=1; %1 for Deterministic Simulation, 0 for Stochastic Simulation
th=24; %Simulation Time (Hours)
STIM=zeros(775,1); %Stimulation vector (adds values to respective species)
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);

4. INITIALIZATION: To re-run initialization (files provided in initialized/ directory), run “RunInitialize.m”. See function header for usage details.

