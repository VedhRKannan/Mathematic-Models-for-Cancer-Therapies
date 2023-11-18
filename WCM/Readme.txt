Readme for:
A Mechanistic Pan-Cancer Pathway Model Informed by Multi-Omics Data Interprets Stochastic Cell Fate Responses to Drugs and Mitogens

Created by Mehdi Bouhaddou (December 2017)
Feel free to contact us regarding any questions related to model execution. We are happy to assist.

-------------------------------------------
To run the model you must first install Sundials CVODES interface for MATLAB (sundialsTB v.2.4.0). Sundials enables your Matlab model to be run on a c-based ODE solver, speeding up model execution.
Here are some useful instructions for installing Sundials CVODES for Windows 7:

1.	Install .NET Framework 4.0 from here: https://www.microsoft.com/en-us/download/details.aspx?id=17851
Make sure no other versions (prior or later) are already installed on your computer. If so, un-install and re-install this version (4.0).

2.	Install Microsoft Windows SDK for Windows 7 and .NET Framework from here: https://www.microsoft.com/en-us/download/details.aspx?id=8279

3.	Download Sundials 2.4.0 (https://computation.llnl.gov/projects/sundials/sundials-software). Follow Sundials documentation to properly install Sundials.
Generally, you must open SundialsTB folder and run install_STB script. 
If you get an error message saying “Error using mex”, one workaround we’ve identified is by following these additional steps:

	1) Navigate to sundials-2.4.0/sundialsTB/cvodes/cvm/src. There should be 4 files in this directory (cvm.c, cvm.h, cvmOpts.c, and cvmWrap.c).
	2) Find and replace all instances of “mxCreateScalarDouble” to “mxCreateDoubleScalar” in cvm.c and cvmWrap.c. (The other two files shouldn’t have any instances of mxCreateScalarDouble).

You must run startup_STB.m upon each Matlab startup to initialize Sundials.

If needed, there are some good troubleshooting tips here:
How do I install Microsoft Windows SDK 7.1?:
https://www.mathworks.com/matlabcentral/answers/101105-how-do-i-install-microsoft-windows-sdk-7-1


-------------------------------------------
Once Sundials CVODES has been installed, you can run the model as follows:

1. TO RUN THE MODEL: To run a single run of the model, run "RunModel.m".

For example, to run model for a deterministic 24-hour simulation of a single serum starved cell without any stimulus, run:

flagD=1; %1 for Deterministic Simulation, 0 for Stochastic Simulation
th=24; %Simulation Time (Hours)
STIM=zeros(775,1); %Stimulation vector (adds values to respective species)
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);

If you want to add some EGF and insulin, you would run it as follows:

flagD=1; %1 for Deterministic Simulation, 0 for Stochastic Simulation
th=24; %Simulation Time (Hours)
STIM=zeros(775,1); %Stimulation vector (adds values to respective species)
STIM(156:162)=[10,0,0,0,0,0,1721]; %Adding EGF and Insulin
[tout_all,xoutG_all,xoutS_all]=RunModel(flagD,th,STIM,[],[],[],[]);

Names Sheet:
- PARCDL tab has indices for each species. For example, EGF is index number 156 (symbol “E” in the model).
- TLC tab has indices for protein conglomerates (see Supplemental Methods for explanation).
- TL tab has indices for all genes in model.
- gm tab has indices for species fed into the gm model (gene switching and mRNA expression).
- d tab has indices for protein complexes or post-translationally modified species that are degraded via first-order degradation kinetics.
- kS tab explains the ordering of the “kS” variable, which holds all the rate constants. First column: how many constants in each sub-variable. Second column: Name of sub-variable. Third column: beginning index of sub-variable. Fourth column: end index of sub-variable.


2. SIMULATIONS: To reproduce simulations used for figures run "Run_Figures_Pub.m". This may take several days depending on the type of simulation. This script generates matfiles that are used to create the figures. We recommend running each section independently (make sure you run the “UBIQUITOUS CODE” section before running). Example matflies have been provided for making figures 2C and S3B.

3. FIGURES: To reproduce figures once matfiles (simulation results) are generated, run "Plot_Figures_Pub.m". This may take several minutes depending on figure. Again, we recommend running each section independently as above. Again, example matfiles have been provided for making figures 2C and S3B.

4. INITIALIZATION: To re-run initialization (files provided in initialized/ directory), run “RunInitialize.m”. See function header for usage details. Files in the "initialized/" directory are specific to MCF10A cells. 

5. TO RUN U87 VERSION OF MODEL: To run U87 cell initialized data, rename "initializedU87" to simply "initialized". You also need to uncomment line 9 in RunPrep.m, which will import the U87 genetic data instead of the MCF10A genetic data.
