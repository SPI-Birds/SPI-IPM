# SPI-IPM
Materials for setting up, running, analysing, and interpreting Integrated Population Models using SPI-Birds data.
Contains code for the different steps/analyses and a code manual providing background information and code documentation. The code manual is published via GitHub pages and accessible here: https://spi-birds.github.io/SPI-IPM/

The repository contains two main folders: "SPI-IPM_Code" and "SPI-IPM_CodeManual".
"SPI-IPM_Code" has subfolders that are labelled to match thematically to the chapters in the code manual:

- **01_Data_Reformatting**
    * DataPrep_PopID.R - *Transform data in SPI-Birds standard format into input for IPM analysis*
    <br/>
    
- **02-04_IPM_Setup&Run**
    * IPMRun_PopID.R - *Implement and run IPM for a number of populations*
    * IPMSetup.R - *Set up IPM for a specified population*
    * InitSim.R - *Simulate initial values*
    * dCJS_CustomDist.R - *Custom distribution for efficient CJS mark-recapture model*
    <br/>
    
- **05_Model_Assessment**
    * Comp_DataVSPred.R - *Compare IPM estimates to data*
    * Comp_IndepVSInteg.R - *Compare posterior estimates from integrated and independent analyses*
    * StochProjections.R - *Run stochastic population projections using IPM posteriors*
    <br/>
    
- **06_Results_Vizualization**
    * PlotResults_CovariateEffects.R - *Predict and plot vital rates as functions of environmental covariates*
    * PlotResults_General.R - *Plot vital rate and population size estimates (averages & over time)*
    * PlotResults_VRForestPlot.R - *Make a forest plot to summarise and compare vital rate estimates across populations*
    <br/>
    
- **07_Follow-up_Analyses**
    * LTRE_FixedDesign
        * fixedLTRE_Plots_ContTime.R - *Plot results from fixed-design LTRE*
        * fixedLTRE_PopID.R - *Run fixed-design LTRE for a pre-defined population (set up in fixedLTRE_RunAll.R)*
        * fixedLTRE_RunAll.R - *Set up for and run fixed-design LTRE for a number of populations*
        
    * LTRE_PeriodDesign
        * periodLTRE_Plots.R - *Plot results from period-design LTRE*
        * periodLTRE_PopID.R - *Run period-design LTRE for a pre-defined population (set up in periodLTRE_RunAll.R)*
        * periodLTRE_RunAll.R - *Set up for and run period-design LTRE for a number of populations*
        
    * LTRE_RandomDesign
        * randomLTRE_Plots_Cont.R - *Plot estimated contributions from random-design LTRE (excl. immigration)*
        * randomLTRE_Plots_ContImm.R - *Plot estimated contributions from random-design LTRE (incl. immigration)*
        * randomLTRE_Plots_SensElas.R - *Plot estimated sensitivities and elasticities*
        * randomLTRE_PopID.R - *Run random-design LTRE for a pre-defined population (set up in randomLTRE_RunAll.R)*
        * randomLTRE_RunAll.R - *Set up for and run random-design LTRE for a number of populations*
        
    * CrossPopCovariation.R - *Test for cross-population covariation in population sizes and vital rates*
    * TimeTrends_Test.R - *Test for time trends in estimated population sizes and vital rates*




The code manual is still under development (sections marked with "TBA" will be added over the course of 2022). 
If you would like assistance using (part of) the code or have other questions/inquiries about any of the materials, please contact Chloé R. Nater (chloe.nater@nina.no).

The latest release is v0.0.1 and is citable via Zenodo: [![DOI](https://zenodo.org/badge/364228589.svg)](https://zenodo.org/badge/latestdoi/364228589)
