To reproduce the main figures of the study, follow these steps:

1. Python Version and Package Requirements:

    Make sure you have Python 2.7 installed on your system.  
    Ensure you have the following packages installed: numpy, scipy, matplotlib, pickle, sklearn, and joblib.
    
2. Data:

    The data/ folder contains the following files:\
    - asymm_ncd_thr_N_74.mat: Directed structural connectivity data.\
    - dcmECfc_singleMouse_sim62131.mat: single subject EC without structural constraint.\
    - dcmECfc_singleMouse_sim92131.mat: single subject EC with 60% structural threshold of strongest kept links.\
    - dcmECfc_singleMouse_sim102131.mat: single subject EC with 40% structural threshold of strongest kept links.\
    - dcmECfc_singleMouse_sim82131.mat: single subject EC with 20% structural threshold of strongest kept links.

3. Running the Code:

    Execute the main.py script in the repository to reproduce the main figures of the study.  
    The script will load the data from the data/ folder and perform the necessary computations.
    
4. Results and Figures:

    The results will be saved in the figure/ folder within the repository.  
    For each nSIM value (62131, 92131, 102131, 82131), two folders will be saved:\
    - ec_sc/: which figures representing the EC-SC coupling for the specific nSIM.\
    - sc_ec/: which figures representing the SC-EC coupling for the specific nSIM.\  
    Specific figures from the manuscript (Figs 2,3) are saved as:  
    - in_out_coupling_perNetwork_bar_iso.svg\
    - incomingVSoutgoing_smaller.svg

5. Interpretation:
    
    The figures saved in the figure/ folder will display the main results of the study.  
    Refer to the manuscript to interpret and understand the meaning of the figures.

Note: If you encounter any issues or errors during execution, ensure that you have correctly set up the Python environment with the required packages and that the data files are in the right location.\
Also, make sure that you are using Python 2.7 to run the code.

