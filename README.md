# Capillary-Kir-Model

Code to accompany the paper titled "Capillary Kir channel as sensor and amplifier of neuronal signals: modeling insights on  K+ -mediated neurovascular communication". Currently under review in PNAS. 

Notes:

1. Run script files within the "scripts" folder to reproduce each figure. The various microvascular networks used in the code are .mat files in the "networks" folder within each script subfolder. Networks are represented as graph objects using MATLAB's "Graph and Network Algorithms" functionality. Some simulation results are saved in the "saved_data" folder within each script subfolder in order to save time when reproducing manuscript figures. Functions used are in the "functions" folder. Some simulations have inherent stochasticity and/or randomness so they may not exactly match data reported in the manuscript, however they should follow the same trend. 

2. Networks used in the simulations for figures 5, 6, and 7 were created using vectorized microvascular network data from the vibrissa primary sensory cortex of a mouse aquired from Blinder et al. Nat. Neuro. 2013 (https://www.nature.com/articles/nn.3426). In-house code was used to create smaller partitions of this network for the simulations. 

3. Functions used from the MATLAB Central File Exchange:
- Harald Hentschke (2020). abfload (https://www.mathworks.com/matlabcentral/fileexchange/6190-abfload), MATLAB Central File Exchange.
- Puck Rombach (2020). Largest Component (https://www.mathworks.com/matlabcentral/fileexchange/30926-largest-component), MATLAB Central File Exchange.
- Tim Franklin (2020). ODE Progress Bar and Interrupt (https://www.mathworks.com/matlabcentral/fileexchange/9904-ode-progress-bar-and-interrupt), MATLAB Central File Exchange. 

4. We have tested this code using MATLAB R2019a on a Windows 10 personal computer. For error reporting and suggestions please contact Arash Moshkforoush (amosh005@fiu.edu) . We welcome your comments and suggestions.

5. Modifying the file structure of the repository while using it may break some of the code.

6. Arash Moshkforoush, Baarbod Ashenagar, and Nikolaos Tsoukias contributed to the development of this method.

7. If this work is used in any way, please cite the paper and provide appropriate acknowledgement.

