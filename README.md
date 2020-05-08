# Capillary-Kir-Model

Code to accompany the paper titled "Capillary Kir channel as sensor and amplifier of neuronal signals: modeling insights on  K+ -mediated neurovascular communication". Currently under review in PNAS. 

*Code is currently being prepared for upload and will be uploaded shortly.*

Notes:

1. Run script files within the "scripts" folder to reproduce each figure. Output .bmp images from each script are placed in appropriate figure-specific subfolders within the "output" folder. The various microvascular networks used in the code are .mat files in the "networks" folder. Networks are represented as graph objects using MATLAB's "Graph and Network Algorithms" functionality. Experimental data used in some figures are included in the "experimental_data" folder. Some simulation results are saved in the "saved_data" folder in order to save time when reproducing manuscript figures. External functions used are in the "functions" folder.

2. Networks used in the simulations were created using vectorized microvascular network data from the vibrissa primary sensory cortex in a mouse aquired from Blinder et al. In-house code was used to create smaller partitions of this network for the simulations. Blinder, P., Tsai, P. S., Kaufhold, J. P., Knutsen, P. M., Suhl, H., & Kleinfeld, D. (2013). The cortical angiome: an interconnected vascular network with noncolumnar patterns of blood flow. Nature neuroscience, 16(7), 889.

3. Functions used from the MATLAB Central File Exchange:
- Harald Hentschke (2020). abfload (https://www.mathworks.com/matlabcentral/fileexchange/6190-abfload), MATLAB Central File Exchange.
- Puck Rombach (2020). Largest Component (https://www.mathworks.com/matlabcentral/fileexchange/30926-largest-component), MATLAB Central File Exchange.
- Tim Franklin (2020). ODE Progress Bar and Interrupt (https://www.mathworks.com/matlabcentral/fileexchange/9904-ode-progress-bar-and-interrupt), MATLAB Central File Exchange. 

4. Some simulations have inherent stochasticity and/or randomness so they may not exactly match data reported in the manuscript. 

5. We have tested this code using MATLAB R2019a on a Windows 10 personal computer. For error reporting and suggestions please contact Arash Moshkforoush (amosh005@fiu.edu) . We welcome your comments and suggestions.

6. Arash Moshkforoush, Baarbod Ashenagar, and Nikolaos Tsoukias contributed to the development of this method.

7. If this work is used in any way, please cite the paper and provide appropriate acknowledgement.

