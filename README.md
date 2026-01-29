# Hierarchy and MIND in schizophrenia spectrum disorders (SSD)

This repository contains code and data created in support of the project García-San-Martín, N.; Bethlehem, R. AI; Sebenius, I. et al. *"Long-term morphometric similarity gradients relate to cortical hierarchy and psychiatric symptoms in schizophrenia"* medRxiv (2026). All code was written in Matlab, R, and Python. Files are described below.

## **Data**

- The [`brain_maps`](Data/brain_maps/) folder contains the hierarchical maps employed in this study, which were obtained using *neuromaps* toolbox, available at https://github.com/netneurolab/neuromaps.
  - The `sensorimotor-association (S-A) regional ranking` file was derived from https://doi.org/10.1016/j.neuron.2021.06.016.
  - The map of `evolutionary expansion`, which represents the ratio of the surface area in humans to that of macaque, was derived from https://www.pnas.org/doi/full/10.1073/pnas.1001229107.
  - The `functional gradient` corresponds to the first gradient of functional connectivity derived from diffusion map embedding, and was obtained from https://www.pnas.org/doi/10.1073/pnas.1608282113.
  * The `functional and structural epicenters in SCZ`, were obtained from an ENIGMA study available at https://www.nature.com/articles/s41380-024-02442-7.
    
- The `perm_sphere_10000_DK.csv` file corresponds to the permutation matrix employed for the statistical analyses (spin test).

## **Code**

The [`Code`](Code/) folder contains all the code required for running the analyses and generate data and figures. All scripts are designed to be sequentially executed. Don't forget to change the location variable regularly. 

-	The code used to compute cortical MIND networks is available at https://github.com/isebenius/MIND and corresponds to [MIND_long_01_MIND.py](Code/MIND_long_01_MIND.py). It computes MIND networks from FreeSurfer directory (by default stored in the surf/ folder) and returns a .csv file for each individual.
  
-	[MIND_long_02_COMBATLS.R](Code/MIND_long_02_COMBATLS.R) - applies the COMBATLS site harmonization method to MIND networks.
  
-	[MIND_long_03_degree_and_edges_PAFIP_COMBATLS.m](Code/MIND_long_03_degree_and_edges_PAFIP_COMBATLS.m) – calculates the regional MIND degree for each HC and SSD individual.
    
-	[MIND_long_04_gradients.m](Code/MIND_long_04_gradients.m) – calculates regional MIND gradients for each HC and SSD individual, generating brain maps of cortical and subcrotical degrees, as well as cortical gradients. The parameters for the gradient decoposition were the following:

```matlab
gm_ref = GradientMaps('approach', 'dm', 'kernel', 'normalized_angle');

```

- [MIND_long_05_longitudinal_dx_MIND.m](Code/MIND_long_05_longitudinal_dx_MIND.m) - performs multiple regression and linear mixed modeling (LMM) on MIND networks.

- [MIND_long_06_longitudinal_BPRS_MIND.m](Code/MIND_long_05_longitudinal_BPRS_MIND.m) - performs generalized linear modeling (GLM) and generalized LMM (GLMM) on psychiatric symptoms and MIND networks. 

-	[MIND_long_07_hierarchy_SCZ_brain_maps.m](Code/MIND_long_07_hierarchy_SCZ_brain_maps.m) – generates regional brain maps of cortical hierarchy and schizophrenia (SCZ) epicenters, and their correlation (co-localization) with MIND associations.

-	[MIND_long_08_cortical_MIND_association_maps.m](Code/MIND_long_08_cortical_MIND_association_maps.m) – generates regional brain maps of cortical MIND associations.

-	[MIND_long_09_subcortical_MIND_association_maps.m](Code/MIND_long_09_subcortical_MIND_association_maps.m) – generates regional brain maps of subcortical MIND associations.


## **License**

This project is licensed under the terms of the [GNU General Public License v3.0 license](LICENSE).


## **Cite**
If you use this software, please cite the following paper and software:

- García-San-Martín, N., Bethlehem, R.A., Sebenius, I. et al. Long-term morphometric similarity gradients relate to cortical hierarchy and psychiatric symptoms in schizophrenia. medRxiv (2026). https://doi.org/
  
- Natalia-García-San-Martín. NeuroimagingBrainNetworks/HierarchyLongitudinalGradientsPsychosis: v1.0.0-alpha. Zenodo https://doi.org/ (2026).
