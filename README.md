# ana_phaseGain

Main analysis code for **Voloh, Oemisch, and Womelsdorf (2020) Phase of Firing Coding of Learning Variables across Prefrontal Cortex, Anterior Cingulate and Striatum during Feature Learning**

**Major functions**
- ana_glm_rate: wrapper to perform regression on spikes
- get_glm_cluster: wrapper to cluster on beta weights from regression
- plot_glm_rate_cluster: Figure 2
- get_glm_phaseBinned: wrapper to perform phase-based regression
- plot_glm_phaseBinned: Figure 3,4,5

**dependencies**
other than MATLAB toolboxs, you need:
- fieldtrip: http://www.fieldtriptoolbox.org/

**analysis steps**
1) Using fieldtrip, get a spike-triggered spectrum structure by calling ft_spiketriggeredspectrum
  - make sure to add some fields: area_lfp (are of each LFP channel), area (area of spike channel), iso (isolation quality of cell).
2) Set paths to code and data using *set_ana_paths.m*
3) Run *ana_glm_rate*
4) Run *get_glm_cluster*
  - plot results using *plot_glm_rate_cluster*
5) Run *get_glm_phaseBinned*
  - plot results using *plot_glm_phaseBinned*
