# OFCGq
This contains code for analyzing raw fiber photometry, RNA-sequencing, TissueCyte whole-brain images, blink tracking and and engram calcium activity, and wth linear mixed-effects modeling. Detailed experimental procedures and analysis methods are described in the Methods section of the associated manuscript.


--FiberPhotometry--
All analyses run on Python ver. 3.13.1 

Fig1f_adra1acKO_eLacco: 
Fiber photometry was used to measure extracellular lactate signals in orbitofrontal cortex (OFC) astrocytes from control mice and mice with conditional knockout (cKO) of Adra1a in OFC astrocytes.

Fig1h_iBARK_eLacco: 
Fiber photometry was used to measure extracellular lactate signals signal in OFC from control mice and mice with iβARK2 in OFC astrocytes.

Fig1j_hM3Dq_eLacco: 
Fiber photometry was used to measure extracellular lactate signals signal in OFC from control mice and mice with hM3Dq in OFC astrocytes.

Fig3a_hM3Dq_iLacco: 
Fiber photometry was used to measure neural intracellular lactate signals in OFC from control mice and mice with hM3Dq in OFC astrocytes.

Fig3b_iBARK_iLacco: 
Fiber photometry was used to measure neural intracellular lactate signals signal in OFC from control mice and mice with iβARK2 in OFC astrocytes.


--RNA-sequencing--
Contains batch scripts for raw FASTQ file processing and differential gene expression analysis in R using edgeR.
Analyses run on R ver. 4.5.3, RStudio ver. 2026.07.1, and command line


--WholeBrainAnalysis--
Includes R code used to analyze regional intensity values extracted from NeuroInfo.
Analyses run on R ver. 4.5.3 and RStudio ver.2026.07.1


--Peri-Blink Calcium--
Includes code for blink detection from DeepLabCut eyelid tracking and analysis of blink-aligned ΔF/F traces in Fos+ and Fos- cells. MATLAB ver. R2024b


--LMM--
Linear mixed-effects modeling (LMM) stats: Includes analysis of nested physiological data (Group, Mouse, and Cell hierarchies) using LMM. MATLAB ver. R2024b
