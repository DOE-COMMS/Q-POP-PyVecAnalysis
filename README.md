![GitHub release version](https://img.shields.io/github/v/release/DOE-COMMS/Q-POP-PyVecAnalysis?color=%2350C878&include_prereleases)
![License](https://img.shields.io/github/license/DOE-COMMS/Q-POP-PyVecAnalysis)
![GitHub Size](https://img.shields.io/github/repo-size/DOE-COMMS/Q-POP-PyVecAnalysis)
![HitCount](https://hits.dwyl.com/DOE-COMMS/Q-POP-PyVecAnalysis.svg?style=flat-square&show=unique)
![HitCount](https://img.shields.io/endpoint?url=https%3A%2F%2Fhits.dwyl.com%2FDOE-COMMS%2FQ-POP-PyVecAnalysis.json&label=total%20hits&color=pink)

# Q-POP-PyVecAnalysis (PyPINS)

**Py**thon **P**rogrammatic **I**dentification of **N**ovel **S**tructures (PyPINS) is a software suite written for identifying various flow structures that can exist in flow fields, such as vortices and eddys. In addition, this suite can be utilized to identify the location and number of novel topological structures like Skyrmions, Monopoles, and similar objects of interest to the condensed physics community. 

## Acknowledgment

Much work for this project was accomplished using the Extreme Science and Engineering Discovery Environment (XSEDE), which was supported by National Science Foundation grant number #1548562. Additionally, partial funding for this software was provided via a grant from the United States Department of Energy Basic Energy Sciences Grant Number: DE-SC0020145. Lastly, I would like to the thank the numerous research group members that also tested the program and suggested ideas for enhancment and usability of the software package.

## Identification Methods

Flow or Vector fields can be decomposed into many different ways. Decomposing the vector fields, whether into different vector fields or into scalar fields, can help researchers identify novel flow and topological structures that may exist in the fields. However, up until now, no such software exists to automate and programmatically identify these structures in a particular flow field. Therefore, using commonly available decomposition methods in the fluid dynamics literature, we have developed a software package for the post-processing of flow fields to understand the time evolution dynamics among other things. After the fields have been decomposed and properly thresholded, feature detection algorithms are applied to identify the interesting flow features that may exist. The list of identification methods is given below and more methods may be added in the future.

1. Vector Curl or Vorticity
2. Vector Divergence
3. Vector Kinetic Energy
4. Galilean Decomposition
5. Reynolds Decomposition
6. Q-Criterion Decomposition
7. Delta-Criterion Decomposition
8. Lambda-Criterion Decomposition
9. Winding Number
10. Topological Charge
11. Swirling Strength

These various identification methods can be applied to both user-specified analytical flow fields or discrete flow fields, such one that may come from Computational Fluid Dynamics Simulations or Experimental Data. Therefore allowing extendable access for the given problem at hand.

## License Information

PyPINS is licensed under GNU GPL v3. More information regarding the licensing can be found in the provided license.txt file that accompanies this software.

## Conflicts of interest

No known conflicts of interest exist for the author and this software. 
