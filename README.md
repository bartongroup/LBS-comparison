# LBS-comparison

This is the GitHub repository underlying the analysis for our paper "_Comparative analysis of methods for the prediction of protein-ligand binding sites_", published in Journal of Cheminformatics and found [here](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00923-z).

[![DOI](https://zenodo.org/badge/826254471.svg)](https://zenodo.org/doi/10.5281/zenodo.13171100)

In this work, we gather 13 ligand binding site predictors, spanning 30 years, focusing on the latest machine learning based methods, such as VN-EGNN, IF-SitePred, GrASP, PUResNet, and DeepPocket and comparing them to established methods such as P2Rank, PRANK and fpocket and earlier methods such PocketFinder, Ligsite and Surfnet. We compare these thirteen methods thoroughly and benchmark them against our curated reference dataset, LIGYSIS, to perform an objective assessment of their prediction capabilities. An informed ranking of the methods, as well as a series of reflections and guidelines to advance this field result as conclusions of this analysis, which represents the most thorough analytical comparison of ligand binding site prediction methods to date, offering a clear framework for future developments in the field of ligand binding site prediction.

![This is the Graphical Abstract](./GRAPHICAL_ABSTRACT.png)


## Dependencies

### Prediction methods

In this work, we compared the performance of eleven ligand binding site prediction methods, which are the following. Instructions to install them can be found in their respective repositories.

1. [VN-EGNN](https://github.com/ml-jku/vnegnn)  [[1]](https://arxiv.org/abs/2404.07194)
2. [IF-SitePred](https://github.com/annacarbery/binding-sites)  [[2]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00821-4)
3. [GrASP](https://github.com/tiwarylab/GrASP)  [[3]](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.3c01698)
4. [PUResNet](https://github.com/jivankandel/PUResNetV2.0)  [[4]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00865-6)
5. [DeepPocket<sub>RESC</sub>](https://github.com/devalab/DeepPocket)  [[5]](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799)
6. [DeepPocket<sub>SEG</sub>](https://github.com/devalab/DeepPocket) [[5]](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799)
7. P2Rank<sub>CONS</sub>  [[6]](https://academic.oup.com/nar/article/50/W1/W593/6591527)
8. [P2Rank](https://github.com/rdk/p2rank)  [[7]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0285-8)
9. fpocket<sub>PRANK</sub> [[8]](https://link.springer.com/article/10.1186/s13321-015-0059-5)
10. [fpocket](https://github.com/Discngine/fpocket) [[9]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-168)
11. [PocketFinder<sup>+</sup>](https://compbio.cs.princeton.edu/concavity/)  [[10]](https://pubmed.ncbi.nlm.nih.gov/15757999/)
12. [Ligsite<sup>+</sup>](https://compbio.cs.princeton.edu/concavity/)  [[11]](https://pubmed.ncbi.nlm.nih.gov/9704298/#:~:text=Using%20a%20set%20of%20receptor,of%20LIGSITE%20is%20its%20speed.)
13. [Surfnet<sup>+</sup>](https://compbio.cs.princeton.edu/concavity/)  [[12]](https://pubmed.ncbi.nlm.nih.gov/8603061/)

P2Rank<sub>CONS</sub> is the same programme as [P2Rank](https://academic.oup.com/nar/article/50/W1/W593/6591527), but passing an extra argument which contains evolutionary conservation information of the target protein. This can be done through the [PRANK web server]() or locally after deployinw P2Rank with Docker. Instructions found [here](https://github.com/cusbg/p2rank-framework/wiki/P2Rank-deploy-with-Docker). 

fpocket<sub>PRANK</sub> are [fpocket](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-168) prdictions re-scored with [PRANK](https://link.springer.com/article/10.1186/s13321-015-0059-5).

DeepPocket<sub>RESC</sub> corresponds to the default [DeepPocket](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799) implementation, i.e., re-scored and re-ranked fpocket predictions. DeepPocket<sub>SEG</sub> represents the pocket shapes extracted from the [DeepPocket](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799) CNN segmentation module. These predictions were obtained by removing the `-r  3` flag, and so pocket shapes were extracted for all pockets, and not just the top-3.

Methods 11-13 represent the implementations by Capra, _et al._, 2009, [[13]](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000585) and are therefore indicated by the <sup>+</sup> superindex.

### Analysis

To run the analysis code, install the conda environment with the following dependencies:

`conda env create -f environment.yml`

or

`conda create -n lbs_comp_env python=3.10 scipy scikit-learn biopython pandas seaborn matplotlib numpy -c conda-forge`

#### Standard Python libraries
- [Scipy](https://scipy.org/)
- [Scikit-learn](https://scikit-learn.org/stable/)
- [Biopython](https://biopython.org/) 
- [Pandas](https://pandas.pydata.org/) 
- [Seaborn](https://seaborn.pydata.org/) 
- [Matplotlib](https://matplotlib.org/) 
- [Numpy](https://numpy.org/)


#### Additional dependencies

ProIntVar is used for some of the data processing. For installation instructions, refer here: [ProIntVar repository](https://github.com/bartongroup/ProIntVar/tree/JSU_branch).

The following are dependencies that were used for preotein and ligand site characterisation and visualisation:
- [POVME](https://github.com/durrantlab/POVME) [[14]](https://pubs.acs.org/doi/10.1021/ct500381c)
- [ProteinVolume](https://gmlab.bio.rpi.edu/PVolume.php)  [[15]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0531-2)
- [freeSASA](https://freesasa.github.io/)  [[16]](https://f1000research.com/articles/5-189/v1)
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)  [[17]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7737788/)
- [PyMol](https://pymol.org/support.html)  [[18]](https://legacy.ccp4.ac.uk/newsletters/newsletter40/11_pymol.pdf)

## References

1. Sestak, F., et al., VN-EGNN: E(3)-Equivariant Graph Neural Networks with Virtual Nodes Enhance Protein Binding Site Identification. [arXiv](https://arxiv.org/abs/2404.07194) [cs.LG], 2024.
2. Carbery A, Buttenschoen M, Skyner R, von Delft F, Deane CM. Learnt representations of proteins can be used for accurate prediction of small molecule binding sites on experimentally determined and predicted protein structures. J Cheminform. 2024 Mar 14;16(1):32. doi: [10.1186/s13321-024-00821-4](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00821-4). PMID: 38486231; PMCID: PMC10941399.
3. Smith Z, Strobel M, Vani BP, Tiwary P. Graph Attention Site Prediction (GrASP): Identifying Druggable Binding Sites Using Graph Neural Networks with Attention. J Chem Inf Model. 2024 Apr 8;64(7):2637-2644. doi: [10.1021/acs.jcim.3c01698](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01698). Epub 2024 Mar 7. PMID: 38453912; PMCID: PMC11182664.
4. Jeevan K, Palistha S, Tayara H, Chong KT. PUResNetV2.0: a deep learning model leveraging sparse representation for improved ligand binding site prediction. J Cheminform. 2024 Jun 7;16(1):66. doi: [10.1186/s13321-024-00865-6](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00865-6). PMID: 38849917; PMCID: PMC11157904.
5. Aggarwal R, Gupta A, Chelur V, Jawahar CV, Priyakumar UD. DeepPocket: Ligand Binding Site Detection and Segmentation using 3D Convolutional Neural Networks. J Chem Inf Model. 2022 Nov 14;62(21):5069-5079. doi: [10.1021/acs.jcim.1c00799](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799). Epub 2021 Aug 10. PMID: 34374539.
6. Jendele L, Krivak R, Skoda P, Novotny M, Hoksza D. PrankWeb: a web server for ligand binding site prediction and visualization. Nucleic Acids Res. 2019 Jul 2;47(W1):W345-W349. doi: [10.1093/nar/gkz424](https://academic.oup.com/nar/article/47/W1/W345/5494740). PMID: 31114880; PMCID: PMC6602436.
7. Krivák R, Hoksza D. P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. J Cheminform. 2018 Aug 14;10(1):39. doi: [10.1186/s13321-018-0285-8](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0285-8). PMID: 30109435; PMCID: PMC6091426.
8. Krivák R, Hoksza D. Improving protein-ligand binding site prediction accuracy by classification of inner pocket points using local features. J Cheminform. 2015 Apr 1;7:12. doi: [10.1186/s13321-015-0059-5](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0059-5). PMID: 25932051; PMCID: PMC4414931.
9. Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand pocket detection. BMC Bioinformatics. 2009 Jun 2;10:168. doi: [10.1186/1471-2105-10-168](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-168). PMID: 19486540; PMCID: PMC2700099.
10. An J, Totrov M, Abagyan R. Pocketome via comprehensive identification and classification of ligand binding envelopes. Mol Cell Proteomics. 2005 Jun;4(6):752-61. doi: [10.1074/mcp.M400159-MCP200](https://www.sciencedirect.com/science/article/pii/S1535947620314742?via%3Dihub). Epub 2005 Mar 9. PMID: 15757999.
11. Hendlich M, Rippmann F, Barnickel G. LIGSITE: automatic and efficient detection of potential small molecule-binding sites in proteins. J Mol Graph Model. 1997 Dec;15(6):359-63, 389. doi: [10.1016/s1093-3263(98)00002-3](https://www.sciencedirect.com/science/article/pii/S1093326398000023?via%3Dihub). PMID: 9704298.
12. Laskowski RA. SURFNET: a program for visualizing molecular surfaces, cavities, and intermolecular interactions. J Mol Graph. 1995 Oct;13(5):323-30, 307-8. doi: [10.1016/0263-7855(95)00073-9](https://www.sciencedirect.com/science/article/pii/0263785595000739?via%3Dihub). PMID: 8603061.
13. Capra JA, Laskowski RA, Thornton JM, Singh M, Funkhouser TA. Predicting protein ligand binding sites by combining evolutionary sequence conservation and 3D structure. PLoS Comput Biol. 2009 Dec;5(12):e1000585. doi: [10.1371/journal.pcbi.1000585](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000585). Epub 2009 Dec 4. PMID: 19997483; PMCID: PMC2777313.
14. Durrant JD, Votapka L, Sørensen J, Amaro RE. POVME 2.0: An Enhanced Tool for Determining Pocket Shape and Volume Characteristics. J Chem Theory Comput. 2014 Nov 11;10(11):5047-5056. doi: [10.1021/ct500381c](https://pubs.acs.org/doi/10.1021/ct500381c). Epub 2014 Sep 29. PMID: 25400521; PMCID: PMC4230373.
15. Chen CR, Makhatadze GI. ProteinVolume: calculating molecular van der Waals and void volumes in proteins. BMC Bioinformatics. 2015 Mar 26;16(1):101. doi: [10.1186/s12859-015-0531-2](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0531-2). PMID: 25885484; PMCID: PMC4379742.
16. Mitternacht S. FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Res. 2016 Feb 18;5:189. doi: [10.12688/f1000research.7931.1](https://f1000research.com/articles/5-189/v1). PMID: 26973785; PMCID: PMC4776673.
17. Pettersen EF, Goddard TD, Huang CC, Meng EC, Couch GS, Croll TI, Morris JH, Ferrin TE. UCSF ChimeraX: Structure visualization for researchers, educators, and developers. Protein Sci. 2021 Jan;30(1):70-82. doi: [10.1002/pro.3943](https://onlinelibrary.wiley.com/doi/10.1002/pro.3943). Epub 2020 Oct 22. PMID: 32881101; PMCID: PMC7737788.
18. Schrödinger, L.L.C., The PyMOL Molecular Graphics System, Version 1.8. 2015. Available at [https://pymol.org/](https://pymol.org/)
