# LBS-comparison

This is the GitHub repository underlying the analysis for our manuscript "_Comparative analysis of methods for the prediction of protein-ligand binding sites_", which can be found [here](insert_ref).

## Dependencies

### Prediction methods

In this work, we compared the performance of eleven ligand binding site prediction methods, which are the following. Instructions to install them can be found in their respective repositories.

1. [VN-EGNN](https://github.com/ml-jku/vnegnn)  [[1]](https://arxiv.org/abs/2404.07194)
2. [IF-SitePred](https://github.com/annacarbery/binding-sites)  [[2]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00821-4)
3. [GrASP](https://github.com/tiwarylab/GrASP)  [[3]](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.3c01698)
4. [PUResNet](https://github.com/jivankandel/PUResNetV2.0)  [[4]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00865-6)
5. [DeepPocket](https://github.com/devalab/DeepPocket)  [[5]](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799)
6. P2Rank<sub>Cons</sub>  [[6]](https://academic.oup.com/nar/article/50/W1/W593/6591527)
7. [P2Rank](https://github.com/rdk/p2rank)  [[7]](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0285-8)
8. [fpocket](https://github.com/Discngine/fpocket) [[8]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-168)
9. [PocketFinder<sup>+</sup>](https://compbio.cs.princeton.edu/concavity/)  [[9]](https://pubmed.ncbi.nlm.nih.gov/15757999/)
10. [Ligsite<sup>+</sup>](https://compbio.cs.princeton.edu/concavity/)  [[10]](https://pubmed.ncbi.nlm.nih.gov/9704298/#:~:text=Using%20a%20set%20of%20receptor,of%20LIGSITE%20is%20its%20speed.)
11. [Surfnet<sup>+</sup>](https://compbio.cs.princeton.edu/concavity/)  [[11]](https://pubmed.ncbi.nlm.nih.gov/8603061/)

P2Rank<sub>Cons</sub> is the same programme as P2Rank, but passing an extra argument which contains evolutionary conservation information of the target protein.

Methods 9-11 represent the implementations by Capra, _et al._, 2009, [[12]](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000585) and are therefore indicated by the <sup>+</sup> superindex.

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
- [POVME](https://github.com/durrantlab/POVME) [[REF]](https://pubs.acs.org/doi/10.1021/ct500381c)
- [ProteinVolume](https://gmlab.bio.rpi.edu/PVolume.php)  [[REF]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0531-2)
- [freeSASA](https://freesasa.github.io/)  [[REF]](https://f1000research.com/articles/5-189/v1)
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)  [[REF]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7737788/)
- [PyMol](https://pymol.org/support.html)  [[REF]](https://legacy.ccp4.ac.uk/newsletters/newsletter40/11_pymol.pdf)

## References
