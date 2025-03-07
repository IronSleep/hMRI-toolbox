# The hMRI-toolbox

A toolbox for quantitative MRI and *in vivo* histology using MRI (hMRI).

## Background

Neuroscience and clinical researchers are increasingly interested in 
quantitative magnetic resonance imaging (qMRI) 
due to its sensitivity to micro-structural properties of brain tissue 
such as axon, myelin, iron and water concentration (Weiskopf et al., 2015).

The hMRI-toolbox is an easy-to-use open-source and flexible tool, for qMRI data handling and processing. 
It allows the estimation of high-quality multi-parameter qMRI maps 
(longitudinal and effective transverse relaxation rates R1 and R2*, proton density PD 
and magnetisation transfer MT saturation) (Weiskopf et al., 2013), 
followed by spatial registration in common space for statistical analysis (Draganski et al., 2011).

Embedded in the [Statistical Parametric Mapping](http://www.fil.ion.ucl.ac.uk/spm) (SPM) framework, 
it can be readily combined with existing SPM toolboxes for estimating diffusion MRI parameter maps, 
and it benefits from the extensive range of established SPM tools for high-accuracy spatial registration and statistical inferences.

The qMRI maps generated by the toolbox can be used for quantitative parameter analysis and accurate delineation of subcortical brain structures. They are key input parameters for biophysical models designed to estimate tissue microstructure properties such as the MR g-ratio and to derive standard and novel MRI biomarkers (Mohammadi et al., 2015). The hMRI toolbox is therefore a first step towards *in vivo* histology using MRI (hMRI) and is being extended further in this direction.

## Wiki-Pages

Online documentation is available as a [Wiki](https://github.com/hMRI-group/Toolbox/wiki). 
It provides guidelines and instructions for installation and usage of the hMRI-toolbox. 
These pages are work-in-progress and updated on a regular basis.

## hMRI-Toolbox and example dataset papers

For a reference on the scientific background, methods and concepts please use this [paper](https://doi.org/10.1016/j.neuroimage.2019.01.029) (formerly available as a [pre-print paper](http://dx.doi.org/10.20347/WIAS.PREPRINT.2527)) and cite it when publishing results compiled with the hMRI-toolbox. 

A full example dataset can be `git clone`ed or downloaded from one of the following links: ([Github](https://github.com/hMRI-group/hMRI-toolbox-demodata) | [OwnCloud](https://owncloud.gwdg.de/index.php/s/znxsMrPUB2TJsLz) | [Mega](https://mega.nz/folder/OnxDALZa#bFe0xKJTuuUeGXrYnlX8dA)). The description of the example dataset is also available in this [paper](https://doi.org/10.1016/j.dib.2019.104132).

## Licence

The hMRI toolbox is free but copyright software, distributed under the terms of the GNU General Public Licence as published by the Free Software Foundation (either version 2, as given in file LICENSE, or at your option, any later version). Further details on "copyleft" can be found at http://www.gnu.org/copyleft/. In particular, the hMRI toolbox is supplied as is. No formal support or maintenance is provided or implied.

## Download

The latest release (and previous releases and pre-releases) of the hMRI-toolbox Matlab code 
can be downloaded as a zip archive from the [releases page](https://github.com/hMRI-group/hMRI-toolbox/releases). 

## E-Mail List

We have created an e-mail list for users of the hMRI-toolbox: HMRI-TOOLBOX@JISCMAIL.AC.UK. 
[Registered users](https://www.jiscmail.ac.uk/cgi-bin/webadmin?SUBED1=HMRI-TOOLBOX&amp;A=1) can login 
to view the message archive on the [list Home Page](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=HMRI-TOOLBOX).

## Developers of the hMRI-toolbox

The development of the hMRI toolbox is an international collaborative effort including the following sites and developers: 

- Luke J Edwards, Tobias Leutritz, Enrico Reimer, Baris Evren Ugurcan, Nikolaus Weiskopf (Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig, Germany) 
- Evelyne Balteau, Christophe Phillips ([University of Liege](https://www.giga.uliege.be), Liege, Belgium) 
- Siawoosh Mohammadi (Medical Center Hamburg-Eppendorf, Hamburg, Germany) 
- Martina F Callaghan, John Ashburner (University College London, London, United Kingdom)
- Karsten Tabelow (Weierstrass Institute for Applied Analysis and Stochastics, Berlin, Germany)
- Bogdan Draganski, Ferath Kerif, Antoine Lutti  (LREN, DNC - CHUV, University Lausanne, Lausanne, Switzerland)
- Maryam Seif (University of Zurich, Zurich, Switzerland) 
- Gunther Helms (Department of Medical Radiation Physics, Lund University, Lund, Sweden)
- Lars Ruthotto (Emory University, Atlanta, GA, United States) 
- Gabriel Ziegler (Otto-von-Guericke-University Magdeburg, Magdeburg, Germany)

## Toolbox Reference

Please cite this key reference when you use the toolbox:
- Tabelow, K., Balteau, E., Ashburner, J., Callaghan, M. F., Draganski, B., Helms, G., Kherif, F., Leutritz, T., Lutti, A., Phillips, C., Reimer, E., Ruthotto, L., Seif, M., Weiskopf, N., Ziegler, G., Mohammadi, S., 2019. hMRI – A toolbox for quantitative MRI in neuroscience and clinical research. Neuroimage 194, 191-210. (https://doi.org/10.1016/j.neuroimage.2019.01.029).

## Additional References:
Please also cite the relevant references for the methods that are implemented in the toolbox detailed below:

### Quantitative map creation
- Helms, Gunther, Henning Dathe, and Peter Dechent. 2008. “Quantitative FLASH MRI at 3T Using a Rational Approximation of the Ernst Equation.” Magnetic Resonance in Medicine 59(3):667–72. (http://www.ncbi.nlm.nih.gov/pubmed/18306368).
- Helms, Gunther, Henning Dathe, Kai Kallenberg, and Peter Dechent. 2008. “High-Resolution Maps of Magnetization Transfer with Inherent Correction for RF Inhomogeneity and T1 Relaxation Obtained from 3D FLASH MRI.” Magnetic Resonance in Medicine 60(6):1396–1407. (http://www.ncbi.nlm.nih.gov/pubmed/19025906).
- Weiskopf, Nikolaus, Martina F. Callaghan, Oliver Josephs, Antoine Lutti, and Siawoosh Mohammadi. 2014. “Estimating the Apparent Transverse Relaxation Time (R2*) from Images with Different Contrasts (ESTATICS) Reduces Motion Artifacts.” Frontiers in Neuroscience 8(September):1–10. (http://www.frontiersin.org/Brain_Imaging_Methods/10.3389/fnins.2014.00278/abstract).
- Weiskopf, N., Suckling, J., Williams, G., Correia, M.M., Inkster, B., Tait, R., Ooi, C., Bullmore, E.T., Lutti, A., 2013. Quantitative multi-parameter mapping of R1, PD*, MT, and R2* at 3T: a multi-center validation. Front. Neurosci. 7, 95. (https://doi.org/10.3389/fnins.2013.00095).
- Leutritz, Tobias et al. 2020. “Multiparameter Mapping of Relaxation (R1, R2*), Proton Density and Magnetization Transfer Saturation at 3 T: A Multicenter Dual-Vendor Reproducibility and Repeatability Study.” Human Brain Mapping 41(15):4232–47. (https://pubmed.ncbi.nlm.nih.gov/32639104/).



### Artefact correction
- Lutti, Antoine et al. 2012. “Robust and Fast Whole Brain Mapping of the RF Transmit Field B1 at 7T.” PloS one 7(3):e32379. (http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3299646).
- Lutti, Antoine, Chloe Hutton, Jürgen Finsterbusch, Gunther Helms, and Nikolaus Weiskopf. 2010. “Optimization and Validation of Methods for Mapping of the Radiofrequency Transmit Field at 3T.” Magnetic Resonance in Medicine 64(1):229–38. (http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3077518).
- Weiskopf, Nikolaus et al. 2011. “Unified Segmentation Based Correction of R1 Brain Maps for RF Transmit Field Inhomogeneities (UNICORT).” NeuroImage 54(3):2116–24. (http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3018573).
- Papp, Daniel, Martina F. Callaghan, Heiko Meyer, Craig Buckley, and Nikolaus Weiskopf. 2016. “Correction of Inter-Scan Motion Artifacts in Quantitative R1 Mapping by Accounting for Receive Coil Sensitivity Effects.” Magnetic Resonance in Medicine 76(5):1478–85. (http://doi.wiley.com/10.1002/mrm.26058).
- Corbin, Nadège and Martina F. Callaghan. 2021. “Imperfect Spoiling in Variable Flip Angle T 1 Mapping at 7T: Quantifying and Minimizing Impact.” Magnetic Resonance in Medicine 86(2):693–708. (https://onlinelibrary.wiley.com/doi/10.1002/mrm.28720).



### Group analysis / statistical methods
- Draganski, B., Ashburner, J., Hutton, C., Kherif, F., Frackowiak, R.S.J., Helms, G., Weiskopf, N., 2011. Regional specificity of MRI contrast parameter changes in normal ageing revealed by voxel-based quantification (VBQ). Neuroimage 55, 1423-1434. (https://doi.org/10.1016/j.neuroimage.2011.01.052).


### Applications
- Mohammadi, S., Carey, D., Dick, F., Diedrichsen, J., Sereno, M.I., Reisert, M., Callaghan, M.F., Weiskopf, N., 2015. Whole-Brain In-vivo Measurements of the Axonal G-Ratio in a Group of 37 Healthy Volunteers. Front Neurosci 9, 441. (https://doi.org/10.3389/fnins.2015.00441).
- Callaghan, Martina F. et al. 2014. “Widespread Age-Related Differences in the Human Brain Microstructure Revealed by Quantitative Magnetic Resonance Imaging.” Neurobiology of Aging 35:1862–72. (http://linkinghub.elsevier.com/retrieve/pii/S0197458014002000).

### Reviews
- Weiskopf, N., Mohammadi, S., Lutti, A., Callaghan, M.F., 2015. Advances in MRI-based computational neuroanatomy: from morphometry to in-vivo histology. Curr. Opin. Neurol. 28, 313-322. (https://doi.org/10.1097/WCO.0000000000000222).
- N. Weiskopf, L. Edwards, G. Helms, S. Mohammadi, and E. Kirilina. 2021. “Quantitative Magnetic Resonance Imaging of Brain Anatomy: Towards in-Vivo Histology.” Nature Reviews Physics. (http://dx.doi.org/10.1038/s42254-021-00326-1).


### Datasets
- Callaghan, M. F., Lutti, A., Ashburner, J., Balteau, E., Corbin, N., Draganski, B., Helms, G., Kherif, F., Leutritz, T., Mohammadi, S., Phillips, C., Reimer, E., Ruthotto, L., Seif, M., Tabelow, K., Ziegler, G., Weiskopf, N., 2019. Example dataset for the hMRI toolbox. Data in Brief 25, 104132. (https://doi.org/10.1016/j.dib.2019.104132).





## Acknowledgments and Funding

- EB received funding from the European Structural and Investment Fund / European Regional Development Fund & the Belgian Walloon Government, project BIOMED-HUB (programme 2014-2020). 
- NW received funding from the European Research Council under the European Union's Seventh Framework Programme (FP7/2007-2013) / ERC grant agreement No 616905. This project has received funding from the European Union's Horizon 2020 research and innovation programme under the grant agreement No 681094, and is supported by the Swiss State Secretariat for Education, Research and Innovation (SERI) under contract number 15.0137. 
- SM received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 658589. 
- NW and SM received funding from the BMBF (01EW1711A and B) in the framework of ERA-NET NEURON. BD is supported by the Swiss National Science Foundation (NCCR Synapsy, project grant Nr 32003B_159780) and Foundation Synapsis. LREN is very grateful to the Roger De Spoelberch and Partridge Foundations for their generous financial support. 
- MFC is supported by the MRC and Spinal Research Charity through the ERA-NET Neuron joint call (MR/R000050/1). 
- The Wellcome Centre for Human Neuroimaging is supported by core funding from the Wellcome [203147/Z/16/Z]. 
- CP is supported by the [F.R.S.-FNRS](http://www.fnrs.be/en/), Belgium. 
- The hMRI Toolbox project is supported by the Max Planck Society.
