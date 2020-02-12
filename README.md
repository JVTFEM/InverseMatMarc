# InverseMatMarc
Inverse material optimisation using Marc from two FE models with remeshing
# About
This respitory contains all the steps, source code and a complete example problem to perform inverse material characterisation from two FE models, with the idea that one can be replaced with digital image correlation (DIC) data. The process make use of MSC Marc to perform the FE analysis, with remeshing as an interesting feature. DOT (Design Optimization Tools from [VR&D](http://www.vrand.com/ "VR&D Home") is the optimiser used within the source code, which provides a library for gradient-based optimisation algorithms.
## Information
This project was setup to work with Windows. The DOT optimisation library is required, with a few environmental variable changes when using the Python wrapper for DOT. If you are not familiar with DOT, the following respitrory from the MODRG could be of help, [DOTWrapper](https://github.com/MODRG/DOTWrapper).
## Requirements
* MSC Marc Mentat to perform the FE analysis. Marc is a commercial software and will need a license to run.
* [PyMentat](https://simcompanion.mscsoftware.com/infocenter/index?page=content&id=DOC11109) and [PyPost](https://simcompanion.mscsoftware.com/infocenter/index?page=content&channel=DOCUMENTATION&cat=MARC_DOCUMENTATION&sort=displaystartdate&dir=descending&max=1000&batch=15&rss=true&itData.offset=0) Python API's to interface with Marc IO files. They are automatically installed with Marc.
* A [Python](https://www.python.org/downloads/) .exe of the same version as the Python API's.
* [Pandas](https://pandas.pydata.org/) Python analysis library
* [NumPy](https://numpy.org/) Python scientific computing library
* [scipy.interpolate](https://docs.scipy.org/doc/scipy/reference/interpolate.html) Sub-package for objects used in interpolation
* [MatPlotlib](https://matplotlib.org/) Python 2D plotting library
* [pyDOE](https://pythonhosted.org/pyDOE/) Python experimental design package
## File description
| File         | Description         |
|:-------------|:--------------------|
|functions.py  | The file that contains all the written functions used to create files and save data for later access.|
|dot11.py      | The Python wrapper for DOT, but edited for use within the pipeline.|
|File_paths.py | The file that contains all the file paths needed within the pipeline to access and save the different files and data created. |
|graphs.py     | The final example post file used to create graphs from the results obtained after optimisation. Optional file. |
|LHC.py        | The file containing the Latent Hypercube (LHC) function specifically for this pipeline. |
|main_code.py  | The file that contains the main code from where all of the other files are called and the whole pipeline is controlled from. This is the main file to run within Python. |
|RBF.py        | The file that does the necessary post processing of interpolating the data such that the objective function can be created.|
|RMS.py        | The file that calculates the Root Mean Square (RMS) which is used as the objective function. |
|RSQ.py        | The file that calculates the R^2 error between the data sets for the graphs.py file. Optional file, but required if graphs.py is used. |
|SEOE.py       | The file that calculates the Standard Error of Estimation for the results. Optional file, but is required if graphs.py is used. |
|sym_test_fem_job1.dat | The input file from the Numerical model for the example. |
|sym_test_exp_job1.dat | The input file from the Experimental model for the example. |
|sym_test_fem_job1.t16 | The output file from the Numerical model for the example. |
|sym_test_exp_job1.t16 | The output file from the Experimental model for the example. |
|sym_test_fem_job1.sts | The Status file from the Numerical model. Used to determine the EXIT CODE. |
|sym_test_fem.mud      | The FE model for the Numerical model for the example. |
|sym_test_exp.mud      | The FE model for the Experimental model for the example. |

## Example problem
This Numerical pipeline was created as part of a Master's thesis. The aim of the project was to investigate if micro-indentation can characterise soft materials. Two FE models were used for the investigation. It was not known if the indentation method would work and hence the use of two FE models instead of a physical experiment and one FE model. Both FE models were thus modelled using the same testing conditions, but with minor changes explained further. The main goal behind the two FE models were that one model represents DIC data, the "Experimental" model, which used a constant Mooney-Rivlin three parameter material model. The second FE model, the "Numerical" model, is the model used and within the inverse FE analysis and the material model here was changed and optimised during the optimisation procedure to determine the feasibility of the indentation method by comparing the displacement results. 

Both FE models used a feature called [remeshing](https://www.mscsoftware.com/product/marc), which in short involves that the mesh changes during every remeshing instance, by replacing it with a better quality mesh with smaller/ larger elements depending on the user criteria and the sample deformation. This causes the mesh to gain more/less nodes, therefore the number of nodes changes during the solving increments and therefore their id's as well. A node tracking method was needed since the built in node tracker only worked on tet4 elements and not tet10 elements. With the aim that the experimental model should represent DIC data, only surface displacements were used and therefore, only one surface on the side of the sample which is identified as the data capturing surface. The node tracking problem was overcome by glueing another surface, meshed with membrane elements with a 0.001 mm thickness as not to change the sample's bending stiffness, to the sample's data capturing surface. During the analysis, these membrane elements do not undergo remeshing and these nodes could be tracked for data capturing. This extra layer was called the "skin_elements" and "skin_nodes", which were the key phrases used to identify the nodes for tracking in the input file. 
Convergence studies showed that the membrane elements need to be meshed at least twice as small as the sample to obtain accurate results. The "Experimental" model was also meshed with a smaller element size as to resemble DIC data through having more data points at different locations, to that of the "Numerical" model.

The example problem considered here, used a 3 mm spherical indenter to apply an indentation, to sample of size 20mm x 20mm x 5mm. A linear position ramp function of 3mm was applied to cause a 3mm indentation in 1s, making the indentation depth a function of time. The interpolation steps uses the indenter's y-displacement. The example problem has to obtain the Mooney-Rivlin three parameter material properties (C10, C01 and C20) for a square sample shown below.

![alt text][example]

[example](https://github.com/Franciena/InverseMatMarc/blob/master/pictures/sph_mid.png)

The exact material properties for the "Experimental" model was thus known as:
| Parameter | Value |
|:----------|:------|
| C10       |0.2605676 MPa |
| C01       |0.0975498 MPa |
| C20       |0.0575007 MPa |

The RMS calculation accounts for the bias due to the different orders in magnitude between the parameters.

# Getting started
Note that [Visual Studio Code](https://code.visualstudio.com/download) was used as the GUI for Python. A quick step to step guide will be given how each software was installed and setup to work with each other. MSC Marc 2019 was used and according to the PyMentat and PyPost API's associated, Python 3.6.4 was used.

## MSC Marc Mentat
1. Install Marc
2. Go to the Program Files and check which Python version works with Marc version, ex.:
```
C:\MSC.Software\Marc\2019.0.0\mentat2019\python\WIN8664\python36.dll
```
3. Download a Python 3.6 version
4. Install Python 3.6
5. From folder:
```
C:\MSC.Software\Marc\2019.0.0\mentat2019\shlib\Win64
```
Copy files:
```
py_mentat.pyd
py_post.pyd
```
And paste them in the Python folder:
```
C:\Program Files (x86)\Python36\DLLs
```
## Visual Studio Code and Python
1. Run the cmd prompt as administrator
2. Upgrade pip in Python:
```
> python -m pip install --upgrade pip setuptools wheel
```
3. Install Python libraries:
```
> python -m pip install --user numpy scipy matplotlib pandas sympy pyDOE
```
1. Install VS Code
2. Open VS Code
3. Install Python extension:
```
File > Preferences > Extensions > Install Python extension
```
4. Select Python interpreter:
```
Ctrl + Shift + P: Select the Python 3.6 version
```
## DOT wrapper
DOT is a commercial software and need a license to work. DOT can be installed through VisualDOC, which is obtainable from VR&D. During the VisualDOC installation, it is optional to install DOT with the whole VisualDOC software package or just DOT on its own. Since VisualDOC was not needed for this pipeline, only DOT was installed from the VisualDOC installation. After DOT is installed the following environmental variables need changing.
1. Create a new environment variable, VRAND_AUT and point it to the license file:
```
VRAND_AUT = C:\...\vrand\licenses\vrand.lic
```
2. Copy the DOT shared library (.dll) in the Python DLLs forder:
```
C:\Program Files (x86)\Python36\DLLs
```
3. Add the path to your DOT shared library (.dll) to your current PATH environment variable:
```
PATH = C:\Program Files (x86)\Python36\DLLs; C:\...\vrand\dot6.0\Win64
```
4. Make sure the dot11.py file is in the same folder as the main_code.py
5. The Python file main_code.py should be ready to run. The [DOTWrapper](https://github.com/MODRG/DOTWrapper) respritory provides an unedited version of the dot11.py file with a test example to test if DOT is working.
## Running the pipeline
The best method to ensure the pipeline is running smoothly, is to download the InverseMatMarc.zip file and running the main_code.py directly from there. But if each code is downloaded separately, the files need to be arranged as follow:
1. All the .py files should be in the same folder, example:
```
C:\...\InverseMatMarc
```
2. All the IO files should be in the same sub-folder named "sph_mid" or any given folder name as long as the "File_paths.py" file is adjusted accordingly:
```
C:\...\InverseMatMarc\sph_mid
```
3. Run the "main_code.py" file in the InverseMatMarc folder.
