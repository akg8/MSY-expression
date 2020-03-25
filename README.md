# MSY-expression

Figures and tables can be reconstructed by running the [Jupyter](https://jupyter-notebook.readthedocs.io/en/stable/) notebook files (`*.ipynb`). These notebook files draw upon functions and variables saved in the Python (v3.6.9) modules in the `msyexp` subdirectory, in addition to the third-party packages named below.

Dependencies: Analyses were conducted in a series of Jupyter notebooks written in Python, using Python packages [numpy](https://numpy.org) (v1.17.2), [scipy](https://scipy.org) (v1.3.1), [pandas](http://pandas.pydata.org) (v0.25.1), [pytables](https://www.pytables.org) (v3.5.2), [scikit-learn](https://scikit-learn.org/stable/) (v0.21.3), [statsmodels](https://www.statsmodels.org/stable/index.html) (v0.10.1), [matplotlib](https://matplotlib.org) (v3.1.1), and [seaborn](https://seaborn.pydata.org) (v0.9.0). All were run on a 2017 MacBook Pro with 16GB of RAM. Python packages can be downloaded and installed individually, or the virtual environment in which these analyses were run can be reconstructed from the environment file `msy3.environment.yml`, e.g.:

    conda env create -f msy3.environment.yml

To run the analyses, clone this repository to a local directory. Then download all files from the associated data repository on Zenodo (DOI: [10.5281/zenodo.3627233](https://doi.org/10.5281/zenodo.3627233), file: `msy_expression.all_input_data.zip`) to a local directory (e.g., `/Users/janedoe/msy_input_files`). 

To make sure the notebook files can find the downloaded data, edit the `DATADIR` variable in the `msyexp.paths.py` module and point it to the relevant local directory. Set the `NBOUTDIR` variable to the directory where you would like all output files to be saved.

    ##### ~~~ EDIT HERE ~~~ #######################################################
    
    # local directory containing all input files downloaded from Zenodo
    # (https://doi.org/10.5281/zenodo.3627110)
    DATADIR = '/Users/janedoe/msy_input_files'
    
    # local directory where you would like all output files to be written;
    # each notebook's files will be saved within a notebook-specific subdirectory
    NBOUTDIR = '/Users/janedoe/msy_output_files'
    
    ##### ~~~ EDIT HERE ~~~ #######################################################
    
Now, run the individual Jupyter notebook files. The approximate order for running these notebooks is given below, however, all but sexbias_analyses.ipynb can be run independently of the others.

1. mappability_analyses.ipynb
2. MSY_overview_analyses.ipynb
3. simulation_analyses.ipynb
4. yxratio_analyses.ipynb
5. coexpression_analyses.ipynb
6. diffexp_analyses.ipynb
7. eif1a_analyses.ipynb
8. eif1a_HPA_replication_analyses.ipynb
9. sexbias_analyses.ipynb
10. massspec_analyses.ipynb

