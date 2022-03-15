Setup and Installation
==========================================================================================

Development and documentation occurs on GitHub_.

PIDGIN is currently only compatible with Python 2.7.x. but it runs under Anaconda 3 installations as long as the environment is set up with Python 2.7 (see  installation instructions below).

It also has the following dependencies:

Required dependencies
~~~~~~~~~~~~~~~~~~~~~

- NumPy_
- SciPy_
- RDKit_
- Scikit-learn_
- Standardiser_
- python_utilities_


Install with Conda
~~~~~~~~~~~~~~~~~~

Follow these steps on Linux/OSX:

1. Download and install Anaconda from https://www.continuum.io/downloads

2. Open terminal in Mac/Linux and run ``conda env create -f pidgin4_env.yml --name pidgin4_env``

* N.B. Rdkit may not import on some systems due to a bug. If this happens upgrade to the latest version of conda before creating the above environment using: ``conda update conda``

* N.B. Installs the IMI eTOX `flatkinson standardiser`_ (replaces ChemAxon's standardizer used in previous PIDGIN versions) and statsmodels for p-value correction in predict_enriched.py

If you encounter an issue (usually occurs when installing the environment on non-Linux systems) try the following:

``conda create -c rdkit -c conda-forge --name pidgin4_env python=2.7 rdkit scikit-learn=0.19.0 pydot graphviz standardiser statsmodels``

3. Now run: ``source activate pidgin4_env`` (This activates the PIDGINv4 virtual environment. N.B This is required for each new terminal session in order to run PIDGIN in the future)

4. Navigate the directory you wish to install PIDGINv4 and in Mac/Linux terminal run ``git clone https://github.com/BenderGroup/PIDGINv4.git`` (recommended) or download/extract the zip from `GitHub`_ webpage (not recommended due to inability to pull updates)

5. Download and unzip the no_ortho_mar22.tar.gz https://doi.org/10.6084/m9.figshare.19108382.v1 (md5sum: dc146e69c8f1638e3741ff7900a97cf3) into the PIDGINv4 main directory to form the no_ortho/ directory (leave all subsequent files compressed)

* N.B Depending on bandwidth, Step 5 may take some time

NOTE: For older models and orthologue models, please contact us.

.. [1] |mervin2018|

.. include:: substitutions.rst
