# Installation

Following are platform dependent installation instructions


## Linux/Mac/Windows (x64)
1. Create a new conda environment.
2. Clone the repository.
3. `conda install --file requirements. txt`

## Jetson Nano (aarch64)

Following installation was tested on `4.9.140-tegra aarch64 GNU/Linux`.

We use a python package manager to isolate the dependencies required for pysero.
conda package manager is not directly supported by above version of Linux.
Instead, pip virtual environment works well for package management.

Major steps in setup are:
1.	Install virtual environment: 
    `sudo apt install -y python3-venv`
    
2.	Create an environment: 
    `python3 -m venv ~/python-envs/pysero`
    
3.	Activate environment: 
    `source ~/python-envs/pysero/bin/activate`
    
4.	Deactivate environment: 
    `deactivate`
    
5. Update pyenv.cfg to use system-wide packages:
    `cd ~/python-envs/env/pysero`
    
    Use `nano pyvenv.cfg to` change:
    `include-system-site-packages = true` 
    (allow packages that require lower-level installation to be accessed in the environment)

Some packages need to be installed in local environment and some system-wide.
Use `pip install -I` to install in the `pysero` site-packages folder, the environment will access the local pacakages first and then search for global packages.

### install dependencies for python packages:

scikit-image dependencies:

        sudo apt-get install python-dev libfreetype6-dev
	sudo apt-get install libfontconfig1-dev
	
	
scipy dependencies:

	sudo apt-get install libblas-dev liblapack-dev libatlas-base-dev gfortran
	sudo -H python3 -m pip install scipy==1.1.0
	
### install python packages:
	
	
	pip install -I scikit-image 
	# compilation of scikit-image takes >5min. May have to try sudo -H pip3 install scikit-image
	pip install -I scipy==1.1.0
	# compilation of scipy can take >5 min.
	pip install -I matplotlib 
	pip install -I pandas==0.24.0 
	# compilation of pandas can take time.
	pip install -I wheel
	pip install -I certifi
	pip install -I cycler
	pip install -I kiwisolver
	pip install -I cython
	pip install -I numpy==1.18.1
	pip install -I python-dateutil
	pip install -I pytz
	pip install -I six
	pip install -I tabulate
	pip install -I tornado
	pip install -I openpyxl
	pip install -I xmltodict
	pip install -I xlrd





