**[Technical Overview](#technical-overview)** |
**[Installation](#installation)** |
**[Contributing](#contributing)** |
**[License](#license)** |
**[Help and Resources](#help-and-resources)** 

# FTIR Data Analytics

With the FTIR data analysis package you can read FTIR data, perform automated
buffer subtractions, spectral offsetting corrections, and deconvolute 
second derivative data.

This project was started as a way to improve spectra deconvolution of FTIR 
data, and only the tools necessary to analyze FTIR absorbance data have been 
built by KBI Biopharma. 

## Technical Overview

## Installation
### Check prerequisites
- [Python](https://www.python.org/downloads/) 3.5 or greater
- Package dependencies are listed below. Validation of package versions has not
been performed, however, only basic functionality from each package is used, 
and it is expected that ftir_data_analytics will work with almost any version 
of these packages. 
    * pandas 
    * numpy
    * matplotlib
    * scipy

### Install packages
Currently, only the source code for this package is provided. The package can 
be installed using python `setuptools` in a virtual environment on your system.
The source code must first be down loaded onto the local users machine. Using
the python environment of your choice, install this FTIR package by running the
following command using bash, the command prompt, or alternative in the root 
folder of the code directory (will contain the file `setup.py`)

```bash
python setup.py install
```

Creating python virtual environments for data analysis can be done using a 
variety of tools. [Conda](https://conda.io/docs/) and the [Enthought Tool 
Suite](http://code.enthought.com/) are the preferred virtual environment 
management tools at KBI Biopharma.

### Check Installation
The validity of the install can be confirmed by importing the FTIR data 
analysis tools into your python environment using `import ftir`. A more robust
test suite may be developed depending on the utilization of these tools. 


## Contributing
If you would like to contribute to the project, please contact Brent Kendrick
or Jared Young at KBI Biopharma. 


## License
We use a shared copyright model that enables all contributors to maintain the
copyright on their contributions.

All code is licensed under the terms of the revised BSD license.


## Help and Resources
- [Reporting Issues](https://github.com/liljerri/ftir_data_analytics/issues)


---

**[Technical Overview](#technical-overview)** |
**[Installation](#installation)** |
**[Contributing](#contributing)** |
**[License](#license)** |
**[Help and Resources](#help-and-resources)** 