# pygplates-tutorials
## Introduction
This repository contains jupyter notebooks which demonstrate how to use pygplates in plate tectonic research.The notebooks are designed to be a companion to the pygplates documentation. The documentation contains sample code illustrating in a general sense how to perform certain tasks; these notebooks contain worked examples on actual data sets.

## Getting Started
### Step 1: Download the repository 

```git clone --depth=1 https://github.com/GPlates/pygplates-tutorials.git```
### Step 2: Pull the gplates/pygplates-notebook docker container 
(install Docker https://www.docker.com/ if you have not done so.)

```docker pull gplates/pygplates-notebook```

## System Requirements:
The fundamental requirement to run these notebooks is to have pygplates installed (compatible with Python 2.7, but not Python 3.X). Installation instructions can be found in the pygplates user documentation.
Other python modules used by the notebooks are:
- numpy
- matplotlib (including the basemap extension)
- pandas

## Docker:
##### build the docker image:
cd docker

docker build -t gplates/pygplates-notebook .

##### run the docker container:
cd pygplates-tutorials

docker run -it --rm -p 18888:8888 -v \`pwd\`:/home/workspace gplates/pygplates-notebook

## Contact and more info:
For general information on GPlates, please see the gplates website:
www.gplates.org

The pygplates documentation can be found here:
www.gplates.org/docs/pygplates/index.html

If you have issues or questions, please consider directing them to the GPlates mailing list:

###### gplates-discuss(at)mailman(dot)sydney(dot)edu(dot)au  

## About 
These tutorials stem from the work of Simon Williams, Michael Tetley, John Cannon and Michael Chin at
EarthByte Group, University of Sydney, 2016-present
