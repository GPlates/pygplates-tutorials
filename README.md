# pygplates-tutorials
Repository for jupyter (formerly ipython) notebooks demonstrating functionality of pygplates

The notebooks are designed to be a companion to the pygplates documentation. The documentation contains sample 
code illustrating in a general sense how to perform certain tasks; these notebooks contain worked examples on actual
data sets.

For general information on GPlates, please see the gplates website:
www.gplates.org

The pygplates documentation can be found here:
www.gplates.org/docs/pygplates/index.html

Data files used within the tutorials is contained within the 'Data' folder. These plate model files (*.gpml and *.rot)
are taken from the sample code released with version 1.5 of the GPlates desktop application.
When running the notebooks, various files will be produced which by default will be created within the /tmp folder.

#### System Requirements:
The fundamental requirement to run these notebooks is to have pygplates installed (compatible with Python 2.7, but not Python 3.X). Installation instructions can be found in the pygplates user documentation.
Other python modules used by the notebooks are:
- numpy
- matplotlib (including the basemap extension)
- pandas

#### Docker:
##### build the docker image:
cd docker

docker build -t gplates/pygplates-notebook .

##### run the docker container:
cd pygplates-tutorials

docker run -it --rm -p 18888:8888 -v \`pwd\`:/home/workspace gplates/pygplates-notebook

#### Contact:
If you have issues or questions, please consider directing them to the GPlates mailing list:

###### gplates-discuss(at)mailman(dot)sydney(dot)edu(dot)au  

#### About 
These tutorials stem from the work of Simon Williams, Michael Tetley, John Cannon, Michael Chin 
EarthByte Group, University of Sydney, 2016-present
