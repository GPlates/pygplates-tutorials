FROM ubuntu:18.04
MAINTAINER Michael Chin

#Install ubuntu libraires and packages
RUN apt-get update -y
RUN apt-get install -y wget python3-pip

#Download pygplates and install it
ADD pygplates_2.0_amd64.deb .
RUN apt-get install -y ./pygplates_2.0_amd64.deb
RUN rm pygplates_2.0_amd64.deb

Env PYTHONPATH ${PYTHONPATH}:/usr/lib:/usr/lib/pygplates/revision26/

# install all the python and ipython notebook requirements
RUN pip3 install --upgrade pip
RUN pip3 install numpy scipy matplotlib jupyter pandas sympy nose
RUN pip3 install ipyparallel pyproj pyshp Pillow
RUN pip3 install moviepy

RUN apt-get install libgeos-dev
RUN wget https://github.com/matplotlib/basemap/archive/v1.2.1rel.tar.gz
RUN tar -vxf v1.2.1rel.tar.gz
RUN cd basemap-1.2.1rel/ && python3 setup.py install_lib
RUN rm -rf basemap-1.2.1rel/ v1.2.1rel.tar.gz

RUN pip3 install Cython
RUN apt-get install -y libproj-dev
RUN wget https://github.com/SciTools/cartopy/archive/v0.18.0b1.tar.gz
RUN tar -vxf v0.18.0b1.tar.gz
RUN cd cartopy-0.18.0b1/ && python3 setup.py install
RUN rm -rf cartopy-0.18.0b1/ v0.18.0b1.tar.gz

RUN pip3 install statistics netCDF4

RUN pip3 uninstall -y shapely
RUN pip3 install shapely --no-binary shapely

ARG CACHE_DATE=2019-02-26

RUN mkdir /home/workspace
WORKDIR /home/workspace

# Add Tini. Tini operates as a process subreaper for jupyter. This prevents
# kernel crashes.
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/v0.18.0/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]

RUN useradd -m -s /bin/bash -N mchin

CMD ["/usr/local/bin/jupyter", "notebook", "--no-browser", "--ip=0.0.0.0", "--NotebookApp.token=''"]

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

EXPOSE 8888

RUN chown mchin /tmp

USER mchin
