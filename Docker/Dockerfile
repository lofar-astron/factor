FROM ubuntu:18.04
LABEL Maintainer="iacobelli@astron.nl;loose@astron.nl"
ENV DEBIAN_FRONTEND=noninteractive
ENV HDF5_USE_FILE_LOCKING=FALSE

# install common dependencies
RUN apt-get update && apt-get --yes install --no-install-recommends bison build-essential cmake emacs flex g++ gcc gettext-base gfortran git libarmadillo-dev libblas-dev libboost-date-time-dev libboost-dev libboost-filesystem-dev libboost-numpy-dev libboost-program-options-dev libboost-python-dev libboost-regex-dev libboost-signals-dev libboost-system-dev libboost-thread-dev libcfitsio-dev libfftw3-dev libgsl-dev libgtkmm-3.0-dev libhdf5-serial-dev liblapacke-dev liblog4cplus-1.1-9 liblog4cplus-dev libncurses5-dev libpng-dev libpython2.7-dev libreadline-dev libxml2-dev openssh-server python python-pip python-tk python-setuptools subversion vim wcslib-dev wget && rm -rf /var/lib/apt/lists/*

# install python2 packages
RUN pip install setuptools wheel Cython
RUN pip install --upgrade aplpy astropy Jinja2 matplotlib numpy==1.16 PySocks python-monetdb scipy shapely wcsaxes xmlrunner

#
RUN mkdir /opt/lofarsoft/

# install boost python 1.63 with python2
RUN cd / && wget https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.bz2/download -O boost_1_63_0.tar.bz2 && tar xvf boost_1_63_0.tar.bz2 && cd boost_1_63_0 && ./bootstrap.sh --with-python=/usr/bin/python --with-libraries=python,date_time,filesystem,system,program_options,test && ./b2 install && cd / && rm -r boost_1_63_0*

# install casacore data
RUN cd / && mkdir -p /opt/lofarsoft/data && cd /opt/lofarsoft/data && wget ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar && tar xvf WSRT_Measures.ztar && rm WSRT_Measures.ztar

# install casacore v2.4.1
RUN cd / && wget https://github.com/casacore/casacore/archive/v2.4.1.tar.gz && tar xvf v2.4.1.tar.gz && cd casacore-2.4.1 && mkdir build && cd build && cmake -DPORTABLE=True -DBUILD_PYTHON=True -DDATA_DIR=/opt/lofarsoft/data -DUSE_OPENMP=ON -DUSE_THREADS=OFF -DUSE_FFTW3=TRUE -DUSE_HDF5=ON -DCXX11=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft/ -DCMAKE_CXX_FLAGS="-fsigned-char -O2 -DNDEBUG -march=native" ../ && make -j4 && make install && cd ../../ && rm -rf casacore-2.4.1 v2.4.1.tar.gz

# install python casacore v3.2.0
RUN cd / && wget https://github.com/casacore/python-casacore/archive/v3.2.0.tar.gz && tar xvf v3.2.0.tar.gz && cd python-casacore-3.2.0 && python setup.py build_ext -I/opt/lofarsoft/include -L/opt/lofarsoft/lib/ && mkdir -p /opt/lofarsoft/lib/python2.7/site-packages/ && export PYTHONPATH=/opt/lofarsoft/lib/python2.7/site-packages/ && python setup.py install --prefix=/opt/lofarsoft && cd ../ && rm -rf v3.2.0.tar.gz

# install AOFlagger 2.14 last py2 version compatible
RUN cd / && wget https://sourceforge.net/projects/aoflagger/files/aoflagger-2.14.0/aoflagger-2.14.0.tar.bz2/download && mv download download.tar && tar xvf download.tar && cd aoflagger-2.14.0 && mkdir build && cd build && cmake ../ -DPORTABLE=True -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft -DCMAKE_PREFIX_PATH=/opt/lofarsoft && make -j4 && make install && cd ../../ && rm -rf download.tar aoflagger-2.14.0

# install losoto 2.0
RUN cd / && git clone https://github.com/revoltek/losoto.git && cd losoto && git checkout c8fbd6194074bef4009cb66dd7ecd59e98664d63 && export PYTHONPATH=/opt/lofarsoft/lib/python2.7/site-packages/ && pip install numexpr==2.7.3 tables==3.5.2 configparser==4.0.2 && python setup.py install --prefix=/opt/lofarsoft && cd / && rm -rf losoto

# install RMextract v0.4
RUN cd / && git clone https://github.com/lofar-astron/RMextract.git && cd RMextract && export PYTHONPATH=/opt/lofarsoft/lib/python2.7/site-packages/ && python setup.py install --prefix=/opt/lofarsoft && cd / && rm -rf RMextract

# install pyBDSF v1.9.2
RUN cd / && git clone https://github.com/lofar-astron/PyBDSF.git && export PYTHONPATH=/opt/lofarsoft/lib/python2.7/site-packages/ && cd PyBDSF && git checkout v1.9.2 && python setup.py install --prefix=/opt/lofarsoft && cd / && rm -rf PyBDSF

# install LSMTool v1.4.2
RUN cd / && git clone https://github.com/darafferty/LSMTool.git && cd LSMTool && git checkout v1.4.2 && export PYTHONPATH=/opt/lofarsoft/lib/python2.7/site-packages/ && python setup.py install --prefix=/opt/lofarsoft && cd / && rm -rf LSMTool

# install the LOFAR beam library
RUN cd / && git clone https://github.com/lofar-astron/LOFARBeam.git && cd LOFARBeam && mkdir build && cd build && cmake -DPORTABLE=True -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft -DCMAKE_PREFIX_PATH=/opt/lofarsoft ../ && make && make install && cd ../../ && rm -rf LOFARBeam

# install IDG 0.7
RUN cd / && git clone https://gitlab.com/astron-idg/idg.git && cd idg && git checkout 0.7 && mkdir build && cd build && cmake -DPORTABLE=True -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft ../ && make && make install && cd / && rm -rf idg

# install WSClean v2.10.1
RUN cd / && git clone https://gitlab.com/aroffringa/wsclean.git && cd wsclean && git checkout v2.10.1 && mkdir -p build && cd build && cmake -DPORTABLE=True -DCASACORE_ROOT_DIR=/opt/lofarsoft -DIDGAPI_LIBRARIES=/opt/lofarsoft/lib/libidg-api.so -DIDGAPI_INCLUDE_DIRS=/opt/lofarsoft/include -DCMAKE_PREFIX_PATH=/opt/lofarsoft -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft ../ && make -j 4 && make install

# install Dysco v1.2
RUN cd / && wget https://github.com/aroffringa/dysco/archive/v1.2.tar.gz && tar xvf v1.2.tar.gz && cd dysco-1.2 && mkdir build && cd build && cmake -DPORTABLE=True -DCASACORE_ROOT_DIR=/opt/lofarsoft -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft ../ && make -j 4 && make install && cd / && rm -rf v1.2.tar.gz dysco-1.2

# install DP3 v4.1 last py2 version compatible
RUN cd / && git clone https://github.com/lofar-astron/DP3.git && cd DP3 && git checkout v4.1 && mkdir build && cd build && cmake -DPORTABLE=True -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft -DCMAKE_PREFIX_PATH=/opt/lofarsoft ../ && make && make install && ln -s /opt/lofarsoft/bin/DPPP /opt/lofarsoft/bin/NDPPP && cd ../../ && rm -rf DP3

# install LOFAR 3.2.1
RUN cd / && svn --non-interactive -q co https://svn.astron.nl/LOFAR/branches/LOFAR-Release-3_2/ source
RUN cd source && mkdir -p build/gnucxx11_optarch && cd build/gnucxx11_optarch && cmake -DPORTABLE=True -DBUILD_PACKAGES="Pipeline ParmDB pyparmdb" -DCMAKE_INSTALL_PREFIX=/opt/lofarsoft -DCASACORE_ROOT_DIR=/opt/lofarsoft -DBUILD_TESTING=OFF -DUSE_OPENMP=True ../../ && make -j1 && make install && cd / && rm -rf source

# install factor 2
RUN cd / && git clone https://github.com/lofar-astron/factor.git && cd factor && export PYTHONPATH=/opt/lofarsoft/lib/python2.7/site-packages/ && python setup.py install --prefix=/opt/lofarsoft && cd / && rm -rf factor

# setup environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONPATH=/opt/lofarsoft/lib/python2.7/site-packages/
ENV LD_LIBRARY_PATH=/opt/lofarsoft/lib:/usr/local/lib/
# help: Docker image containing all required software to run factor 2.0. This image is built on Ubuntu 18.04.
