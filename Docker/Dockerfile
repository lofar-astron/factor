#############################################################################
#  First stage of the build: compile and install the required software.     #
#  This stage requires quite a lot of packages from the repository.         #
#############################################################################

FROM ubuntu:18.04 as builder

LABEL Maintainer="loose@astron.nl"

ARG BUILD_DIR=/build
ARG INSTALL_DIR=/opt/lofarsoft
ARG PORTABLE=ON

# Install all build-time dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update && \
    apt-get install -y \
        bison \
        build-essential \
        cmake \
        flex \
        gfortran \
        git \
        libarmadillo-dev \
        libboost-date-time-dev \
        libboost-filesystem-dev \
        libboost-numpy-dev \
        libboost-program-options-dev \
        libboost-python-dev \
        libboost-system-dev \
        libboost-test-dev \
        libcfitsio-dev \
        libfftw3-dev \
        libgsl-dev \
        libgtkmm-3.0-dev \
        libhdf5-serial-dev \
        liblua5.3-dev \
        libpng-dev \
        libxml2-dev \
        ninja-build \
        python-dev \
        python-matplotlib \
        python-numpy \
        python-pip \
        wcslib-dev \
        wget

WORKDIR ${BUILD_DIR}

#===========================================================================#
#                          C + +   P A C K A G E S                          #
#===========================================================================#

# Install casacore data
RUN mkdir -p ${INSTALL_DIR}/share/casacore/data
RUN wget ftp://anonymous@ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar
RUN tar -C ${INSTALL_DIR}/share/casacore/data -x -f WSRT_Measures.ztar

# install casacore v2.4.1, last release using casa namespace
RUN git clone --depth 1 --branch v2.4.1 https://github.com/casacore/casacore.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_PYTHON=ON \
    -DBUILD_TESTING=OFF \
    -DUSE_FFTW3=ON \
    -DUSE_HDF5=ON \
    -DUSE_OPENMP=ON \
    -DUSE_THREADS=ON \
    -H${BUILD_DIR}/casacore \
    -B${BUILD_DIR}/casacore/build
RUN make -j $(nproc) -C ${BUILD_DIR}/casacore/build install

# Install AOFlagger 2.14, last py2 compatible release
RUN git clone --depth 1 --branch v2.14.0 https://gitlab.com/aroffringa/aoflagger.git
# Remove optional (unused) Fortran support to support building wiht Ninja.
RUN sed -i '/^enable_language(Fortran OPTIONAL)$/d' ${BUILD_DIR}/aoflagger/CMakeLists.txt
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_TESTING=OFF \
    -H${BUILD_DIR}/aoflagger \
    -B${BUILD_DIR}/aoflagger/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/aoflagger/build install

# Install the LOFAR beam library (not sure why we use HEAD of master here?)
RUN git clone --depth 1 https://github.com/lofar-astron/LOFARBeam.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_TESTING=OFF \
    -H${BUILD_DIR}/LOFARBeam \
    -B${BUILD_DIR}/LOFARBeam/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/LOFARBeam/build install

# Install IDG 0.8, last py2 compatible release
RUN git clone --depth 1 --branch 0.8 https://gitlab.com/astron-idg/idg.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_TESTING=OFF \
    -H${BUILD_DIR}/idg \
    -B${BUILD_DIR}/idg/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/idg/build install

# Install WSClean v2.10.1, last py2 compatible release
RUN git clone --depth 1 --branch v2.10.1 https://gitlab.com/aroffringa/wsclean.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_TESTING=OFF \
    -H${BUILD_DIR}/wsclean \
    -B${BUILD_DIR}/wsclean/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/wsclean/build install

# Install Dysco v1.2
RUN git clone --depth 1 --branch v1.2 https://github.com/aroffringa/dysco.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_TESTING=OFF \
    -H${BUILD_DIR}/dysco \
    -B${BUILD_DIR}/dysco/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/dysco/build install

# Install DP3 v4.1, last py2 compatible release
RUN git clone --depth 1 --branch v4.1 https://github.com/lofar-astron/DP3.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_TESTING=OFF \
    -H${BUILD_DIR}/DP3 \
    -B${BUILD_DIR}/DP3/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/DP3/build install

# Install LOFAR 3.2.18, last release using casa namespace
RUN git clone --depth 1 --branch LOFAR-Release-3_2_18 \
    https://git.astron.nl/ro/lofar.git
# Remove unnecessary requirement, we're building without tests
RUN sed -i '/^find_python_module(xmlrunner REQUIRED)$/d' \
    ${BUILD_DIR}/lofar/CEP/Pipeline/framework/CMakeLists.txt
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_PACKAGES="Pipeline ParmDB pyparmdb" \
    -DBUILD_TESTING=OFF \
    -DUSE_LOG4CPLUS=OFF \
    -DUSE_OPENMP=ON \
    -H${BUILD_DIR}/lofar \
    -B${BUILD_DIR}/lofar/build/gnucxx11_optarch \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/lofar/build/gnucxx11_optarch install


#===========================================================================#
#                       P Y T H O N   P A C K A G E S                       #
#===========================================================================#

# Set-up environment for python
ENV PATH="${INSTALL_DIR}/bin:${PATH}"
ENV PYTHONUSERBASE=${INSTALL_DIR}

# Install python2 base packages
RUN pip install --user --upgrade pip
RUN pip install --user --upgrade setuptools wheel

RUN pip install --user --upgrade \
    numexpr==2.7.3 \
    python-monetdb \
    shapely \
    wcsaxes==0.7

# Install python-casacore v3.2.0
RUN pip install --user --upgrade \
    --global-option build_ext \
    --global-option -I${INSTALL_DIR}/include \
    --global-option -L${INSTALL_DIR}/lib \
    git+https://github.com/casacore/python-casacore.git@v3.2.0

# Install losoto 2.1 pre-release, last py2 compatible commit
RUN pip install --user --upgrade \
    git+https://github.com/revoltek/losoto.git@d35aa97

# Install RMextract v0.4
RUN pip install --user --upgrade \
    git+https://github.com/lofar-astron/RMextract.git@v0.4

# Install pyBDSF v1.9.2
RUN pip install --user --upgrade \
    git+https://github.com/lofar-astron/PyBDSF.git@v1.9.2

# Install LSMTool v1.4.2
RUN pip install --user --upgrade \
    git+https://github.com/darafferty/LSMTool.git@v1.4.2

# Install factor 2
RUN pip install --user --upgrade \
    git+https://github.com/lofar-astron/factor.git

# Cleanup header files, we don't need them at run-time
RUN rm -rf ${INSTALL_DIR}/include


#############################################################################
#  Second stage of the build. Install minimal set of packages, only those   #
#  required at run-time, to reduce the size of the final Docker image.      #
#############################################################################

FROM ubuntu:18.04 as runner

# Arguments are not copied over between stages
ARG INSTALL_DIR=/opt/lofarsoft

# Copy the installed software from the builder stage to the runner stage.
COPY --from=builder ${INSTALL_DIR} ${INSTALL_DIR}

# Install all run-time dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update && \
    apt-get install -y \
        libarmadillo8 \
        libatkmm-1.6-1v5 \
        libboost-date-time1.65.1 \
        libboost-filesystem1.65.1 \
        libboost-numpy1.65.1 \
        libboost-program-options1.65.1 \
        libboost-python1.65.1 \
        libboost-system1.65.1 \
        libcairomm-1.0-1v5 \
        libcfitsio5 \
        libfftw3-double3 \
        libfftw3-single3 \
        libgfortran4 \
        libglibmm-2.4-1v5 \
        libgomp1 \
        libgsl23 \
        libgslcblas0 \
        libgtkmm-3.0-1v5 \
        libhdf5-100 \
        libhdf5-cpp-100 \
        liblapack3 \
        libpng16-16 \
        libpython2.7 \
        libsigc++-2.0-0v5 \
        libstdc++6 \
        libwcs5 \
        libxml2 \
        python \
        python-matplotlib \
        python-numpy

RUN rm -rf /var/lib/apt/lists/*

# Avoid concurrency issues when accessing HDF5 files from Python
ENV HDF5_USE_FILE_LOCKING=FALSE

# Create entry point for container
COPY entrypoint.sh /
ENTRYPOINT [ "/entrypoint.sh" ]
