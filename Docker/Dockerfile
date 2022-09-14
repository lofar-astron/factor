#############################################################################
#  First stage of the build: compile and install the required software.     #
#  This stage requires quite a lot of packages from the repository.         #
#############################################################################

FROM ubuntu:18.04 as builder

LABEL Maintainer="loose@astron.nl"

ARG BUILD_DIR=/build
ARG INSTALL_DIR=/opt/lofarsoft
ARG PORTABLE=TRUE

# Install all build-time dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update && \
    apt-get install -y \
        bison \
        build-essential \
        casacore-data \
        casacore-dev \
        casacore-tools \
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
        liblog4cplus-dev \
        liblua5.3-dev \
        libpng-dev \
        libxml2-dev \
        ninja-build \
        python-casacore \
        python-dev \
        python-matplotlib \
        python-pip \
        subversion \
        wcslib-dev \
        wget

WORKDIR ${BUILD_DIR}

#===========================================================================#
#                          C + +   P A C K A G E S                          #
#===========================================================================#

# install AOFlagger 2.14 last py2 version compatible
RUN git clone --depth 1 --branch v2.14.0 https://gitlab.com/aroffringa/aoflagger.git
# Remove optional (unused) Fortran support to support building wiht Ninja.
RUN sed -i '/^enable_language(Fortran OPTIONAL)$/d' ${BUILD_DIR}/aoflagger/CMakeLists.txt
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -H${BUILD_DIR}/aoflagger \
    -B${BUILD_DIR}/aoflagger/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/aoflagger/build install

# install the LOFAR beam library (not sure why we use HEAD of master here?)
RUN git clone --depth 1 https://github.com/lofar-astron/LOFARBeam.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -H${BUILD_DIR}/LOFARBeam \
    -B${BUILD_DIR}/LOFARBeam/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/LOFARBeam/build install

# install IDG 0.7
RUN git clone --depth 1 --branch 0.7 https://gitlab.com/astron-idg/idg.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -H${BUILD_DIR}/idg \
    -B${BUILD_DIR}/idg/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/idg/build install

# install WSClean v2.10.1
RUN git clone --depth 1 --branch v2.10.1 https://gitlab.com/aroffringa/wsclean.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -H${BUILD_DIR}/wsclean \
    -B${BUILD_DIR}/wsclean/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/wsclean/build install

# install Dysco v1.2
RUN git clone --depth 1 --branch v1.2 https://github.com/aroffringa/dysco.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -H${BUILD_DIR}/dysco \
    -B${BUILD_DIR}/dysco/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/dysco/build install

# install DP3 v4.1 last py2 version compatible
RUN git clone --depth 1 --branch v4.1 https://github.com/lofar-astron/DP3.git
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -H${BUILD_DIR}/DP3 \
    -B${BUILD_DIR}/DP3/build \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/DP3/build install

# install LOFAR 3.2.1
RUN svn --non-interactive --quiet checkout \
    https://svn.astron.nl/LOFAR/branches/LOFAR-Release-3_2
# Remove unnecessary requirement, we're building without tests
RUN sed -i '/^find_python_module(xmlrunner REQUIRED)$/d' \
    ${BUILD_DIR}/LOFAR-Release-3_2/CEP/Pipeline/framework/CMakeLists.txt
RUN cmake \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DCMAKE_PREFIX_PATH=${INSTALL_DIR} \
    -DPORTABLE=${PORTABLE} \
    -DBUILD_PACKAGES="Pipeline ParmDB pyparmdb" \
    -DBUILD_TESTING=OFF \
    -DUSE_OPENMP=TRUE \
    -H${BUILD_DIR}/LOFAR-Release-3_2 \
    -B${BUILD_DIR}/LOFAR-Release-3_2/build/gnucxx11_optarch \
    -G Ninja
RUN ninja -C ${BUILD_DIR}/LOFAR-Release-3_2/build/gnucxx11_optarch install

# Cleanup stuff we don't need at run-time
RUN rm -rf ${INSTALL_DIR}/include ${INSTALL_DIR}/sbin

#===========================================================================#
#                       P Y T H O N   P A C K A G E S                       #
#===========================================================================#

# Set-up environment for python
ENV PATH="${INSTALL_DIR}/bin:${PATH}"
ENV PYTHONUSERBASE=${INSTALL_DIR}

# install python2 base packages
RUN pip install --user --upgrade pip
RUN pip install --user --upgrade setuptools wheel

# # install python-casacore
# RUN pip install --user git+https://github.com/casacore/python-casacore.git@v3.2.0

# install losoto 2.0
RUN pip install --user numexpr==2.7.3  # avoids install of incompatible 2.8
RUN pip install --user git+https://github.com/revoltek/losoto.git@2.0

# install RMextract v0.4
RUN pip install --user git+https://github.com/lofar-astron/RMextract.git@v0.4

# install pyBDSF v1.9.2
RUN pip install --user git+https://github.com/lofar-astron/PyBDSF.git@v1.9.2

# install LSMTool v1.4.2
RUN pip install --user git+https://github.com/darafferty/LSMTool.git@v1.4.2

# install factor 2
RUN pip install --user git+https://github.com/lofar-astron/factor.git


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
        libcasa-casa2 \
        libcasa-fits2 \
        libcasa-measures2 \
        libcasa-ms2 \
        libcasa-python2 \
        libcasa-scimath2 \
        libcasa-tables2 \
        libcfitsio5 \
        libfftw3-double3 \
        libfftw3-single3 \
        libgcc1 \
        libgfortran4 \
        libglibmm-2.4-1v5 \
        libgomp1 \
        libgsl23 \
        libgslcblas0 \
        libgtkmm-3.0-1v5 \
        libhdf5-100 \
        libhdf5-cpp-100 \
        liblapack3 \
        liblog4cplus-1.1-9 \
        libpng16-16 \
        libpython2.7 \
        libsigc++-2.0-0v5 \
        libstdc++6 \
        libxml2 \
        python \
        python-casacore \
        python-matplotlib \
        wget

RUN rm -rf /var/lib/apt/lists/*

# # Install WSRT Measures (extra casacore data). We purposely do this in the
# # `runner` stage and not in the `builder` stage, because otherwise the
# # previous `apt-get install` would clobber the files we had installed in
# # in the `builder` stage (as `casacore-data` is installed as a dependency).
# # Note: The file on the ftp site is updated daily. When warnings regarding
# # leap seconds appear, ignore them or regenerate the docker image.
# RUN wget -q -O /WSRT_Measures.ztar \
#         ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar && \
#     cd /var/lib/casacore/data && \
#     tar xfz /WSRT_Measures.ztar && \
#     rm /WSRT_Measures.ztar

# Avoid concurrency issues when accessing HDF5 files from Python
ENV HDF5_USE_FILE_LOCKING=FALSE

# Create entry point for container
COPY entrypoint.sh /
ENTRYPOINT [ "/entrypoint.sh" ]