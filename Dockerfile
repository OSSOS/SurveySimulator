# this builds a container that can be used to run the SurveySimulator (python and fortran)
# This container is loaded into the canfar Science Portal for use/execution.
# Can also be used directly with docker.
FROM condaforge/miniforge3:latest as base
# FROM jupyter/scipy-notebook as base
# USER root
# ENV DEBIAN_FRONTEND="noninteractive"
RUN apt -y -q update 
RUN apt -y -q install curl wget man man-db git build-essential zip unzip xdg-utils less emacs nano xterm vim rsync tree gfortran
RUN apt -y install python3-numpy
# Python package build uses setuptools + f90wrap (pip install . → setup.py → make -C F95).
# meson/ninja are not required for the SurveySimulator build.



# SKAHA system settings and permissions
RUN apt install -y -q sssd libnss-sss libpam-sss
COPY etc/nofiles.conf /etc/security/limits.d/
COPY etc/nsswitch.conf /etc/
## see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN touch /etc/sudo.conf && echo "Set disable_coredump false" > /etc/sudo.conf
# generate missing dbus uuid (issue #47)
RUN dbus-uuidgen --ensure


# setup this container for skaha launching
COPY etc/startup.sh /skaha/startup.sh
RUN chmod +x /skaha/startup.sh


# setup a the needed python environment
# RUN apt-get update && yes | apt-get install python3.11 pip
# RUN yes | apt install python3.12-venv
# RUN python3 -m venv /opt/SSim/venv
RUN pip install jupyter
RUN pip install cadctap
RUN pip install vos
RUN pip install scipy
RUN pip install astropy
RUN pip install astroquery
RUN pip install matplotlib
RUN pip install f90wrap
# RUN pip install git+https://github.com/jameskermode/f90wrap
RUN pip install rebound
RUN pip3 install astroplan
RUN pip install Deprecated
RUN pip install canfar


# Build the SSim
RUN mkdir -p /opt/SSim
COPY ./ /opt/SSim/
# COPY src /opt/SSim/python
WORKDIR /opt/SSim/F95

# install Fortran based binary of SSim
RUN make clean && ls && make Driver GIMEOBJ=ReadModelFromFile
RUN cp Driver /usr/local/bin/SSim
# RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections 
# RUN apt-get install -y ttf-mscorefonts-installer

# install the Python based version of SSim
FROM base as deploy
WORKDIR /opt/SSim/
RUN pip install .
# RUN python setup.py install

# Two build sets, deploy and test
FROM base as test

RUN mkdir -p /arc/home
RUN groupadd -g 1001 testuser
RUN useradd -u 1001 -g 1001 -s /bin/bash -d /arc/home/testuser -m testuser
RUN chown -R testuser /opt/SSim
WORKDIR /opt/SSim/
# RUN pip3 install -e .
USER testuser
WORKDIR /arc/home/testuser
COPY etc/ReadModelFromFile.in ./
ENTRYPOINT ["/skaha/startup.sh"]
