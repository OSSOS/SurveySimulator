# Unified CANFAR/skaha + Cursor/dev container for the OSS Survey Simulator.
#
# One image serves every runtime:
#   * CANFAR/skaha  - notebook (JupyterLab), desktop-app (xterm) and headless sessions
#   * Cursor cloud  - .cursor/environment.json builds and tests the /workspace checkout
#   * Local dev     - .devcontainer/devcontainer.json builds the same image on a laptop
#
# Base: CANFAR's astroml image is explicitly a "same container, different interfaces"
# image - it already ships conda Python 3.12, JupyterLab, xterm, the CADC client tools
# and the SSS user mapping skaha relies on, and runs in all three skaha session types.
# We extend it rather than rebuild that stack (CANFAR guidance: always extend the base).
# https://www.opencadc.org/canfar/latest/platform/containers/
ARG ASTROML_TAG=26.06
FROM images.canfar.net/skaha/astroml:${ASTROML_TAG}

LABEL maintainer="J.J. Kavelaars <jjkavelaars@gmail.com>"

USER root

# conda is the runtime Python for every session type; put it first on PATH for all
# users and non-login shells (skaha, the ubuntu dev user, and the Cursor agent).
ENV PATH=/opt/conda/bin:${PATH}

# Fortran toolchain for the F95 detection engine, plus tini for signal/zombie handling.
# astroml may already provide some of these; apt is idempotent and only adds what is missing.
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential gfortran make tini \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /var/tmp/*

# Bake the stable Survey Simulator into the image for the CANFAR runtime
# (CANFAR convention: ship tested code in the image). This installs the Python
# package (ossssim + the f90wrap ossssimlib extension) into the conda env and builds
# the Fortran Driver, exposed on PATH as `SSim`.
RUN mkdir -p /opt/SSim
COPY . /opt/SSim/
WORKDIR /opt/SSim
RUN pip install . \
 && make -C F95 clean \
 && make -C F95 Driver GIMEOBJ=ReadModelFromFile \
 && cp F95/Driver /usr/local/bin/SSim

# Development user for Cursor cloud and the local devcontainer. skaha injects the real
# CADC user at runtime via SSS, so this `ubuntu` user is only used by Cursor/devcontainer.
# Passwordless sudo lets the editable (-e) install write into the root-owned conda env
# without a multi-GB `chown -R /opt/conda` layer.
RUN if ! id -u ubuntu >/dev/null 2>&1; then useradd -m -s /bin/bash ubuntu; fi \
 && echo "ubuntu ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-ubuntu \
 && chmod 0440 /etc/sudoers.d/90-ubuntu

# Our own explicit, auditable launch path. skaha overrides CMD per session type
# (notebook/desktop-app/headless); tini reaps processes and startup.sh execs that command.
RUN mkdir -p /skaha
COPY etc/startup.sh /skaha/startup.sh
RUN chmod +x /skaha/startup.sh

WORKDIR /
ENTRYPOINT ["tini", "-g", "--", "/skaha/startup.sh"]
