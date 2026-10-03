# Building the SurveySimulator image

A single Docker image serves every runtime:

- **CANFAR / skaha** — the same image runs as a **notebook** (JupyterLab), a
  **desktop-app** (`xterm` with `SSim`/`Driver` on `PATH`), and a **headless**
  batch session.
- **Cursor cloud** — [`.cursor/environment.json`](../.cursor/environment.json)
  builds this image and does an editable install of the `/workspace` checkout.
- **Local development** — [`.devcontainer/devcontainer.json`](../.devcontainer/devcontainer.json)
  builds the same image on a laptop (Dev Containers).

The image is built by extending CANFAR's `astroml` base
(`images.canfar.net/skaha/astroml`), which already provides conda Python 3.12,
JupyterLab, `xterm`, the CADC client tools, and the SSS user mapping skaha
relies on. The [`Dockerfile`](../Dockerfile) adds the Fortran toolchain, bakes
the Survey Simulator, and sets an explicit `tini` + [`etc/startup.sh`](../etc/startup.sh)
entrypoint (skaha overrides the container `CMD` per session type).

## Build

```
make build
```

Override the base image snapshot with `make build ASTROML_TAG=26.04` (defaults
to a pinned monthly tag).

## Run locally

Run the image as a desktop-app `xterm` (the simulator binary is on `PATH` as `SSim`):

```
make dev
docker run --rm -it images.canfar.net/uvickbos/ssim:<version> xterm
```

Or start a shell to run the Fortran `Driver` / Python `ossssim` directly:

```
docker run --rm -it images.canfar.net/uvickbos/ssim:<version> bash
```

## Production

Build and push to `images.canfar.net`. Run `docker login images.canfar.net`
first (see https://github.com/opencadc/science-containers for Harbor access).

```
make deploy
```

### Tag on images.canfar.net

Once loaded, log into `images.canfar.net` and tag the image with the session
types it supports (`notebook`, `desktop-app`, and `headless`).
