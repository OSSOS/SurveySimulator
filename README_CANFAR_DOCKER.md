# Building the SurveySimulator image

To run the SSim on CANFAR we need to make a `docker` container.  There are `dev` and `deploy` targets in the Makefile. 


## Testing/Development build and run
To build a development version of the container that you will run locally to ensure everything is working.  
You can also use this conatiner to run SSim analysis locally.
```
make dev
docker run --user testuser --interactive --tty --ip 0.0.0.0 --rm --env DISPLAY=host.docker.internal:0 images.canfar.net/uvickbos/ssim:python xterm -fg white -bg black -title ssim:python
```
This will launch an xterm running as testuser with access to the Survey Simulator.

### Note on macOS and X11  
To use the above testuser on OS-X requires having X11 running. On OS-X do the following:
- Install XQuartz (likely you already have)
- start and `xterm`
- set XQuartz->Preferences->Security : Allow conections from network clients.
- type `xhost +` in your xterm window to allow open connections to X11

## Production

This is the version we will load to `images.canfar.net`.  The `make` command builds the production versionof the container and pushes to `CANFAR` 
You may need to do a `docker login` before running this build step, see https://github.com/opencadc/skaha/tree/master/skaha-containers#publishing

To make a release to `canfar.net` run in deploy mode.

```
make deploy 
```

### Tag on images.canfar.net ###
Once you have loaded the images log into `images.canfar.net` and tagged them as `desktop-app` images.



