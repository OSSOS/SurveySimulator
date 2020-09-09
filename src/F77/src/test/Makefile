
all: test

CMDS := Parametric ReadModelFromFile


test: 
	@$(MAKE) -C .. $(CMDS)
	@./test.sh ../Parametric
	@./test.sh ../ReadModelFromFile

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG
