
all: test

test: clean
	@$(MAKE) -C .. Driver GIMEOBJ=InnerHotModel
	@./test.sh ../InnerHotModel
	@$(MAKE) -C .. Driver GIMEOBJ=ReadModelFromFile
	@./test.sh ../ReadModelFromFile

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG
