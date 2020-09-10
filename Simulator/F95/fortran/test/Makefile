
all: test

test: 
	@$(MAKE) -C .. Driver GIMEOBJ=InnerBeltModel
	@./test.sh ../InnerBeltModel
	@$(MAKE) -C .. Driver GIMEOBJ=ReadModelFromFile
	@./test.sh ../ReadModelFromFile

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG
