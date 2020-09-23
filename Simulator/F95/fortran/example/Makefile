MODELS = InnerBeltModel ReadModelFromFile

test: clean $(MODELS)

$(MODELS):
	$(MAKE) -C .. Driver GIMEOBJ=$@
	@./test.sh ../$@

.PHONY: $(MODELS)

.DEFAULT:

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG

