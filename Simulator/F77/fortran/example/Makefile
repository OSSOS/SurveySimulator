MODELS = InnerHotModel ReadModelFromFile

test: clean $(MODELS)

$(MODELS):
	$(MAKE) -C .. Driver GIMEOBJ=$@
	@./test.sh ../$@

.PHONY: $(MODELS)

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG
