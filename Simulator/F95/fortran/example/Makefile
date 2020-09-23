MODELS = InnerBeltModel ReadModelFromFile

%: %.sh

test: clean $(MODELS)

$(MODELS):
	$(MAKE) -C ../ Driver GIMEOBJ=$@
	@./test.sh ../$@

.PHONY: $(MODELS)

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG

