MODELS = InnerHotModel ReadModelFromFile
LANGUAGE = F95

%: %.sh

test: clean $(MODELS)

$(MODELS):
	$(MAKE) -C ../ Driver GIMEOBJ=$@
	@./test.sh ../$@ $(LANGUAGE)

.PHONY: $(MODELS)

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test

