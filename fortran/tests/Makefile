MODELS = InnerHotModel ReadModelFromFile HST
LANGUAGE = F95

%: %.sh

test: clean $(MODELS)

$(MODELS): clean binEphem.DE
	$(MAKE) -s -C ../$(LANGUAGE) Driver GIMEOBJ=$@ 
	@./test.sh ../$(LANGUAGE)/$@ $(LANGUAGE) 
	$(MAKE) -s -C ../$(LANGUAGE) clean

.PHONY: $(MODELS)

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test

binEphem.DE:
	\rm -f binEphem.DE
	$(MAKE) -s -C ../DE440 binEphem.DE
	ln -s ../DE440/binEphem.DE binEphem.DE
