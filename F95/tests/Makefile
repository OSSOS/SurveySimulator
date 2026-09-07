MODELS = InnerHotModel ReadModelFromFile
# MODELS = ReadModelFromFile
LANGUAGE = F95

%: %.sh

test: clean $(MODELS) HST test_eta_rate

$(MODELS): clean
	$(MAKE) -s -C ../ Driver GIMEOBJ=$@ 
	@./test.sh ../$@ $(LANGUAGE) 
	$(MAKE) -s -C ../ clean

HST: clean
	$(MAKE) -s -C ../ Driver GIMEOBJ=ReadModelFromFile
	@./test.sh ../$@ $(LANGUAGE) 
	$(MAKE) -s -C ../ clean

# Unit test for optional rate-dependent m0 in single_param / eta_raw
test_eta_rate:
	$(MAKE) -s -C ../ libsrc.a
	gfortran -O0 -I../ -o test_eta_rate test_eta_rate.f95 -L../ -lsrc
	./test_eta_rate
	\rm -f test_eta_rate
	$(MAKE) -s -C ../ clean

.PHONY: $(MODELS) test_eta_rate

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test test_eta_rate
