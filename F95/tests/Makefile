MODELS = InnerHotModel ReadModelFromFile
# MODELS = ReadModelFromFile
LANGUAGE = F95

%: %.sh

test: clean $(MODELS) HST test_eta_rate test_polygon_area test_pointing_geom

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

# Unit test for spherical polygon_area
test_polygon_area:
	$(MAKE) -s -C ../ libsrc.a
	gfortran -O0 -I../ -o test_polygon_area test_polygon_area.f95 -L../ -lsrc
	./test_polygon_area
	\rm -f test_polygon_area
	$(MAKE) -s -C ../ rmtmp

# Unit test for survey_load / pointing_geom against CFEPS
test_pointing_geom:
	$(MAKE) -s -C ../ libsrc.a
	gfortran -O0 -I../ -o test_pointing_geom test_pointing_geom.f95 -L../ -lsrc
	./test_pointing_geom
	\rm -f test_pointing_geom
	$(MAKE) -s -C ../ rmtmp

.PHONY: $(MODELS) test_eta_rate test_polygon_area test_pointing_geom

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test test_polygon_area test_pointing_geom
