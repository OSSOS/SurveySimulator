M	= gimeobjut
SRCM	= $(addsuffix .f95,${M})

SSIM	= SurveySubsF95
GOBJ	= GiMeObjF95

MODEL	= Res32Model.f95
#MODEL	= Res32Model.f95

clean:
	rm -f *~ $(SSIM).py $(GOBJ).py _$(SSIM)*.so _$(GOBJ)*.so

Modules:
	\rm -f ../srcF95/$(SRCM)
	\ln -s ../PythonF95/$(MODEL) ../srcF95/$(SRCM)
	(cd ../srcF95/; make Modules)
	\mv ../srcF95/$(SSIM).py ../srcF95/_$(SSIM).so ../srcF95/$(GOBJ).py ../srcF95/_$(GOBJ).so .
	\rm ../srcF95/$(SRCM)
	echo "Modules have been built"

test: Modules
	./test.sh
