
README4 for the OSSOSv4 release 
--------------------------------

Web link to the release is at:
https://wiki.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/ossos/index.php/OSSOSv4_data_release

The README files within this release tree are NOT a substitute for the
extensive documentation on the wiki pages!

--------------------------
HISTORY
19 Sep 2013	Preliminary V1 release
10 Mar 2014	V2 release, E block stage 2 release + O block stage 'almost 1'
 7 May 2014	V3 release, both E block and O block are 'stage 2'
28 Jan 2015     V4 release, E block and O block at 'stage 3'
?? Feb 2015     coming V5 release.  L block added stage 2
?? Mar 2015     coming V6 release.  H block added stage 1
--------------------------

Contents
---------------------------------------------------------------------------

Main release files (will always be present) -----

READMEv4.txt   		this file
OSSOSv4.detections	with orbital elements and classes. See below for format.
                        Contains only E and O block
OSSOSv4ScatterPlots.pdf	Plot of uncertainties and orbital elements
mpc             subdir	MPC format files for E+O detections
ast             subdir	OSSOS AST format files for E+O (see below)      **New**
SurvSimv4       subdir  Survey Simulator for OSSOS E+O, and/or CFEPS

Non-standard files in this release ---------
    none, note that the "ast" subdirectory is new, but will be standard

Reminders:
E block : Straddles "E"cliptic near RA=14h
O block : "O"ff ecliptic block near RA=15h

------------------------------------------------------------------------------

File Formats

The *detections* file gives an averaged magnitude (and band, and uncertainty) for each
    object DURING THE DISCOVERY NIGHT'S TRIPLE (this is the only magnitude that
    matters for the characterization).  The 'surmised' absolute magnitude (that
    is, the H_r magnitude surmised using the average m_r and the discovery geometry)
    is also given (and has the same error). The dynamical classification of the objects
    is given in the first several columns; see the wiki page for details.

Caveats for V4 release
------------------------------------------------------------
- The 'characterization mag' cut is r=24.04 for E block, 24.40 for O block
- There are 50/36 characterized objects in the E/O block, respectively.  The
    unclassified uo3eXX and uo3oXX objects are not part of the release but
    can be obtained by request to ossoscore@gmail.com
- Orbital classifications are mostly now 'secure' according to SSBN08 
    nomenclature. The class file gives an 'S'=secure or 'I'=insecure rating
- The 'ast' subdirectory is a set of files which contain information
  useful to those wanting photometric information.  The meanings of
  these colums are documented at
https://wiki.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/ossos/index.php/MPC_codes_used_in_flagging_astrometry#.ast_suffix_files
