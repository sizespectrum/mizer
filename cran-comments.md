## Test environments
  
### local Windows 10, R 4.2.2

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓
  
### Winbuilder: R Under development (unstable) (2022-12-21 r83491 ucrt)

OK

### R-hub builder: Windows Server 2022, R-devel, 64 bit

OK

### GitHub actions R-CMD-check

OK on all 5 platforms

### RHub builder: Ubuntu Linux 20.04.1 LTS, R-release, GCC
  
NOTES
  
* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Gustav Delius <gustav.delius@york.ac.uk>’

  Found the following (possibly) invalid DOIs:
    DOI: 10.1111/2041-210X.12256
      From: inst/CITATION
      Status: Service Unavailable
      Message: 503
    
* checking examples ... [36s/146s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                              user system elapsed
  plotBiomassObservedVsModel 2.464  0.004   9.891
  plotYieldObservedVsModel   1.928  0.008   7.906
  addSpecies                 1.658  0.011   7.011
  newCommunityParams         1.528  0.004   6.249
  initialN-set               1.262  0.015   5.592

I have checked the DOI at https://www.doi.org/ and it is correct.

### R-hub builder: Fedora Linux, R-devel, clang, gfortran

* checking examples ... [36s/144s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                              user system elapsed
  plotBiomassObservedVsModel 2.502  0.000   9.975
  plotYieldObservedVsModel   1.999  0.002   7.954
  addSpecies                 1.574  0.010   6.384
  newCommunityParams         1.562  0.000   6.439
  
* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable
  


