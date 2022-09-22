## Test environments
  
* local Windows 10, R 4.2.1

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓
  

* rhub::check_for_cran(platforms = "macos-highsierra-release-cran")
  
  SUCCESS

* RStudio Cloud with x86_64-pc-linux-gnu (64-bit), R 4.2.1
* Win builder devel

  0 errors ✔ | 0 warnings ✔ | 1 note ✖
  
  checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Gustav Delius <gustav.delius@york.ac.uk>’
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1111/2041-210X.12256
      From: inst/CITATION
      Status: Service Unavailable
      Message: 503
      
  I have checked the DOI at https://www.doi.org/ and it is correct
