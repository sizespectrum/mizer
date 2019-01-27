## Test environments
* ubuntu 18.04 bionic, R 3.5.2
* win-builder (devel and release)
* R-hub

## Errors
None.

## Warnings
* r-hub Fedora Linux, R-devel, clang, gfortra: re-building the
  vignette with knitr fails with message "there is no package called codetools".
  I think this must be a bug with how knitr is set up on r-hub with Fedora.
  
## Notes
* r-hub Windows Server 2008 R2 SP1, R-devel, 32/64 bit: can not check the
  sizes of PDF files because GhostScript executables were missing.
  
* r-hub Ubuntu Linux 16.04 LTS, R-release, GCC: flags a
  misspelling which is however intentional


