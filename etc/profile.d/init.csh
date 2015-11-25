if ( -f /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/etc/profile.d/dependencies-setup.csh ) source /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/etc/profile.d/dependencies-setup.csh; endif
set HECTOR_ROOT="/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed"
set HECTOR_VERSION="1.3.4-ikhhed"
set HECTOR_REVISION="1"
set HECTOR_CATEGORY="external"
if ( -d /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/bin ) then
  if ( ${?PATH} ) then
    setenv PATH "/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/bin:$PATH"
  else
    setenv PATH "/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/bin"
  endif
endif
if ( -d /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/lib ) then
  if ( ${?LD_LIBRARY_PATH} ) then
    setenv LD_LIBRARY_PATH "/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/lib:$LD_LIBRARY_PATH"
  else
    setenv LD_LIBRARY_PATH "/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/lib"
  endif
endif

