if [ -f /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/etc/profile.d/dependencies-setup.sh ]; then . /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/etc/profile.d/dependencies-setup.sh; fi
HECTOR_ROOT="/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed"
HECTOR_VERSION="1.3.4-ikhhed"
HECTOR_REVISION="1"
HECTOR_CATEGORY="external"
[ ! -d /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/bin ] || export PATH="/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/bin${PATH:+:$PATH}";
[ ! -d /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/lib ] || export LD_LIBRARY_PATH="/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/hector/1.3.4-ikhhed/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}";

