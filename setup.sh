# setting up neccessary libraries
command -v scram > /dev/null 2>&1
if [ "X$SCRAM_ARCH" -ne "slc6_amd64_gcc462" -a "X$SCRAM_ARCH" -ne "slc5_amd64_gcc462" ]; then
  echo "wrong architecture or none given SCRAM_ARCH: $SCRAM_ARCH"
  exit 1
fi
if [ $? -ne 0 ]; then
  echo "scarm doesn't exist, run cmsenv"
  exit 1
fi
if [ -z "$LHAPDF_ROOT" -o "$1" == "-f" ]; then
  echo "sourcing lhapdf from scram tool info lhapdf"
  LHA_BASE=`scram tool info lhapdf | grep LHAPDF_BASE | sed 's/.*=\(.*\)/\1/'`
  source $LHA_BASE/etc/profile.d/init.sh
else
  echo "already existing, doing nothing \n resetting it by calling with force -f "
fi
# main MCatNLO source
export LOCAL_MCatNLODIR=$PWD
# make scripts available
export PATH=$PWD/bin:$PATH
