# setting up neccessary libraries
LHA_BASE=`scram tool info lhapdf | grep LHAPDF_BASE | sed 's/.*=\(.*\)/\1/'`
source $LHA_BASE/etc/profile.d/init.sh
