#! /bin/bash
mainDir=$PWD
installDir=$mainDir/MCatNLO_4_1_0_Base
tarBallDir=$mainDir/mcatnlo_410_tarBall
if [ "X$1" == "X" ]; then 
  mkdir $tarBallDir
  cd $tarBallDir
  wget http://www.hep.phy.cam.ac.uk/theory/webber/MCatNLO/Package4.10_dist.tar.gz
  mkdir $installDir
  cd $installDir
  tar -xzf $tarBallDir/Package4.10_dist.tar.gz
  patch < $mainDir/cmsswMCatNLO_4_1_0.patch
  chmod u+x MCatNLO.inputs
  chmod -R a-w $installDir
else if [  "X$1" == "Xclean" ]; then
  chmod -R u+w $installDir
  rm -rf $tarBallDir
  rm -rf $installDir
else
  echo "option $1 not supported"
  exit 1
fi
fi  
