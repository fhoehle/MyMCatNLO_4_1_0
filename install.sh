#! /bin/bash
mainDir=$PWD
installDir=$mainDir/MCatNLO_4_1_0_Base
tarBallDir=$mainDir/mcatnlo_410_tarBall
mkdir $tarBallDir
cd $tarBallDir
wget http://www.hep.phy.cam.ac.uk/theory/webber/MCatNLO/Package4.10_dist.tar.gz
mkdir $installDir
cd $installDir
tar -xzf $tarBallDir/Package4.10_dist.tar.gz
patch < $mainDir/cmsswMCatNLO_4_1_0.patch
