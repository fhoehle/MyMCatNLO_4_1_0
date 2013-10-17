#! /bin/bash
installDir=$PWD/MCatNLO_4_1_0_Base
tarBallDir=$PWD/mcatnlo_410_tarBall
mkdir $tarBallDir
cd $tarBallDir
wget http://www.hep.phy.cam.ac.uk/theory/webber/MCatNLO/Package4.10_dist.tar.gz
mkdir $installDir
cd $installDir
tar -xzf $tarBallDir/Package4.10_dist.tar.gz

