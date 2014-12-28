#!/bin/sh

# Install OpenMM to /usr/local/openmm
cwd=`pwd`

cd `dirname $0`
unzip OpenMM-6.2.0-Linux.zip
cd OpenMM-6.2.0-Linux/
sudo mkdir -p /usr/local/openmm

sudo ./install.sh << EOF


EOF

cd $cwd
