#!/bin/sh

#cd /etc/yum.repos.d
#wget http://people.centos.org/tru/devtools-2/devtools-2.repo

#yum install devtoolset-2-gcc 
#yum install devtoolset-2-binutils 
#yum install devtoolset-2-gcc-gfortran 
#yum install devtoolset-2-gcc-c++

#source /opt/rh/devtoolset-2/enable
#scl enable devtoolset-2 bash
yum -y install centos-release-scl
yum -y install devtoolset-7-gcc devtoolset-8-gcc-c++ devtoolset-8-binutils
yum -y install devtoolset-7-gcc-gfortran 
yum -y install devtoolset-7-gcc-c++
echo "source /opt/rh/devtoolset-7/enable" >> /etc/profile
#echo "export PATH=/opt/rh/devtoolset-7/enable :$PATH">> /etc/profile
source /etc/profile
gcc --version

echo "gcc update successfully"
