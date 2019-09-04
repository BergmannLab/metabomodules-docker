#! /bin/bash -

# Copyright 2019 Bergmann's Lab UNIL <mattia.tomasoni@unil.ch>
#
#    This file is part of METABOMODULES Tool.
#
#    METABOMODULES Tool is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    METABOMODULES Tool is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    GNU General Public License: <https://www.gnu.org/licenses/>.

echo 
echo "Installing METABOMODULES"
echo

# ask for root password (needed to make monet command available from all
# locations by adding link in /usr/local/bin/monet)
echo "Superuser rights are required."
sudo ls > /dev/null
echo

# check whether metabomodules is already installed
if [ -f /usr/local/bin/metabomodules ]; then
  read -p "metabomodules is already installed. Would you like to overwrite? [y|n] " -n 1 -r
  echo ""
  if [[ $REPLY =~ ^[y]$ ]]; then
    ./uninstall.sh > /dev/null 2>&1
  elif [[ $REPLY =~ ^[n]$ ]]; then
    echo "EXITING: metabomodules WAS NOT RE-INSTALLED."
    exit 0
  else
    echo "invalid option selected"
    echo "EXITING: metabomodules WAS NOT RE-INSTALLED."
    exit 0
  fi  
fi

# check that either docker or singularity are installed
echo "- Checking that docker and/or singularity are installed..."
sudo bash -c "docker --help > /tmp/docker_test 2>&1" 2>&1 /dev/null
if [ $? -eq "0" ]; then docker_installed=true; else docker_installed=false; fi
sudo bash -c "singularity --help > /tmp/singularity_test"  2>&1 /dev/null
if [ $? -eq "0" ]; then singularity_installed=true; else singularity_installed=false; fi

if ! $docker_installed && ! $singularity_installed; then
  echo "  WARNING: Neither docker nor singularity seem to be installed."
  echo "    Install either of the two to successfully run metabomodules."
  echo "    Please visit https://www.docker.com or http://singularity.lbl.gov"
else
  echo "  ...OK"
fi

# store metabomodules code in the home directory
echo "- Copying files..."
mkdir ~/.metabomodules
cp -r ./* ~/.metabomodules
chmod -R 750 ~/.metabomodules
echo "  ...OK"

# make metabomodules command available
echo "- Updating operating system..."
sudo ln -s ~/.metabomodules/metabomodules /usr/local/bin/metabomodules
echo "  ...OK"

echo "" && echo "FINISHED: metabomodules WAS INSTALLED SUCCESSFULLY."
echo "Invoke metabomodules in a bash shell from any location."
echo 

# reload current shell so changes to ~/.bashrc will become available
exec bash
