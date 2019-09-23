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

# check whether metabomodules is already installed
if [ -d ~/.metabomodules ]; then
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
echo "- Are docker or singularity installed?..."
docker --help > /tmp/docker_test 2>&1
if [ $? -eq "0" ]; then docker_installed=true; else docker_installed=false; fi
singularity --help > /tmp/singularity_test 2>&1
if [ $? -eq "0" ]; then singularity_installed=true; else singularity_installed=false; fi

if ! $docker_installed && ! $singularity_installed; then
  echo "  ERROR: Neither docker nor singularity seem to be installed."
  echo "    Install either of the two to successfully run metabomodules."
  echo "    Please visit https://www.docker.com or http://singularity.lbl.gov"
  echo "" && echo "ABORTING: metabomodules WAS NOT INSTALLED."
  exit 1
else
  echo "  ...YES"
fi

# store metabomodules code in the home directory
echo "- Copying files..."
mkdir ~/.metabomodules
cp -r ./* ~/.metabomodules
chmod -R 750 ~/.metabomodules
echo "  ...DONE"

# make monet command available: add location where metabomodules is installed to $PATH
case ":$PATH:" in
  *"metabomodules"*)
    # already in path due to previous installation, nothing to do
    echo "" && echo "SUCCESS: installation completed."
    ;;
  *)
    echo ""; echo export PATH=\"'$PATH':$HOME/.metabomodules\" >>  ~/.bash_profile
    echo ""
    echo "SUCCESS: please provide your password to finalize the installation."
    echo "You are being REDIRECTED TO YOUR HOME DIRECTORY"
    echo "INVOKE metabomodules FROM ANY LOCATION"
    username=$(whoami)
    su - $username #re-login and execute the profile script to make metabomodules command available
    ;;
esac
