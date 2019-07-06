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
echo "Uninstalling metabomodules"
echo 

# ask for root password (needed to remove /usr/local/bin/monet)
echo "Superuser rights are required."
sudo ls > /dev/null
echo


sudo rm -f /usr/local/bin/metabomodules
rm -rf ~/.metabomodules
echo FINISHED: metabomodules was UNINSTALLED SUCCESSFULLY.
echo
