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

# load parameters from file
. /metabomodules_code/parameters.txt

# run metabomodules
cd /metabomodules_code
#Rscript --vanilla ./Final.R "$filename" "$b" "$c" "$i" "$filter" "$threshold" "$interWeight" "$weighted" "$dir" "$post" "$smallest" "$largest" "$b2" "$c2" "$i2"
python3 metabomodules.py --octave --input ./input.csv --output ./output \
	-N "$number_pseudospectra" \
	--ACP_OffDiagDist "$ACP_OffDiagDist"

# docker generates output files owned by root: make them read/writable
chmod 777 -R /metabomodules_code


