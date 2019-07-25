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
python3 metabomodules.py --octave --input ./input.csv --output ./output \
	-N "$number_pseudospectra" \
	--number_permutations "$number_permutations" \
	--OffDiagDist "$ACP_OffDiagDist" \
	--remNeigbPairsFlag "$ACP_remNeigbPairsFlag" \
	--dsame "$ISA_dsame" \
	--dconv "$ISA_dconv" \
	--nseed "$ISA_nseed" \
	--seedsparsity "$ISA_seedsparsity" \
	--maxiter "$ISA_maxiter" \
	--sgc "$ISA_sgc" \
	--sgr "$ISA_sgr" \
	--thc "$ISA_thc" \
	--thr "$ISA_thr" \
	--norm "$ISA_norm" \
	--nt "$ISA_nt" \
	--inputhasheaders "$ISA_inputhasheaders" \
	--inputhaslabels "$ISA_inputhaslabels" \
	--nopurge "$ISA_nopurge" \
	--quiet "$ISA_quiet" \
	--nosweep "$ISA_nosweep" \
	--onefile "$ISA_onefile" \
	--gopseudo "$ISA_gopseudo" \
	--FILTER_z_score_th "$FILTER_z_score_th" \
	--FILTER_adj_score_th "$FILTER_adj_score_th" \
	--FILTER_redo "$FILTER_redo"


#--ISA_seedmatrix "$ISA_seedmatrix" \


# docker generates output files owned by root: make them read/writable
chmod 777 -R /metabomodules_code
