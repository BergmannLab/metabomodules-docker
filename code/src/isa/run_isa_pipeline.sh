input=$1 #e.g. /p-sara1/mirjam/ISApaper/data/colaus1.focus.mr.decov.20180815.csv
o=$2 #e.g. /p-sara1/mirjam/ISApaper/data.out/2018-09-19__120252/ps.isa.colaus1.focus.mr.decov.20180815/isa.colaus1.focus.mr.decov.20180815.pseudospectrum.tsv
python3=$3 #e.g. /usr/local/bin/
base=$4 # location of the current script

#move input to script location 
cp $input $base
cd $base

# run ISA
f=$(basename $input)
"$python3"python3 ./isawrp.py -i $f -o $f \
--inputhasheaders --inputhaslabels --gopseudo \
--seedsparsity 3 --nt \
--thc 1.:1.:7. --thr 1.:1.:7. \
--sgc 0 --sgr 1 --dconv 0.99 \
--dsame 0.50 --nseed 250 

for i in $(seq 11 30)
  do
  "$python3"python3 ./isawrp.py -i $f -o ${f/.csv/.attract$i.csv} \
  --inputhasheaders --inputhaslabels --gopseudo \
  --seedsparsity 3 --nt \
  --thc 1.:1.:7. --thr 1.:1.:7. \
  --sgc 0 --sgr 1 --dconv 0.99 \
  --nseed 500 --nosweep --nopurge
done

"$python3"python3 ./attractor.py ${f/csv/colscore.tsv}


# move output to user-defined locaiton
mv isa.${f/csv/}pseudospectrum.tsv $o
mv isa.${f/csv/}info.tsv ${o/.pseudospectrum.tsv}.info.tsv

#rm temporary files
#rm *.tsv *.csv
