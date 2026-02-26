n=5 # segment count
seg=22222111122222 # segment pattern

NA=100 # A-chain count
NB=100 # B-chain count 
#L=300 # 
L=400
L2=$((2*$L))

pck_inp='populate_tmp.inp'

pck_out='IC_tmp.xyz' 

sysName="b70_N200_L$L.lt" # double quote is needed for bash 

dataFile="b70_N200_L$L.data"

lmp_input='Template_input.in'

python3 LT_writer.py $n $seg 

python3 writePackmolInput.py $n $NA $NB $L $pck_inp $pck_out

python3 writeSysLT.py $n $NA $NB $L $sysName

packmol < $pck_inp

moltemplate.sh -xyz $pck_out $sysName -nocheck
python3 updateInput.py $lmp_input $L
python3 fix_datafiles.py $dataFile

done 

