n=1 # segment count
seg=2222222222222222222222222111111111111111111112222222222222222222222222 # segment pattern

NA=100 # A-chain count
NB=100 # B-chain count 
#L=300 # 
L=300
L2=$((2*$L))

pck_inp='populate_tmp.inp'

pck_out='IC_tmp.xyz' 

cvfile='N400_Rg_L700.colvars'

sysName="b70_N200_L$L.lt" # double quote is needed for bash 

dataFile="b70_N200_L$L.data"

lmp_input='Template_input.in'

python3 LT_writer.py $n $seg 

python3 writePackmolInput.py $n $NA $NB $L $pck_inp $pck_out

python3 writeSysLT.py $n $NA $NB $L $sysName

packmol < $pck_inp

moltemplate.sh -xyz $pck_out $sysName -nocheck
python3 updateColVar.py $pck_out $cvfile $L $n $NA $NB $seg
python3 updateInput.py $lmp_input $L
python3 fix_datafiles.py $dataFile

done 

