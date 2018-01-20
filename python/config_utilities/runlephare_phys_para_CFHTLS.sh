#! /bin/sh
#Script to run lephare
#usage: runlephare.sh file

LEPHAREDIR=$HOME/lephare_dev
PARA=$HOME/lephare_dev/config/lephare_CFHTLenS.para

LEPHAREWORK=$HOME/lephare_dev/work

export LEPHAREDIR
export LEPHAREWORK

LIBR=MAG_BC03_I09

FILEIN=$1
FILEOUT=${FILEIN%.lepharein}.${LIBR}.lephareout

#Select the correct options depending on 
#the "i" band. i:old one(4th filter), y:new one(5th filter).
if [[ "${FILEIN}" =~ "izrgu" ]]
then
    OPTIONS="-NZ_PRIOR 4,2,4 -MAG_REF 4 -MABS_REF 4 -ZMAX_FILT 4"
fi
if [[ "${FILEIN}" =~ "yzrgu" ]]
then
   OPTIONS="-NZ_PRIOR 5,2,5 -MAG_REF 5 -MABS_REF 5 -ZMAX_FILT 5"
fi

echo "Computing physical parameters"

rm -f $FILEOUT

$LEPHAREDIR/source/zphota -c $PARA $OPTIONS -CAT_IN $FILEIN -CAT_OUT $FILEOUT -ZPHOTLIB $LIBR,MAG_STAR,MAG_QSO -ZFIX yes -ZMAX_MAGLIM  24.0


