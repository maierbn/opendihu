mkdir -p tmp
MAXP=20
# choose an even divider, as the time integration for the reaction term is split in half
DTSS=("3e-3")
DTSS=("6e-3")
DTSS=("3e-3" "6e-3")
DTSS=("6e-3" "3e-3" "3e-3*2**-1" "3e-3*2**-2" "3e-3*2**-3" "3e-3*2**-4")
MS=(64 32 16 8 4 2) # dt_3D/M for reaction
NS=(64 32 16 8 4 2 1) # dt_3D/N for diffusion
# NS=(1)
# MS=(2)
for TEST in "short"; do
  for DTS in ${DTSS[@]}; do
    for M in ${MS[@]}; do  
      for N in ${NS[@]}; do
        if test $TEST = "short"; then
          ID="short"_$DTS"_"$N"_"$M
          ENDTIME="16*"$DTS
          OUTINTERVAL="1"
          CHECK="0000016"
        elif test $TEST = "long"; then
          ID=$DTS"_"$N"_"$M
          ENDTIME="20.0"
          OUTINTERVAL="64"
          CHECK="0000050"
        elif test $TEST = "medium"; then
          ID="medium_"_$DTS"_"$N"_"$M
          ENDTIME="10.0"
          OUTINTERVAL="64"
          CHECK="0000025"
        fi
        
        if test ! -e "out/fibre_dt_"$ID"_"$CHECK".vtp"; then
          echo "tmp/_settings_"$ID".py"
          sed ' 
            s|dt_0D = 3e-3|dt_0D = '$DTS'/'$M'|; 
            s|dt_1D = 1e-3|dt_1D = '$DTS'/'$N'|; 
            s|dt_3D = 3e-3|dt_3D = '$DTS'|; 
            s|end_time = 20.0|end_time = '$ENDTIME'|;
            s|#nInstances = 1|nInstances = 1|; 
            s|"filename": "out/fibre_"+str(i)|"filename": "out/fibre_dt_'$ID'"|; 
            s|"filename": "out/fibre_r_"+str(i)|"filename": "out/fibre_dt_r_'$ID'"|; 
            s|"outputInterval": int(1./dt_0D\*output_timestep)|"outputInterval": '$M'*'$OUTINTERVAL'|; 
            s|"outputInterval": int(1./dt_1D\*output_timestep)|"outputInterval": '$N'*'$OUTINTERVAL'|; 
            ' ../settings_multiple_fibers.py > tmp/_settings_$ID.py
          while test $(ps -e | grep multiple_fibers | grep -v grep | wc -l) -ge $MAXP; do
            sleep 5
          done
          #mpirun -np 4 -bind-to core   
          ./multiple_fibers tmp/_settings_$ID.py &
          sleep 2
        else
          echo tmp/_settings_$ID.py skipped
        fi
      done
    done
  done
done