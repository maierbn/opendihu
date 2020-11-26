# run tensile tests with different number of elements, to be called inside build_release, i.e.
# cd build_release
# ../run.sh

rm result_convergence.csv

#for scenario in compressible_mooney_rivlin compressible_mooney_rivlin_decoupled incompressible_mooney_rivlin nearly_incompressible_mooney_rivlin nearly_incompressible_mooney_rivlin_decoupled linear nearly_incompressible_mooney_rivlin_febio; do
for scenario in nearly_incompressible_mooney_rivlin_febio nearly_incompressible_mooney_rivlin; do
#scenario=compressible_mooney_rivlin 
for nx in `seq 1 1 10`; do

./$scenario ../settings_convergence.py  $scenario $nx

done

done
