# run shear tests with different force, to be called inside build_release, i.e.
# cd build_release
# ../run_force.sh


#rm result.csv

for scenario in compressible_mooney_rivlin compressible_mooney_rivlin_decoupled incompressible_mooney_rivlin nearly_incompressible_mooney_rivlin nearly_incompressible_mooney_rivlin_decoupled linear nearly_incompressible_mooney_rivlin_febio; do
#for scenario in nearly_incompressible_mooney_rivlin_decoupled nearly_incompressible_mooney_rivlin_febio; do
#for scenario in nearly_incompressible_mooney_rivlin_decoupled; do
#scenario=compressible_mooney_rivlin 
for force in `seq 0 2 50`; do

./$scenario ../settings_force.py  $scenario $force

done

done
