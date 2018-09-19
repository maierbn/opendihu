

# settings
n_rings=43   # this has to be odd for quadratic elements
#n_rings=3
min_z=40
max_z=260

triangulation_type=2  # triangulation_type:  0 = scipy, 1 = triangle, 2 = custom with CoG, 3 = custom with minimized distance
parametric_space_shape=4   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid, 4 = unit circle with optimized positions
n_points_x=10   # has to be even for quadratic elements
#n_points_x=4



input_file=../meshes/biceps_full.stl

#./remove_inside_triangles.py ${input_file} ${input_file}_repaired.stl    # remove non-surface triangles (takes ~1h)

mkdir -p out

cp $input_file out/mesh_00.stl    # cp input file

# cut surface mesh at specified z positions and create rings from it. Write result to `rings_created`, debugging output to out/mesh_01.stl (takes <1min)
# usage: ./create_rings.py <input file> [<n rings> [<min z> <max z>]]
./create_rings.py $input_file $n_rings $min_z $max_z    

# extract existing rings from surface mesh   (this is the alternative to create_rings.py)   Write result to `rings_extracted`
#./extract_rings.py ../biceps.stl     

# rename ring output file to `rings`
mv rings_created rings      
#mv rings_extracted rings

# create a mesh, reads input from `rings`, write output to `mesh`, debugging output to out/mesh_02* - out/mesh_09* (takes ~2min)
# usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x>]]]]"

./create_mesh.py $triangulation_type $parametric_space_shape $n_points_x

