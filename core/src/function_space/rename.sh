for file in `ls *basis_on_mesh*`; do
echo "$file"
echo "$(echo $file | sed 's/basis_on_mesh/function_space/g')"
mv $file "$(echo $file | sed 's/basis_on_mesh/function_space/g')"
done
