for file in `find . -type f \( -name "*.cpp" -o -name "*.h" -o -name "*.tpp" \)`; do
  sed -i -e '$a\' $file
  #echo $file
done
