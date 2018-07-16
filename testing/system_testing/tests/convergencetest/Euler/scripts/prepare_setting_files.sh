for N in {1..8}
do
  yes | cp -rf ../settings/settings_TBR.py ../settings/settings_$N.py
  sed -i "s@TBR@$N@g" ../settings/settings_$N.py
done

