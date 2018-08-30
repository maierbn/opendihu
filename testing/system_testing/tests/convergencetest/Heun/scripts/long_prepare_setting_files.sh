for N in {1..8}
do
  yes | cp -rf ../settings/long_settings_TBR.py ../settings/long_settings_$N.py
  sed -i "s@TBR@$N@g" ../settings/long_settings_$N.py
done

