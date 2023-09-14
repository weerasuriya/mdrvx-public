#!/usr/bin/env bash

# Copy main file
cp ../2_main.R 2_main_rebuild.R
file=2_main_rebuild.R
# List of source files

sourcefiles=($(sed -n -E 's/source\("(.*)".*\)/\1/p' $file))
for i in "${sourcefiles[@]}"
do
	echo "$i"
#	sed -i "/source("$i")/r ../$i" $1
#	echo -e "\n" >> $i
	sed -E -i -e '/source\("'"$i"'".*\)/{
	p
	r '../"${i}"'
	G
	}' $file
done

sed -i -e '/source(/d' $file
