#!/bin/bash
mkdir -v /home/anthony/Downloads/0036500097
cp -rv {CMakeLists.txt,copy.sh,env,include,octave,readme.txt,src} /home/anthony/Downloads/0036500097/
cd /home/anthony/Downloads
zip -rv 0036500097.zip 0036500097
rm -rf 0036500097
