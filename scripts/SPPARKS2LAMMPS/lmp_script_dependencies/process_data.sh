
myvar=$(sed -n '1p' sites.txt)

if [ "$myvar" = "ITEM: TIMESTEP" ]; then 
sed -i '1,3d;5d;9d' sites.txt
sed -i "1s/$/ atoms \n&/" sites.txt
sed -i '3s/$/ xlo xhi &/' sites.txt
sed -i '4s/$/ ylo yhi &/' sites.txt
sed -i '5s/$/ zlo zhi &/' sites.txt
sed -i '6i\ \nAtoms \n' sites.txt
sed -i '1i\ ' sites.txt
fi
