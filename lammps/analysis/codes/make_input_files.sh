####################
#User defined inputs
pos=$1       #User needs to input the filename of lammpstrj file
nequil=0     #Number of initial steps not included for analysis
stride=1     #Run analysis at each $stride frames
nhistrdf=400 #Number of points in the histogram for radial distribution
ntype=3      #Number of atom types
zdir=3       #Direction normal to the TiO2 surface (x:1,y:2,z:3)
zoffset=0    #Offset number density distribution by a constant

###########################################
#Variables obtained from the lammpstrj file
nat=`head $pos | awk '(NR==4){print $1}'`
nframes=`grep TIMESTEP $pos | wc -l | awk '{print $1}'`
nattype=`tail -n $nat $pos | awk -v"a=$atypemsd" 'BEGIN{c=0} ($2==a){c+=1} END{print c}'`
a=`tail -n $((nat+9)) $pos | head -n 10 | awk '(NR==6){print $2-$1}'`
b=`tail -n $((nat+9)) $pos | head -n 10 | awk '(NR==7){print $2-$1}'`
c=`tail -n $((nat+9)) $pos | head -n 10 | awk '(NR==8){print $2-$1}'`

#######################
#Making the input files
cat << EOF > index_rdf
$nat $nframes $nequil $stride $nhistrdf
$a $b $c
EOF

cat << EOF > index_zdf
$nat $ntype $nframes $nequil $stride $nhistrdf $zdir
$zoffset
EOF
