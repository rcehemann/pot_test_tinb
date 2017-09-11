#!/bin/bash

source ~/.bashrc

# tests binary spline MEAM potentials
ps2png=~/utils/ps2pnghr
fixbb=~/utils/fixbb
BMFit="/n/jww-1/ehemann.2/testingScripts/BMFit.py"
cspline=/n/jww-1/ehemann.2/testingScripts/cspline.py
genstruct="perl /n/jww-1/ehemann.2/TestingScripts/binary/randbin.pl" 	# perl script for generating random alloy structures
lammps="mpirun -np 4 /n/jww-1/ehemann.2/lammps_test_17/src/lmp_mpicc"
meamz="mpirun -np 4 /n/jww-1/ehemann.2/Meamzilla/build/meamz"
genstruct="/n/jww-1/ehemann.2/testingScripts/genstruct.py"
fit="/n/jww-1/ehemann.2/testingScripts/linFit.py"
ecgen="/n/jww-1/ehemann.2/testingScripts/LEC.py"
gammaTilt="/n/jww-1/ehemann.2/testingScripts/gammaTilt.py"
betaApp="/n/jww-1/ehemann.2/testingScripts/betaAlphappAlloy.py"
betaOmegaAlloy="/n/jww-1/ehemann.2/testingScripts/betaOmegaAlloy.py"
sqsgen="python /n/jww-1/ehemann.2/testingScripts/sqsgen.py"

readhcpphon="/n/jww-1/ehemann.2/testingScripts/read_hcp_phonons.py"
readbccphon="/n/jww-1/ehemann.2/testingScripts/read_bcc_phonons.py"

#################################
## DEFINE COLORS FOR OUTPUT #####
grn='\e[0;32m'			#
non='\e[0m'			#
blu='\e[0;36m'			#
yel='\e[0;33m'			#
dbl='\e[1;34m'			#
red='\e[1;31m'			#
prp='\e[1;35m'			#
cya='\e[1;36m'			#
#################################

#################################
## TEST FLAGS			#
SS_FLAG="TRUE"
ALPHA_PRIME_FLAG="TRUE"
phonon_flag="TRUE"
interstitial_flag="TRUE"
vacancy_flag="TRUE"
ideal_flag="TRUE"
SQS_FLAG="TRUE"
gamma_flag="TRUE"
gamma_lammps_flag="TRUE"
twin_flag="TRUE"
transition_flag="TRUE"
thermo_flag="FALSE"
betaOmegaLammpsFlag="TRUE"
hcp_gamma_flag="TRUE"
lammps_evsv_flag="TRUE"
lammps_elcon_flag="TRUE"
lammps_elcon_g1="TRUE"
lammps_elcon_APP="TRUE"
lammps_SQS_flag="TRUE"
lammps_AP_flag="TRUE"
#################################

#################################
## TEST PARAMETERS		#
PRESMAX=100
PRESSEQ=$(seq 0 10 $PRESMAX)
PRESSLAT=$(seq 0 5 100)
stpct=50	# strain % for E-V
nevpt=1000	# number of points for E-V
ecstr=1		# % strain for EC calcs
MDIN=20		# inner point for Take[]
MDOU=80		# outer poutn for Take[] (BM fit)
gampts=100
BETA_RELAX_ETOL="1e-6"
META_RELAX_ETOL="1e-4"
PCONV=160.2176487
declare -A masses=( ["V"]='50.94' ["Nb"]='92.91' ["Mo"]='95.96' ["Ta"]='180.95' ["W"]='183.84' ["Ti"]='47.867' ) # masses of your favorite elements
declare -a pureTiOMG=( ["0"]='4.58028' ["25"]='4.347233' ["50"]='4.17987' ["75"]='4.0571646' ["100"]='3.96762' )
#################################

elem1=$1
elem2=$2
dftdat="/n/jww-1/ehemann.2/testingScripts/DFT_DATA"
typ="GMEAM"
#typ="LREP"
#typ="EAM"
#pairstyle="meam/alloy/spline"
pairstyle="gmeam/spline"
#pairstyle="lrep"
#pairstyle="eam/alloy"
mass1=${masses[$elem1]}
mass2=${masses[$elem2]}

declare -A idx=( ["Ti"]='1' ["Nb"]='2' )

###########################
#	STORE DFT VALUES
#

elem10=$elem1
elem20=$elem2

if [ "$elem1" == "Nb" ]; then
	elem1="Ti"
	elem2="Nb"
fi

D03lat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/D03/values.dat"`
G1lat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/G1/values.dat"`
SQS7525lat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/SQS7525/values.dat"`

SQS5050lat=3.25
SQS2575lat=3.30

D03ep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/D03/values.dat"`
G1ep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/G1/values.dat"`
SQS7525ep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/SQS7525/values.dat"`

D03bulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/D03/values.dat"`
G1bulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/G1/values.dat"`
SQS7525bulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/SQS7525/values.dat"`

APPlat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
APlat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/AP/values.dat"`
omglat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`

L60Ti3Nblat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/L60Ti3Nb/values.dat"`
L60TiNb3lat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/L60TiNb3/values.dat"`

L60Ti3Nbep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/L60Ti3Nb/values.dat"`
L60TiNb3ep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/L60TiNb3/values.dat"`

L60Ti3Nbbulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/L60Ti3Nb/values.dat"`
L60TiNb3bulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/L60TiNb3/values.dat"`

fccTiep=`awk 'NR==1{print $2}' "$dftdat/Ti/fcc/values.dat"`
fccTilat=`awk 'NR==3{print $2}' "$dftdat/Ti/fcc/values.dat"`
fccTibulkd=`awk 'NR==4{print $2}' "$dftdat/Ti/fcc/values.dat"`

A15Tiep=`awk 'NR==1{print $2}' "$dftdat/Ti/A15/values.dat"`
A15Tilat=`awk 'NR==3{print $2}' "$dftdat/Ti/A15/values.dat"`
A15Tibulkd=`awk 'NR==4{print $2}' "$dftdat/Ti/A15/values.dat"`

A15Ti3Nbep=`awk 'NR==1{print $2}' "$dftdat/Ti-Nb/A15Ti3Nb/values.dat"`
A15Ti3Nblat=`awk 'NR==3{print $2}' "$dftdat/Ti-Nb/A15Ti3Nb/values.dat"`
A15Ti3Nbbulkd=`awk 'NR==4{print $2}' "$dftdat/Ti-Nb/A15Ti3Nb/values.dat"`

A15TiNb3ep=`awk 'NR==1{print $2}' "$dftdat/Ti-Nb/A15TiNb3/values.dat"`
A15TiNb3lat=`awk 'NR==3{print $2}' "$dftdat/Ti-Nb/A15TiNb3/values.dat"`
A15TiNb3bulkd=`awk 'NR==4{print $2}' "$dftdat/Ti-Nb/A15TiNb3/values.dat"`
A15TiNb3C11d=`awk 'NR==5{print $2}' "$dftdat/Ti-Nb/A15TiNb3/values.dat"`
A15TiNb3C12d=`awk 'NR==6{print $2}' "$dftdat/Ti-Nb/A15TiNb3/values.dat"`
A15TiNb3C44d=`awk 'NR==7{print $2}' "$dftdat/Ti-Nb/A15TiNb3/values.dat"`

omgTi2Nbep=`awk 'NR==1{print $2}' "$dftdat/Ti-Nb/omgTi2Nb/values.dat"`
omgTi2Nblat=`awk 'NR==3{print $2}' "$dftdat/Ti-Nb/omgTi2Nb/values.dat"`
omgTi2Nbcoa=`awk 'NR==4{print $2}' "$dftdat/Ti-Nb/omgTi2Nb/values.dat"`

omgTiNb2ep=`awk 'NR==1{print $2}' "$dftdat/Ti-Nb/omgTiNb2/values.dat"`
omgTiNb2lat=`awk 'NR==3{print $2}' "$dftdat/Ti-Nb/omgTiNb2/values.dat"`
omgTiNb2coa=`awk 'NR==4{print $2}' "$dftdat/Ti-Nb/omgTiNb2/values.dat"`

fccNbep=`awk 'NR==1{print $2}' "$dftdat/Nb/fcc/values.dat"`
fccNblat=`awk 'NR==3{print $2}' "$dftdat/Nb/fcc/values.dat"`
fccNbbulkd=`awk 'NR==4{print $2}' "$dftdat/Nb/fcc/values.dat"`

A15Nbep=`awk 'NR==1{print $2}' "$dftdat/Nb/A15/values.dat"`
A15Nblat=`awk 'NR==3{print $2}' "$dftdat/Nb/A15/values.dat"`
A15Nbbulkd=`awk 'NR==4{print $2}' "$dftdat/Nb/A15/values.dat"`

hcpNbep=`awk 'NR==1{print $2}' "$dftdat/Nb/hcp/values.dat"`
hcpNblat=`awk 'NR==3{print $2}' "$dftdat/Nb/hcp/values.dat"`
hcpNbcoa=`awk 'NR==4{print $2}' "$dftdat/Nb/hcp/values.dat"`
hcpNbbulkd=`awk 'NR==5{print $2}' "$dftdat/Nb/hcp/values.dat"`

L10lat=4.489
L10coa=0.80

D019Ti3Nblat=5.900
D019Ti3Nbcoa=0.816

APPep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
APep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/AP/values.dat"`
omgep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`

APPbulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
APbulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/AP/values.dat"`
omgbulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`

APPboa=`awk 'NR==5{print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
APPcoa=`awk 'NR==6{print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
omgcoa=`awk 'NR==5{print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`
APcoa=`awk 'NR==5{print $2}' "$dftdat/$elem1-$elem2/AP/values.dat"`

B2lat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/TiNb/B2/values.dat"`
B2bulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/TiNb/B2/values.dat"`
B2ep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/TiNb/B2/values.dat"`

bcc110lat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/TiNb/bcc110/values.dat"`
bcc110bulkd=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/TiNb/bcc110/values.dat"`
bcc110ep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/TiNb/bcc110/values.dat"`

A3lat=`awk 'NR==3{print $2}' "$dftdat/$elem1-$elem2/TiNb/A3/values.dat"`
A3coa=`awk 'NR==4{print $2}' "$dftdat/$elem1-$elem2/TiNb/A3/values.dat"`
A3bulkd=`awk 'NR==5{print $2}' "$dftdat/$elem1-$elem2/TiNb/A3/values.dat"`
A3ep=`awk 'NR==1{print $2}' "$dftdat/$elem1-$elem2/TiNb/A3/values.dat"`

C11B2d=`awk 'NR==5{print $2}' "$dftdat/$elem1-$elem2/TiNb/B2/values.dat"`
C12B2d=`awk 'NR==6{print $2}' "$dftdat/$elem1-$elem2/TiNb/B2/values.dat"`
C44B2d=`awk 'NR==7{print $2}' "$dftdat/$elem1-$elem2/TiNb/B2/values.dat"`

bccTilat=`awk 'NR==3{print $2}' "$dftdat/$elem1/bcc/values.dat"`
bccNblat=`awk 'NR==3{print $2}' "$dftdat/$elem2/bcc/values.dat"`
hcpTilat=`awk 'NR==3{print $2}' "$dftdat/$elem1/hcp/values.dat"`
hcpTicoa=`awk 'NR==4{print $2}' "$dftdat/$elem1/hcp/values.dat"`
omgTilat=`awk 'NR==3{print $2}' "$dftdat/$elem1/omega/values.dat"`
omgTicoa=`awk 'NR==4{print $2}' "$dftdat/$elem1/omega/values.dat"`

omgNblat=`awk 'NR==3{print $2}' "$dftdat/$elem2/omega/values.dat"`
omgNbcoa=`awk 'NR==4{print $2}' "$dftdat/$elem2/omega/values.dat"`

bccTiep=`awk 'NR==1{print $2}' "$dftdat/$elem1/bcc/values.dat"`
bccNbep=`awk 'NR==1{print $2}' "$dftdat/$elem2/bcc/values.dat"`
hcpTiep=`awk 'NR==1{print $2}' "$dftdat/$elem1/hcp/values.dat"`
omgTiep=`awk 'NR==1{print $2}' "$dftdat/$elem1/omega/values.dat"`

omgNbep=`awk 'NR==1{print $2}' "$dftdat/$elem2/omega/values.dat"`

A15nblat=`awk 'NR==3{print $2}' "$dftdat/$elem2/A15/values.dat"`
A15nbeqp=`awk 'NR==1{print $2}' "$dftdat/$elem2/A15/values.dat"`

C11D03d=`awk 'NR==5{print $2}' "$dftdat/$elem1-$elem2/D03/values.dat"`
C12D03d=`awk 'NR==6{print $2}' "$dftdat/$elem1-$elem2/D03/values.dat"`
C44D03d=`awk 'NR==7{print $2}' "$dftdat/$elem1-$elem2/D03/values.dat"`

C11G1d=`awk 'NR==5{print $2}' "$dftdat/$elem1-$elem2/G1/values.dat"`
C12G1d=`awk 'NR==6{print $2}' "$dftdat/$elem1-$elem2/G1/values.dat"`
C44G1d=`awk 'NR==7{print $2}' "$dftdat/$elem1-$elem2/G1/values.dat"`

C11APPd=`awk 'NR==7  {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C12APPd=`awk 'NR==8  {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C13APPd=`awk 'NR==9  {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C22APPd=`awk 'NR==10 {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C23APPd=`awk 'NR==11 {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C33APPd=`awk 'NR==12 {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C44APPd=`awk 'NR==13 {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C55APPd=`awk 'NR==14 {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`
C66APPd=`awk 'NR==15 {print $2}' "$dftdat/$elem1-$elem2/APP/values.dat"`

C11omgd=`awk 'NR==6  {print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`
C12omgd=`awk 'NR==7  {print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`
C13omgd=`awk 'NR==8  {print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`
C33omgd=`awk 'NR==9  {print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`
C44omgd=`awk 'NR==10 {print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`
C66omgd=`awk 'NR==11 {print $2}' "$dftdat/$elem1-$elem2/omega/values.dat"`

C11BCCtid=`awk 'NR==5{print $2}' "$dftdat/$elem1/bcc/values.dat"`
C12BCCtid=`awk 'NR==6{print $2}' "$dftdat/$elem1/bcc/values.dat"`
C44BCCtid=`awk 'NR==7{print $2}' "$dftdat/$elem1/bcc/values.dat"`

C11BCCnbd=`awk 'NR==5{print $2}' "$dftdat/$elem2/bcc/values.dat"`
C12BCCnbd=`awk 'NR==6{print $2}' "$dftdat/$elem2/bcc/values.dat"`
C44BCCnbd=`awk 'NR==7{print $2}' "$dftdat/$elem2/bcc/values.dat"`

C11HCPtid=`awk 'NR==6{print $2}' "$dftdat/$elem1/hcp/values.dat"`
C12HCPtid=`awk 'NR==7{print $2}' "$dftdat/$elem1/hcp/values.dat"`
C13HCPtid=`awk 'NR==8{print $2}' "$dftdat/$elem1/hcp/values.dat"`
C33HCPtid=`awk 'NR==9{print $2}' "$dftdat/$elem1/hcp/values.dat"`
C44HCPtid=`awk 'NR==10{print $2}' "$dftdat/$elem1/hcp/values.dat"`

C11OMGtid=`awk 'NR==6{print $2}' "$dftdat/$elem1/omega/values.dat"`
C12OMGtid=`awk 'NR==7{print $2}' "$dftdat/$elem1/omega/values.dat"`
C13OMGtid=`awk 'NR==8{print $2}' "$dftdat/$elem1/omega/values.dat"`
C33OMGtid=`awk 'NR==9{print $2}' "$dftdat/$elem1/omega/values.dat"`
C44OMGtid=`awk 'NR==10{print $2}' "$dftdat/$elem1/omega/values.dat"`

hcpPrisEasyd=`grep '^0.50' $dftdat/$elem1/hcp/gamma_eh.dat | awk '{print $2}'`
hcpPrisHardd=`grep '^0.50' $dftdat/$elem1/hcp/gamma_eh.dat | awk '{print $3}'`

elem1=$elem10
elem2=$elem20

#pureTiOMG=`awk 'NR==3 {print $2}' "$dftdat/$elem1/omega/values.dat"`

#
#
########################

# decide how many pots to test
if [ "$4" == "" ]; then
	if [ "$3" == "test" ]; then
		potz=1
	else
		npots=10
		potz=$(seq 1 $npots)
	fi
else
	potz=$4
fi

# main loop
for setnum in $3; do
for potnum in $potz; do

echo -e "${blu}set $setnum pot $potnum ${non}"

# prepare HTML directories
mkdir -p ~/public_html/POTS/"$elem1"-"$elem2"/$typ/"set.$setnum"/
mkdir -p ~/public_html/POTS/"$elem1"-"$elem2"/$typ/"set.$setnum"/"pot.$potnum"/

# prepare test directories
dir="/n/jww-1/ehemann.2/pots/$typ/$elem1-$elem2/set.$setnum/pot.$potnum"
wd=`pwd`
rm -f *.pckl
cd $dir
	mkdir -p tests
	cd tests
		mkdir -p HCPti
			rm -f ./HCPti/*
		mkdir -p BCCti
			rm -f ./BCCti/*
		mkdir -p A15ti
			rm -f ./A15ti/*
		mkdir -p FCCti
			rm -f ./FCCti/*
		mkdir -p OMGti
			rm -f ./OMGti/*
		mkdir -p BCCnb
			rm -f ./BCCnb/*
		mkdir -p FCCnb
			rm -f ./FCCnb/*
		mkdir -p HCPnb
			rm -f ./HCPnb/*
		mkdir -p OMGnb
			rm -f ./OMGnb/*
		mkdir -p A15nb
			rm -f ./A15nb/*
		mkdir -p B2
			rm -f ./B2/*
		mkdir -p bcc110
			rm -f ./bcc110/*
		mkdir -p A3
			rm -f ./A3/*
		mkdir -p D03
			rm -f ./D03/*
		mkdir -p G1
			rm -f ./G1/*
		mkdir -p SQS7525
			rm -f ./SQS7525/*
		mkdir -p SQS5050
			rm -f ./SQS5050/*
		mkdir -p SQS2575
			rm -f ./SQS2575/*
		mkdir -p L10
			rm -f ./L10/*
		mkdir -p L60Ti3Nb
			rm -f ./L60Ti3Nb/*
		mkdir -p L60TiNb3
			rm -f ./L60TiNb3/*
		mkdir -p D019Ti3Nb
			rm -f ./D019Ti3Nb/*
		mkdir -p A15Ti3Nb
			rm -f ./A15Ti3Nb/*
		mkdir -p A15TiNb3
			rm -f ./A15TiNb3/*
		mkdir -p omgTi2Nb
			rm -f ./omgTi2Nb/*
		mkdir -p omgTiNb2
			rm -f ./omgTiNb2/*
		mkdir -p APP
			rm -f ./APP/*
		mkdir -p AP
			rm -f ./AP/*
		mkdir -p omega
			rm -f ./omega/*
		mkdir -p transitions
			rm -f ./transitions/*
		mkdir -p thermo
			rm -f ./thermo/*
		mkdir -p SS
			rm -f ./SS/*
	cd ../
cd $wd


# paths to potential files
potPath="$dir/lammps.pt"
mmzpot="$dir/end"

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#		T E S T I N G    P O R T I O N
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#################################################################
#	single element tests
#################################################################
echo -e "${grn}TESTING SINGLE ELEMENTS... ${non}"

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HCPti >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t HCP-Ti... ${non}"

cat > $dir/tests/HCPti/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid
variable boxa equal $hcpTilat
variable boxb equal $hcpTilat*sqrt(3)/2
variable boxc equal $hcpTilat*$hcpTicoa
variable boxxy equal $hcpTilat*(-0.5)

lattice		custom $hcpTilat a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $hcpTicoa &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
		
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0 scalexy yes scalexz yes scaleyz yes
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/HCPti/eq.in > $dir/tests/HCPti/eq.out

HCPtieqp=`tail -1 $dir/tests/HCPti/eq.out | awk '{print $1}'`
HCPticoap=`tail -1 $dir/tests/HCPti/eq.out | awk '{print $2}'`
HCPtivop=`tail -1 $dir/tests/HCPti/eq.out | awk '{print $3}'`
HCPtipote=`tail -1 $dir/tests/HCPti/eq.out | awk '{print $4}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for HCPti -----------------------------------------

cat > $dir/tests/HCPti/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR HCPti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $HCPtieqp
variable boxb equal $HCPtieqp*sqrt(3)/2
variable boxc equal $HCPtieqp*$HCPticoap
variable boxxy equal $HCPtieqp*(-0.5)

lattice		custom $HCPtieqp a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $HCPticoap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
		
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

dump D all xyz 1 $dir/tests/HCPti/evsv.xyz
fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/HCPti/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/HCPti/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd  equal  \${boxa}*(1-v_dmax/2)
variable xf  equal  \${boxa}*(1+v_dmax/2)
variable yd  equal  \${boxb}*(1-v_dmax/2)
variable yf  equal  \${boxb}*(1+v_dmax/2)
variable zd  equal  \${boxc}*(1-v_dmax/2)
variable zf  equal  \${boxc}*(1+v_dmax/2)
variable xyd equal -0.5*\${xd} 
variable xyf equal -0.5*\${xf} 
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/HCPti/evsv.in loop
!!
$lammps < $dir/tests/HCPti/evsv.in > $dir/tests/HCPti/evsv.out

sed -i '1d' $dir/tests/HCPti/evsv.dat
line=`python $BMFit $dir/tests/HCPti/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/HCPti/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	HCPtibulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	HCPtibulkp=0
fi

sed -i '1d' $dir/tests/HCPti/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/HCPti/pvsv.dat | sort -k2,2 -g | head -1`
	HCPtiPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${HCPtiPLAT[0]} $HCPtieqp"
		#HCPtieqp=`echo $line | awk '{print $3}'`
		#HCPtivop=`echo $line | awk '{print $1}'`
		#HCPtipote=`echo $line | awk '{print 1000*$4}'`
	fi
done
HCPtiPLAT["0"]=$HCPtieqp
awk -v Voo=$HCPtivop '{print $1/Voo, $2, $3}' $dir/tests/HCPti/pvsv.dat > tmp; mv tmp $dir/tests/HCPti/pvsv.dat

# ---------------------- pressure dependence of lattice constants for HCPti---------
echo -e "${prp}Pressure dependence of lattice constants...${non}"

rm -f $dir/tests/HCPti/abcvp.dat; echo "#a, c/a, gamma" > $dir/tests/HCPti/abcvp.dat
for PRESS in $PRESSLAT; do
PRES=`echo "10000*$PRESS" | bc -l`
cat > $dir/tests/HCPti/abcvp.lin << __
## LAMMPS script
units metal
boundary p p p
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable a1	equal	1.00
variable a2	equal	sqrt(3)
variable a3	equal	$HCPticoap
variable boxa equal $HCPtieqp*\${a1}
variable boxb equal $HCPtieqp*\${a2}
variable boxc equal $HCPtieqp*\${a3}

lattice custom $HCPtieqp &
	a1 \${a1}	0.0	0.0 &
	a2 0.0		\${a2}	0.0 &
	a3 0.0		0.0	\${a3} &
	basis 0.0	0.0	0.0 basis 0.5	0.5	0.0 &
	basis 0.0	0.33333 0.5 basis 0.5   0.83333 0.5
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#region frz1 sphere 0.0 0.0 0.0 0.1 units box
#region frz2 sphere 0.5 0.5 0.0 0.1 units box
#region frz union 2 frz1 frz2
#group frz region frz
#
#fix FREEZE freeze frz

thermo 1000
thermo_style custom step temp press pe ke etotal lx ly lz
timestep 0.001

variable a equal lx
variable coa equal lz/lx
variable Gamma equal 180-2*180*atan(ly/lx)/3.141592654

fix P all press/berendsen aniso $PRES $PRES 1000.0 couple xy
run 10000

dump D all custom 1 $dir/tests/HCPti/abcvp.coords xs ys zs

print "$PRESS \${a} \${coa} \${Gamma}"
__

$lammps < $dir/tests/HCPti/abcvp.lin > $dir/tests/HCPti/abcvp.out
tail -1 $dir/tests/HCPti/abcvp.out >> $dir/tests/HCPti/abcvp.dat

done

#---------------------------- stacking faults for hcp ---------------------------------------
if [ "$hcp_gamma_flag" == "TRUE" ]; then
echo "Stacking Faults..."
if [ "$gamma_lammps_flag" == "TRUE" ]; then
for PRESS in 0; do
for jj in $(seq 1 2); do
case $jj in
	1) XLAT="11-20"; ZLAT="01-10e"; FILE=$dir/tests/HCPti/"$PRESS"SF_011b0e_112b0.dat; RILE=$dir/tests/HCPti/"$PRESS"SF_011b0e_112b0_relaxed.dat;;		# prismatic easy
	2) XLAT="11-20"; ZLAT="01-10h"; FILE=$dir/tests/HCPti/"$PRESS"SF_011b0h_112b0.dat; RILE=$dir/tests/HCPti/"$PRESS"SF_011b0h_112b0_relaxed.dat;;		# prismatic hard
esac
printf "$jj "

python $gammaTilt HCP ${HCPtiPLAT["$PRESS"]} $ZLAT $XLAT 0.0 lmp $HCPticoap > $dir/tests/HCPti/read_gamma.dat
declare -a tilts=( `python $gammaTilt HCP ${HCPtiPLAT["$PRESS"]} $ZLAT $XLAT 0.5 lmp $HCPticoap | awk 'NR==7{print 2*$1, 2*$2, 2*$3}'` )
area=`awk 'NR==1{print $7}' $dir/tests/HCPti/read_gamma.dat`
#echo ${tilts[@]}
for i in $(seq 0 $gampts); do
deli=`echo "$i/$gampts" | bc -l`
python $gammaTilt HCP ${HCPtiPLAT["$PRESS"]} $ZLAT $XLAT $deli lmp $HCPticoap > $dir/tests/HCPti/read_gamma.dat
sed -i 's/1 atom types/2 atom types/g' $dir/tests/HCPti/read_gamma.dat
cat > lammps_in << __
units metal
boundary p p p
atom_style atomic
box tilt large
read_data $dir/tests/HCPti/read_gamma.dat
mass 1 $mass1
mass 2 $mass2
pair_style $pairstyle
pair_coeff * * $potPath $elem1 $elem2
thermo 1
thermo_style custom pe ke etotal
timestep 0.01
variable EATOM equal pe
run 0
print "EINIT: \${EATOM}"
fix frz all setforce 0.0 NULL 0.0
min_style fire
minimize 1e-8 1e-8 5000 50000
run 0
print "EATOM: \${EATOM}"
__
$lammps < $dir/tests/HCPti/gamma_in > $dir/tests/HCPti/gamma_out
if [ "$i" == "0" ]; then e0=`grep 'EATOM:' $dir/tests/HCPti/gamma_out | awk '{print $2}'`; mevpaa=0;
			eu0=`grep 'EINIT:' $dir/tests/HCPti/gamma_out | awk '{print $2}'`; mevini=0;
else
mevpaa=`grep 'EATOM:' $dir/tests/HCPti/gamma_out | awk -v a=$area -v e0=$e0 '{print 1000*($2-e0)/a}'`
mevini=`grep 'EINIT:' $dir/tests/HCPti/gamma_out | awk -v a=$area -v e0=$eu0 '{print 1000*($2-e0)/a}'`
fi
echo $deli $mevpaa >> $RILE
echo $deli $mevini >> $FILE

done
done
done

echo ""
else

cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 1

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/HCPti/gamma.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

for PRESS in $(seq 0 25 0); do
P=`echo "10000*$PRESS" | bc -l`

for jj in $(seq 1 2); do
case $jj in
	1) XLAT="11-20"; ZLAT="01-10e"; FILE=$dir/tests/HCPti/"$PRESS"SF_011b0e_112b0.dat;;		# prismatic easy
	2) XLAT="11-20"; ZLAT="01-10h"; FILE=$dir/tests/HCPti/"$PRESS"SF_011b0h_112b0.dat;;		# prismatic hard
esac
printf "$jj "

echo "# DEL, GAMMA" > $FILE
rm -f $dir/tests/HCPti/gamma.conf; touch $dir/tests/HCPti/gamma.conf
for d in $(seq 0 $gampts); do

	DEL=`echo "scale=6; $d/$gampts" | bc -l`
	python $gammaTilt HCP ${HCPtiPLAT["$PRESS"]} $ZLAT $XLAT $DEL conf $HCPticoap >> $dir/tests/HCPti/gamma.conf

done

$meamz -p $dir/tests/meamz_params > $dir/tests/HCPti/meamzilla.out

e0=`awk 'NR==3{print $6}' data.energy`
area=`awk 'NR==2{print $3}' $dir/tests/HCPti/gamma.conf`
tail -n +3 data.energy | awk -v npts=$gampts -v e0=$e0 -v a=$area '{print $1/npts, 1000*$4*($6-e0)/a}' >> $FILE

rm data.*

done
done

echo ""
fi

hcpPrisEasyp=`awk 'NR>1{print ($1-0.5)**2, $2}' $dir/tests/HCPti/0SF_011b0e_112b0.dat | sort -gk1 | awk 'NR==1{print $2}'`
hcpPrisHardp=`awk 'NR>1{print ($1-0.5)**2, $2}' $dir/tests/HCPti/0SF_011b0h_112b0.dat | sort -gk1 | awk 'NR==1{print $2}'`

fi

#---------------------------- elastic constants for HCP -----------------------------------
echo -e "${prp}Elastic constants:${non}"

echo "# pressure, c11, c12, c13, c33, c44" > $dir/tests/HCPti/C_VS_P.dat
for PRESS in 0 10; do
echo "$PRESS GPa..."

for jj in $(seq 1 7); do
printf "\t $jj "

if [ $lammps_elcon_flag == "TRUE" ]; then
P=`echo $PRESS*10000 | bc -l`

# NOTE: THIS NEEDS TO BE ADJUSTED IF SUPERCELLS LARGER THAN 1X1X1 ARE TO BE USED
#case $jj in
#	1) XD="1+v_d"; YD="sqrt(3)*(1+v_d)/2"; ZD="1+v_d"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	2) XD="1+v_d"; YD="sqrt(3)*(1-v_d)/2"; ZD="1/(1-v_d2)"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	3) XD="1/(1-v_d2)"; YD="sqrt(3)*(1+v_d)/2"; ZD="1-v_d"; XYD="-0.5/(1-v_d2)"; XZD="0"; YZD="0" ;;
#	4) XD="1-v_d"; YD="sqrt(3)/(2-2*v_d2)"; ZD="1+v_d"; XYD="-0.5*(1-v_d)"; XZD="0"; YZD="0" ;;
#	5) XD="1/(1-(v_d2)/4)"; YD="sqrt(3)/2"; ZD="1"; XYD="-0.5/(1-(v_d2)/4)"; XZD="0"; YZD="v_d" ;;
#	6) XD="1"; YD="sqrt(3)/((1-(v_d2)/4)*2)"; ZD="1"; XYD="-0.5"; XZD="v_d"; YZD="0" ;;
#	7) XD="1"; YD="sqrt(3)/2"; ZD="1/(1-(v_d2)/4)"; XYD="-0.5+v_d"; XZD="0"; YZD="0" ;;
#esac
case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac
cat > $dir/tests/HCPti/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms
variable boxa equal ${HCPtiPLAT["$PRESS"]}
variable boxb equal ${HCPtiPLAT["$PRESS"]}*sqrt(3)/2
variable boxc equal ${HCPtiPLAT["$PRESS"]}*$HCPticoap
variable boxxy equal ${HCPtiPLAT["$PRESS"]}*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom ${HCPtiPLAT["$PRESS"]} a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $HCPticoap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
create_box 2 mybox
create_atoms ${idx["Ti"]} box
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

dump D all xyz 1 elcon$jj.xyz

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

# minimize
#fix 1 all box/relax aniso $P fixedpoint 0 0 0
#min_style	cg
#minimize	0.0 1e-10 100 1000
#
#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz xy xz yz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/HCPti/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)
variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/HCPti/elcon.lin loop 
!!

$lammps < $dir/tests/HCPti/elcon.lin > $dir/tests/HCPti/elcon$jj.out
sed -i '1d' $dir/tests/HCPti/stresses$jj.dat

else	# compute stress-strain curves with meamzilla
cat > $dir/tests/meamz_params <<@@
ngroups 1

optstyle powell
num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/HCPti/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
rm -f $dir/tests/HCPti/elcon$jj.conf
strain=`echo "$ecstr/100" | bc -l`
for i in $(seq -$DELPT $DELPT); do

	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del HCP ${HCPtiPLAT["$PRESS"]} conf $HCPticoap >> $dir/tests/HCPti/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/HCPti/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=pc*$5; getline	
					syy=pc*$5; getline	
					szz=pc*$5; getline	
					sxy=pc*$5; getline	
					syz=pc*$5; getline	
					szx=pc*$5;

					print del, sxx, syy, szz, syz, szx, sxy

				     }' data.stress >> $dir/tests/HCPti/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""
rm -f $dir/tests/HCPti/EC_fits.dat; touch $dir/tests/HCPti/EC_fits.dat


# first row fits
awk '{print $1, $2}' $dir/tests/HCPti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/HCPti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/HCPti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/HCPti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/HCPti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/HCPti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/HCPti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/HCPti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/HCPti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/HCPti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/HCPti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/HCPti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/HCPti/stresses5.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/HCPti/stresses6.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/HCPti/stresses7.dat > tmp; python $fit tmp >> $dir/tests/HCPti/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/HCPti/EC_fits.dat` )
HCPtiC11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
HCPtiC12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
HCPtiC13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
HCPtiC22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
HCPtiC23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
HCPtiC33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
HCPtiC44i=${coup[12]}
HCPtiC55i=${coup[13]}
HCPtiC66i=${coup[14]}

HCPtiC11p=`python -c "print ($HCPtiC11i + $HCPtiC22i)/2"`
HCPtiC13p=`python -c "print ($HCPtiC13i + $HCPtiC23i)/2"`
HCPtiC44p=`python -c "print ($HCPtiC44i + $HCPtiC55i)/2"`

echo "$PRESS $HCPtiC11p $HCPtiC12i $HCPtiC13p $HCPtiC33i $HCPtiC44p" >> $dir/tests/HCPti/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	HCPtiC11=`printf '%3.f' $HCPtiC11p`
	HCPtiC12=`printf '%3.f' $HCPtiC12i`
	HCPtiC13=`printf '%3.f' $HCPtiC13p`
	HCPtiC33=`printf '%3.f' $HCPtiC33i`
	HCPtiC44=`printf '%3.f' $HCPtiC44p`
	
	C11HCPtid=`printf '%3.f' $C11HCPtid`
	C12HCPtid=`printf '%3.f' $C12HCPtid`
	C13HCPtid=`printf '%3.f' $C13HCPtid`
	C33HCPtid=`printf '%3.f' $C33HCPtid`
	C44HCPtid=`printf '%3.f' $C44HCPtid`
fi
done

## PHONONS FOR HCP TI
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/HCPti/OPT.POSCAR << -
HCP ti
$HCPtieqp
1.000	0.000		0.000
-0.50	0.86602540378 	0.00
0.000	0.00		$HCPticoap
Ti
2
Direct
0.000000000000	0.000000000000	0.000000000000
0.333333333000	0.666666667000	0.500000000000
-

cat > $dir/tests/HCPti/__init__.py << !!
!!

cp $readhcpphon $dir/tests/HCPti/

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/HCPti/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.emt import EMT
from ase.units import _hbar, _e
from ase.io import read
import read_hcp_phonons as RHP
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/HCPti/OPT.POSCAR')

atoms.set_calculator(calc)

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

## High-symmetry points in the Brillouin zone
G = [0, 0, 0] 
H = [1./2, 0, 1./2]
K = [1./3, 1./3, 0]
M = [1./2, 0, 0]
A = [0, 0, 1./2]
L = [1./3, 1./3, 1./2]

point_names = ['\$\Gamma\$', '\$K\$', '\$M\$', '\$\Gamma\$', '\$A\$']
dirs = ['\$[\\\xi\\\xi0]\$', '' ,'\$[\\\xi00]\$', '\$[00\\\xi]\$']
path = [G, K, M, G, A]

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

exper = RHP.read_hcp_phonons('/n/jww-1/ehemann.2/testingScripts/EXP_DATA/Ti/hcp/phonons.dat', Q/np.pi)
dft = np.loadtxt('$dftdat/Ti/hcp/phonons.dat')
b = 2*np.pi/2.9392072
# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    omega_nd= dft[:, n+1]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_thisd= np.max(omega_nd)
    min_thisd= np.min(omega_nd)
    max_band = np.max([max_band, max_thisd, max_this])
    min_band = np.min([min_band, max_thisd, min_this])
    plt.plot(q, omega_n, 'r-', lw=2)
    plt.plot(b*dft[:,0], omega_nd, color='gray', linestyle='--', lw=2)

plt.errorbar(np.pi*exper[:,0], exper[:,1], yerr=exper[:,2], fmt='.', color='black')

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('hcp Ti')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='salmon', edgecolor='r', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/hcpTi_phonons.png')
ph.clean
!

python $dir/tests/HCPti/phonons.py > $dir/tests/HCPti/phonon_log
rm -f *.pckl
fi

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BCCti >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t BCC-Ti... ${non}"

cat > $dir/tests/BCCti/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

# define box variables
variable boxa equal $bccTilat
variable boxb equal $bccTilat
variable boxc equal $bccTilat
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		bcc $bccTilat

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/BCCti/eq.in > $dir/tests/BCCti/eq.out

BCCtieqp=`tail -1 $dir/tests/BCCti/eq.out | awk '{print $1}'`
BCCticoap=`tail -1 $dir/tests/BCCti/eq.out | awk '{print $2}'`
BCCtivop=`tail -1 $dir/tests/BCCti/eq.out | awk '{print $3}'`
BCCtipote=`tail -1 $dir/tests/BCCti/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for BCCti -----------------------------------------

cat > $dir/tests/BCCti/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCCti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $BCCtieqp
variable boxb equal $BCCtieqp
variable boxc equal $BCCtieqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		bcc $BCCtieqp
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/BCCti/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/BCCti/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/BCCti/evsv.in loop
!!
$lammps < $dir/tests/BCCti/evsv.in > $dir/tests/BCCti/evsv.out

sed -i '1d' $dir/tests/BCCti/evsv.dat
line=`python $BMFit $dir/tests/BCCti/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/BCCti/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	BCCtibulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	BCCtibulkp=0
fi

sed -i '1d' $dir/tests/BCCti/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/BCCti/pvsv.dat | sort -k2,2 -g | head -1`
	BCCtiPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${BCCtiPLAT[0]} $BCCtieqp"
		#BCCtieqp=`echo $line | awk '{print $3}'`
		#BCCtivop=`echo $line | awk '{print $1}'`
		#BCCtipote=`echo $line | awk '{print 1000*$4}'`
	fi
done
BCCtiPLAT["0"]=$BCCtieqp
awk -v Voo=$BCCtivop '{print $1/Voo, $2, $3}' $dir/tests/BCCti/pvsv.dat > tmp; mv tmp $dir/tests/BCCti/pvsv.dat

#---------------------------- elastic constants for bcc -----------------------------------
echo -e "${prp}Elastic constants:${non}"

echo "# pressure, c11, c12, c44" > $dir/tests/BCCti/C_VS_P.dat
for PRESS in $PRESSEQ; do 
echo "$PRESS GPa..."
for jj in $(seq 1 7); do

printf "\t $jj "

if [ $lammps_elcon_flag == "TRUE" ]; then
P=`echo $PRESS*10000 | bc -l`

case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac
cat > $dir/tests/BCCti/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms
variable boxa equal ${BCCtiPLAT["$PRESS"]}
variable boxb equal ${BCCtiPLAT["$PRESS"]}
variable boxc equal ${BCCtiPLAT["$PRESS"]}
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice bcc ${BCCtiPLAT["$PRESS"]} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region	box prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box 2 box
create_atoms ${idx["Ti"]} box
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/BCCti/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)
variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})


change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/BCCti/elcon.lin loop 
!!

$lammps < $dir/tests/BCCti/elcon.lin > $dir/tests/BCCti/elcon.out
sed -i '1d' $dir/tests/BCCti/stresses$jj.dat

else	# compute stress-strain curves with meamzilla
cat > $dir/tests/meamz_params <<@@
ngroups 1

optstyle powell
num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/BCCti/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
rm -f $dir/tests/BCCti/elcon$jj.conf
strain=`echo "$ecstr/100" | bc -l`
for i in $(seq -$DELPT $DELPT); do

	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del BCC ${BCCtiPLAT["$PRESS"]} conf >> $dir/tests/BCCti/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/BCCti/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=pc*$5; getline	
					syy=pc*$5; getline	
					szz=pc*$5; getline	
					sxy=pc*$5; getline	
					syz=pc*$5; getline	
					szx=pc*$5;

					print del, sxx, syy, szz, syz, szx, sxy

				     }' data.stress >> $dir/tests/BCCti/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""
rm -f $dir/tests/BCCti/EC_fits.dat; touch $dir/tests/BCCti/EC_fits.dat


# first row fits
awk '{print $1, $2}' $dir/tests/BCCti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/BCCti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/BCCti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/BCCti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/BCCti/stresses5.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/BCCti/stresses6.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/BCCti/stresses7.dat > tmp; python $fit tmp >> $dir/tests/BCCti/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/BCCti/EC_fits.dat` )
BCCtiC11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
BCCtiC12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
BCCtiC13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
BCCtiC22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
BCCtiC23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
BCCtiC33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
BCCtiC44i=${coup[12]}
BCCtiC55i=${coup[13]}
BCCtiC66i=${coup[14]}

BCCtiC11p=`python -c "print ($BCCtiC11i+$BCCtiC22i+$BCCtiC33i)/3"`
BCCtiC12p=`python -c "print ($BCCtiC12i+$BCCtiC23i+$BCCtiC13i)/3"`
BCCtiC44p=`python -c "print ($BCCtiC44i+$BCCtiC55i+$BCCtiC66i)/3"`

echo $PRESS $BCCtiC11p $BCCtiC12p $BCCtiC44p >> $dir/tests/BCCti/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	BCCtiC11=`printf '%3.f' $BCCtiC11p`
	BCCtiC12=`printf '%3.f' $BCCtiC12p`
	BCCtiC44=`printf '%3.f' $BCCtiC44p`

	C11BCCtid=`printf '%3.f' $C11BCCtid`
	C12BCCtid=`printf '%3.f' $C12BCCtid`
	C44BCCtid=`printf '%3.f' $C44BCCtid`
fi
done

#----------------------------- bct deformation for BCCti -----------------------------------------
echo "BCT deformation..."
cat > $dir/tests/BCCti/bct.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCCti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmin equal 0.80 
variable	dmax equal 1.60 
variable	jmax equal $nevpt

variable boxa equal $BCCtieqp
variable boxb equal $BCCtieqp
variable boxc equal $BCCtieqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		bcc $BCCtieqp
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable E0 equal \${Epa}
variable DE equal v_Epa-\${E0}
variable coa equal lz/\${boxa}

fix P all print 1 "\${coa} \${DE}" file $dir/tests/BCCti/bct.dat screen no title "# V/atom | Energy"

variable zd equal \${boxc}*\${dmin}
variable zf equal \${boxc}*\${dmax}
change_box all z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/BCCti/bct.in loop
!!
$lammps < $dir/tests/BCCti/bct.in > $dir/tests/BCCti/bct.out

## PHONONS FOR BCC ti
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/BCCti/OPT.POSCAR << -
HCP ti
$BCCtieqp
-0.5000	0.50000	0.50000
0.50000	-0.5000	0.50000
0.50000	0.50000	-0.5000
Ti
1
Direct
0.000000000000	0.000000000000	0.000000000000
-

cp $readbccphon $dir/tests/BCCti/

cat > $dir/tests/BCCti/__init__.py << --
--

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/BCCti/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.emt import EMT
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/BCCti/OPT.POSCAR')

atoms.set_calculator(calc)

# Phonon calculator
N = 10
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.0001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

points = ibz_points['bcc']
G = points['Gamma']
H = points['H']
N = points['N']
P = points['P']

path = [G, H, P, G, N]
point_names = ['\$\Gamma\$', '\$H\$', '\$P\$', '\$\Gamma\$', '\$N\$']
dirs = ['\$[00\\\xi]\$', '\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi0]\$']

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)
#omega_kn_dft = np.loadtxt("$dftdat/Ti/bcc/omega_kn.dat")
#qd           = np.loadtxt("$dftdat/Ti/bcc/q.dat")

dft = np.loadtxt("$dftdat/Ti/bcc/phonons.dat")
b = 2*np.pi/3.25414794
# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)

import read_bcc_phonons as RBP
exper = RBP.read_bcc_phonons("/n/jww-1/ehemann.2/testingScripts/EXP_DATA/Ti/bcc/phonons.dat", Q/np.pi)

np.savetxt('/n/jww-1/ehemann.2/testingScripts/binary/bcc_exp_phons.dat', exper)

# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    omega_nd= dft[:, n+1]
    max_this = max(np.max(omega_n), np.max(omega_nd))

    min_this = min(np.min(omega_n), np.min(omega_nd))
    max_band = np.max([max_band, max_this])
    min_band = np.min([min_band, min_this])
    plt.plot(q, omega_n, color='red', linestyle='-', lw=2)
    plt.plot(b*dft[:,0], omega_nd, color='gray', linestyle='--', lw=2)

plt.errorbar(np.pi*exper[:,0], conv*exper[:,1]/1000, xerr=exper[:,2], yerr=conv*exper[:,3]/1000, fmt='.', color='black')

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('bcc Ti')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=-1000, color='salmon', edgecolor='r', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/bccTi_phonons.png')
ph.clean
!

python $dir/tests/BCCti/phonons.py > $dir/tests/BCCti/phonon_log
rm -f *.pckl
fi
#--------------------------------------------------------------------------------------------


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FCCti >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t FCC-Ti... ${non}"

cat > $dir/tests/FCCti/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

# define box variables
variable boxa equal $fccTilat
variable boxb equal $fccTilat
variable boxc equal $fccTilat
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		fcc $fccTilat

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/FCCti/eq.in > $dir/tests/FCCti/eq.out

FCCtieqp=`tail -1 $dir/tests/FCCti/eq.out | awk '{print $1}'`
FCCticoap=`tail -1 $dir/tests/FCCti/eq.out | awk '{print $2}'`
FCCtivop=`tail -1 $dir/tests/FCCti/eq.out | awk '{print $3}'`
FCCtipote=`tail -1 $dir/tests/FCCti/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for FCCti -----------------------------------------

cat > $dir/tests/FCCti/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR FCCti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $FCCtieqp
variable boxb equal $FCCtieqp
variable boxc equal $FCCtieqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		fcc $FCCtieqp
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/FCCti/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/FCCti/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/FCCti/evsv.in loop
!!
$lammps < $dir/tests/FCCti/evsv.in > $dir/tests/FCCti/evsv.out

sed -i '1d' $dir/tests/FCCti/evsv.dat
line=`python $BMFit $dir/tests/FCCti/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/FCCti/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	FCCtibulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	FCCtibulkp=0
fi

sed -i '1d' $dir/tests/FCCti/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/FCCti/pvsv.dat | sort -k2,2 -g | head -1`
	FCCtiPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${FCCtiPLAT[0]} $FCCtieqp"
		#FCCtieqp=`echo $line | awk '{print $3}'`
		#FCCtivop=`echo $line | awk '{print $1}'`
		#FCCtipote=`echo $line | awk '{print 1000*$4}'`
	fi
done
FCCtiPLAT["0"]=$FCCtieqp
awk -v Voo=$FCCtivop '{print $1/Voo, $2, $3}' $dir/tests/FCCti/pvsv.dat > tmp; mv tmp $dir/tests/FCCti/pvsv.dat
#<><><><><><><><><><><>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< OMGti >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t omega-Ti... ${non}"

cat > $dir/tests/OMGti/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid
variable boxa equal $omgTilat
variable boxb equal $omgTilat*sqrt(3)/2
variable boxc equal $omgTilat*$omgTicoa
variable boxxy equal $omgTilat*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $omgTilat a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $omgTicoa &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large

create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0 scalexy yes scalexz yes scaleyz yes
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/OMGti/eq.in > $dir/tests/OMGti/eq.out

OMGtieqp=`tail -1 $dir/tests/OMGti/eq.out | awk '{print $1}'`
OMGticoap=`tail -1 $dir/tests/OMGti/eq.out | awk '{print $2}'`
OMGtivop=`tail -1 $dir/tests/OMGti/eq.out | awk '{print $3}'`
OMGtipote=`tail -1 $dir/tests/OMGti/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for OMGti -----------------------------------------

cat > $dir/tests/OMGti/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR OMGti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

#initialization variables
variable boxa equal $OMGtieqp
variable boxb equal $OMGtieqp*sqrt(3)/2
variable boxc equal $OMGtieqp*$OMGticoap
variable boxxy equal $OMGtieqp*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $OMGtieqp a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $OMGticoap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/OMGti/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/OMGti/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd  equal  \${boxa}*(1-v_dmax/2)
variable xf  equal  \${boxa}*(1+v_dmax/2)
variable yd  equal  \${boxb}*(1-v_dmax/2)
variable yf  equal  \${boxb}*(1+v_dmax/2)
variable zd  equal  \${boxc}*(1-v_dmax/2)
variable zf  equal  \${boxc}*(1+v_dmax/2)
variable xyd equal -0.5*\${xd}
variable xyf equal -0.5*\${xf}
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/OMGti/evsv.in loop
!!
$lammps < $dir/tests/OMGti/evsv.in > $dir/tests/OMGti/evsv.out

sed -i '1d' $dir/tests/OMGti/evsv.dat
line=`python $BMFit $dir/tests/OMGti/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/OMGti/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	OMGtibulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	OMGtibulkp=0
fi

sed -i '1d' $dir/tests/OMGti/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/OMGti/pvsv.dat | sort -k2,2 -g | head -1`
	OMGtiPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${OMGtiPLAT[0]} $OMGtieqp"
		#OMGtieqp=`echo $line | awk '{print $3}'`
		#OMGtivop=`echo $line | awk '{print $1}'`
		#OMGtipote=`echo $line | awk '{print 1000*$4}'`
	fi
done
OMGtiPLAT["0"]=$OMGtieqp
awk -v Voo=$OMGtivop '{print $1/Voo, $2, $3}' $dir/tests/OMGti/pvsv.dat > tmp; mv tmp $dir/tests/OMGti/pvsv.dat

# ---------------------- pressure dependence of lattice constants for omg ti ----------
echo -e "${prp}Pressure dependence of lattice constants...${non}"

rm -f $dir/tests/OMGti/abcvp.dat; echo "#a, c/a, gamma" > $dir/tests/OMGti/abcvp.dat
for PRESS in $PRESSLAT; do
PRES=`echo "10000*$PRESS" | bc -l`
cat > $dir/tests/OMGti/abcvp.lin << __
## LAMMPS script
units metal
boundary p p p
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable a1	equal	1.00
variable a2	equal	sqrt(3)
variable a3	equal	$OMGticoap
variable boxa equal $OMGtieqp*\${a1}
variable boxb equal $OMGtieqp*\${a2}
variable boxc equal $OMGtieqp*\${a3}

lattice custom $OMGtieqp &
	a1 \${a1}	0.0	0.0 &
	a2 0.0		\${a2}	0.0 &
	a3 0.0		0.0	\${a3} &
	basis 0.0	0.0	0.0 basis 0.5	0.5	0.0 &
	basis 0.0	0.33333 0.5 basis 0.5   0.83333 0.5 &
	basis 0.0	0.66667 0.5 basis 0.5   0.16667 0.5
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#region frz1 sphere 0.0 0.0 0.0 0.1 units box
#region frz2 sphere 0.5 0.5 0.0 0.1 units box
#region frz union 2 frz1 frz2
#group frz region frz
#
#fix FREEZE freeze frz

thermo 1000
thermo_style custom step temp press pe ke etotal lx ly lz
timestep 0.001

variable a equal lx
variable coa equal lz/lx
variable Gamma equal 180-2*180*atan(ly/lx)/3.141592654

fix P all press/berendsen aniso $PRES $PRES 1000.0
run 10000

dump D all custom 1 $dir/tests/OMGti/abcvp.coords xs ys zs

print "$PRESS \${a} \${coa} \${Gamma}"
__

$lammps < $dir/tests/OMGti/abcvp.lin > $dir/tests/OMGti/abcvp.out
tail -1 $dir/tests/OMGti/abcvp.out >> $dir/tests/OMGti/abcvp.dat

done

#---------------------------- elastic constants for omega ---------------------------------
echo -e "${prp}Elastic constants:${non}"

echo "# pressure, c11, c12, c13, c33, c44" > $dir/tests/OMGti/C_VS_P.dat
for PRESS in $PRESSEQ; do
echo "$PRESS GPa..."

for jj in $(seq 1 7); do
printf "\t $jj "

if [ $lammps_elcon_flag == "TRUE" ]; then
P=`echo $PRESS*10000 | bc -l`

# NOTE: THIS NEEDS TO BE ADJUSTED IF SUPERCELLS LARGER THAN 1X1X1 ARE TO BE USED
#case $jj in
#	1) XD="1+v_d"; YD="sqrt(3)*(1+v_d)/2"; ZD="1+v_d"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	2) XD="1+v_d"; YD="sqrt(3)*(1-v_d)/2"; ZD="1/(1-v_d2)"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	3) XD="1/(1-v_d2)"; YD="sqrt(3)*(1+v_d)/2"; ZD="1-v_d"; XYD="-0.5/(1-v_d2)"; XZD="0"; YZD="0" ;;
#	4) XD="1-v_d"; YD="sqrt(3)/(2*(1-v_d2))"; ZD="1+v_d"; XYD="-0.5*(1-v_d)"; XZD="0"; YZD="0" ;;
#	5) XD="1/(1-(v_d2)/4)"; YD="sqrt(3)/2"; ZD="1"; XYD="-0.5/(1-(v_d2)/4)"; XZD="0"; YZD="v_d" ;;
#	6) XD="1"; YD="sqrt(3)/(2*(1-(v_d2)/4))"; ZD="1"; XYD="-0.5"; XZD="v_d"; YZD="0" ;;
#	7) XD="1"; YD="sqrt(3)/2"; ZD="1/(1-(v_d2)/4)"; XYD="-0.5+v_d"; XZD="0"; YZD="0" ;;
#esac
case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac
cat > $dir/tests/OMGti/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms
variable boxa equal ${OMGtiPLAT["$PRESS"]}
variable boxb equal ${OMGtiPLAT["$PRESS"]}*0.86602540378
variable boxc equal ${OMGtiPLAT["$PRESS"]}*$OMGticoap
variable boxxy equal ${OMGtiPLAT["$PRESS"]}*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom ${OMGtiPLAT["$PRESS"]} a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $OMGticoap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
create_box 2 mybox
create_atoms ${idx["Ti"]} box
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel


# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-v_Sxx0)/vol/10000
variable Syy equal (c_sigyy-v_Syy0)/vol/10000
variable Szz equal (c_sigzz-v_Szz0)/vol/10000
variable Sxy equal (c_sigxy-v_Sxy0)/vol/10000
variable Sxz equal (c_sigxz-v_Sxz0)/vol/10000
variable Syz equal (c_sigyz-v_Syz0)/vol/10000

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/OMGti/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)
variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 1e-20 1e-20 100000 100000

run 1

next j
jump $dir/tests/OMGti/elcon.lin loop 
!!

$lammps < $dir/tests/OMGti/elcon.lin > $dir/tests/OMGti/elcon.out
sed -i '1d' $dir/tests/OMGti/stresses$jj.dat

else	# compute stress-strain curves with meamzilla
cat > $dir/tests/meamz_params <<@@
ngroups 1

optstyle powell
num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/OMGti/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
rm -f $dir/tests/OMGti/elcon$jj.conf
strain=`echo "$ecstr/100" | bc -l`
for i in $(seq -$DELPT $DELPT); do

	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del OMG ${OMGtiPLAT["$PRESS"]} conf $OMGticoap >> $dir/tests/OMGti/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/OMGti/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=pc*$5; getline	
					syy=pc*$5; getline	
					szz=pc*$5; getline	
					sxy=pc*$5; getline	
					syz=pc*$5; getline	
					szx=pc*$5;

					print del, sxx, syy, szz, syz, szx, sxy

				     }' data.stress >> $dir/tests/OMGti/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""
rm -f $dir/tests/OMGti/EC_fits.dat; touch $dir/tests/OMGti/EC_fits.dat


# first row fits
awk '{print $1, $2}' $dir/tests/OMGti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/OMGti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/OMGti/stresses1.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/OMGti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/OMGti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/OMGti/stresses2.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/OMGti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/OMGti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/OMGti/stresses3.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/OMGti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/OMGti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/OMGti/stresses4.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/OMGti/stresses5.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/OMGti/stresses6.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/OMGti/stresses7.dat > tmp; python $fit tmp >> $dir/tests/OMGti/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/OMGti/EC_fits.dat` )
OMGtiC11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
OMGtiC12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
OMGtiC13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
OMGtiC22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
OMGtiC23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
OMGtiC33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
OMGtiC44i=${coup[12]}
OMGtiC55i=${coup[13]}
OMGtiC66i=${coup[14]}

OMGtiC11p=`python -c "print ($OMGtiC11i + $OMGtiC22i)/2"`
OMGtiC13p=`python -c "print ($OMGtiC13i + $OMGtiC23i)/2"`
OMGtiC44p=`python -c "print ($OMGtiC44i + $OMGtiC55i)/2"`

echo "$PRESS $OMGtiC11p $OMGtiC12i $OMGtiC13p $OMGtiC33i $OMGtiC44p" >> $dir/tests/OMGti/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	OMGtiC11=`printf '%3.f' $OMGtiC11p`
	OMGtiC12=`printf '%3.f' $OMGtiC12i`
	OMGtiC13=`printf '%3.f' $OMGtiC13p`
	OMGtiC33=`printf '%3.f' $OMGtiC33i`
	OMGtiC44=`printf '%3.f' $OMGtiC44p`

	C11OMGtid=`printf '%3.f' $C11OMGtid`
	C12OMGtid=`printf '%3.f' $C12OMGtid`
	C13OMGtid=`printf '%3.f' $C13OMGtid`
	C33OMGtid=`printf '%3.f' $C33OMGtid`
	C44OMGtid=`printf '%3.f' $C44OMGtid`
fi
done

## PHONONS FOR OMEGA TI
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/OMGti/OPT.POSCAR << -
HCP ti
$OMGtieqp
1.000	0.000		0.000
-0.50	0.86602540378 	0.00
0.000	0.00		$OMGticoap
Ti
3
Direct
0.000000000	0.000000000	0.000000000
0.333333333	0.666666667	0.500000000
0.666666667	0.333333333	0.500000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/OMGti/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.emt import EMT
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/OMGti/OPT.POSCAR')

atoms.set_calculator(calc)



# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

## High-symmetry points in the Brillouin zone
G = [0, 0, 0] 
H = [1./2, 0, 1./2]
K = [1./3, 1./3, 0]
M = [1./2, 0, 0]
A = [0, 0, 1./2]
L = [1./3, 1./3, 1./2]

point_names = ['\$\Gamma\$', '\$K\$', '\$M\$', '\$\Gamma\$', '\$A\$']
dirs = ['\$[\\\xi\\\xi0]\$', '' ,'\$[\\\xi00]\$', '\$[00\\\xi]\$']
path = [G, K, M, G, A]

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

dft = np.loadtxt('$dftdat/Ti/omega/phonons.dat')
b = 2*np.pi/4.58018151

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    omega_nd= dft[:, n+1]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_thisd= np.max(omega_nd)
    min_thisd= np.min(omega_nd)
    max_band = np.max([max_band, max_thisd, max_this])
    min_band = np.min([min_band, max_thisd, min_this])
    plt.plot(q, omega_n, 'r-', lw=2)
    plt.plot(b*dft[:,0], omega_nd, color='gray', linestyle='--', lw=2)

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('\$\omega\$ Ti')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='salmon', edgecolor='r', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/omgTi_phonons.png')
ph.clean
!

python $dir/tests/OMGti/phonons.py > $dir/tests/OMGti/phonon_log
rm -f *.pckl
fi

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< A15 Ti >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t A15 Ti ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/A15ti/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM A15ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

variable SIZE equal $A15Tilat

lattice		custom $A15Tilat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} & 
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]}
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax iso 0.0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

run 0

print '\${a} \${vat} \${spe}'
!!

$lammps < $dir/tests/A15ti/eq.in > $dir/tests/A15ti/eq.out

A15tieqp=`tail -1 $dir/tests/A15ti/eq.out | awk '{print $1}'`
A15tivop=`tail -1 $dir/tests/A15ti/eq.out | awk '{print $2}'`
A15tipote=`tail -1 $dir/tests/A15ti/eq.out | awk '{print $3}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for A15ti ------------------------------------------

cat > $dir/tests/A15ti/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR A15ti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable SIZE equal $A15tieqp

lattice		custom $A15tieqp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} & 
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]}

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/A15ti/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/A15ti/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/A15ti/evsv.in loop
!!
$lammps < $dir/tests/A15ti/evsv.in > $dir/tests/A15ti/evsv.out

sed -i '1d' $dir/tests/A15ti/evsv.dat
line=`python $BMFit $dir/tests/A15ti/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15ti/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15tibulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15tibulkp=0
fi

else # A15ti lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/A15ti/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

A15timmz=`echo "$A15tilat" | bc -l`
echo $A15timmz
rm -f $dir/tests/A15ti/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$A15timmz
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0,	0,	0,	0,	0,	0,	0
	print 0,	0.5,	0.5,	0.5,	0,	0,	0
	print 0,	0.5,	0,	0.25,	0,	0,	0
	print 0,	0.5,	0,	0.75,	0,	0,	0
	print 0,	0,	0.25,	0.5,	0,	0,	0
	print 0,	0,	0.75,	0.5,	0,	0,	0
	print 0,	0.25,	0.5,	0,	0,	0,	0
	print 0,	0.75,	0.5,	0,	0,	0,	0
" > $dir/tests/A15ti/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15ti/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15ti/pvsv.dat
awk -v a0=$A15timmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol;
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15ti/pvsv.dat 

echo "#v, a" > $dir/tests/A15ti/evsv.dat
awk -v a0=$A15timmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/8, $6;
				      }' data.energy >> $dir/tests/A15ti/evsv.dat 


sed -i '1d' $dir/tests/A15ti/evsv.dat
line=`python $BMFit $dir/tests/A15ti/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15ti/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15tibulkp=`echo $line | awk '{print $1}'`
	A15tivop=`echo $line | awk '{print $2}'`
	A15tipote=`echo $line | awk '{print $3}'`
	A15tieqp=`python -c "print (8.*$A15tivop)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15tibulkp=0
	A15tivop=0
	A15tipote=0
	A15tieqp=0
fi

echo -e "${prp}E-V curve...${non}"

A15timmz=`echo "$A15tieqp" | bc -l`
rm -f $dir/tests/A15ti/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $A15timmz*math.pow(a,1./3)
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0,	0,	0,	0,	0,	0,	0
	print 0,	0.5,	0.5,	0.5,	0,	0,	0
	print 0,	0.5,	0,	0.25,	0,	0,	0
	print 0,	0.5,	0,	0.75,	0,	0,	0
	print 0,	0,	0.25,	0.5,	0,	0,	0
	print 0,	0,	0.75,	0.5,	0,	0,	0
	print 0,	0.25,	0.5,	0,	0,	0,	0
	print 0,	0.75,	0.5,	0,	0,	0,	0

" > $dir/tests/A15ti/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15ti/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15ti/pvsv.dat
awk -v a0=$A15timmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15ti/pvsv.dat 

echo "#v, a" > $dir/tests/A15ti/evsv.dat
awk -v a0=$A15timmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/8.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/A15ti/evsv.dat 
fi # lammps_evsv_flag A15ti

sed -i '1d' $dir/tests/A15ti/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/A15ti/pvsv.dat | sort -k2,2 -g | head -1`
	A15tiPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${A15tiPLAT[0]} $A15tieqp"
		#A15tieqp=`echo $line | awk '{print $3}'`
		#A15tivop=`echo $line | awk '{print $1}'`
		#A15tipote=`echo $line | awk '{print 1000*$4}'`
	fi
done
A15tiPLAT["0"]=$A15tieqp

awk -v Voo=$A15tivop '{print $1/Voo, $2, $3}' $dir/tests/A15ti/pvsv.dat > tmp; mv tmp $dir/tests/A15ti/pvsv.dat




#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BCCnb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t BCC-Nb... ${non}"

cat > $dir/tests/BCCnb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

variable boxa equal $bccNblat
variable boxb equal $bccNblat
variable boxc equal $bccNblat
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		bcc $bccNblat

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/BCCnb/eq.in > $dir/tests/BCCnb/eq.out

BCCnbeqp=`tail -1 $dir/tests/BCCnb/eq.out | awk '{print $1}'`
BCCnbcoap=`tail -1 $dir/tests/BCCnb/eq.out | awk '{print $1}'`
BCCnbvop=`tail -1 $dir/tests/BCCnb/eq.out | awk '{print $3}'`
BCCnbpote=`tail -1 $dir/tests/BCCnb/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for BCCnb -----------------------------------------

cat > $dir/tests/BCCnb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCCnb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $BCCnbeqp
variable boxb equal $BCCnbeqp
variable boxc equal $BCCnbeqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0


lattice		bcc $BCCnbeqp
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/BCCnb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/BCCnb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/BCCnb/evsv.in loop
!!
$lammps < $dir/tests/BCCnb/evsv.in > $dir/tests/BCCnb/evsv.out

sed -i '1d' $dir/tests/BCCnb/evsv.dat
line=`python $BMFit $dir/tests/BCCnb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/BCCnb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	BCCnbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	BCCnbbulkp=0
fi

sed -i '1d' $dir/tests/BCCnb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/BCCnb/pvsv.dat | sort -k2,2 -g | head -1`
	BCCnbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${BCCnbPLAT[0]} $BCCnbeqp"
		#BCCnbeqp=`echo $line | awk '{print $3}'`
		#BCCnbvop=`echo $line | awk '{print $1}'`
		#BCCnbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
BCCnbPLAT["0"]=$BCCnbeqp
awk -v Voo=$BCCnbvop '{print $1/Voo, $2, $3}' $dir/tests/BCCnb/pvsv.dat > tmp; mv tmp $dir/tests/BCCnb/pvsv.dat

#---------------------------- stacking faults for bcc nb-------------------------------------
if [ "$hcp_gamma_flag" == "TRUE" ]; then
echo "Stacking Faults..."
if [ "$gamma_lammps_flag" == "TRUE" ]; then
for PRESS in 0; do
for jj in $(seq 1 2); do
case $jj in
	1) XLAT="111"; ZLAT="1-10"; FILE=$dir/tests/BCCnb/"$PRESS"SF_110_111.dat; t=0;
	   RILE=$dir/tests/BCCnb/"$PRESS"SF_110_111_relaxed.dat; lw="ly"; tf="xy"; reg="INF INF \${BOT} \${TOP} INF INF";; # symmetric
	2) XLAT="111"; ZLAT="11-2"; FILE=$dir/tests/BCCnb/"$PRESS"SF_112_111.dat; t=1;
	   RILE=$dir/tests/BCCnb/"$PRESS"SF_112_111_relaxed.dat; lw="lz"; tf="xz"; reg="INF INF INF INF \${BOT} \${TOP}";; # asymmetric
esac
printf "$jj "

rm -f $FILE; rm -f $RILE
python $gammaTilt BCC ${BCCnbPLAT["$PRESS"]} $ZLAT $XLAT 0.0 lmp > $dir/tests/BCCnb/read_gamma.dat
declare -a tilts=( `python $gammaTilt BCC ${BCCnbPLAT["$PRESS"]} $ZLAT $XLAT 0.5 lmp | awk 'NR==7{print 2*$1, 2*$2, 2*$3}'` )
for i in $(seq 0 $gampts); do
deli=`echo "$i/$gampts" | bc -l`
python $gammaTilt BCC ${BCCnbPLAT["$PRESS"]} $ZLAT $XLAT $deli lmp > $dir/tests/BCCnb/read_gamma.dat.dat
sed -i 's/1 atom types/2 atom types/g' $dir/tests/BCCnb/read_gamma.dat
if [ "$i" == "0" ]; then area=`awk 'NR==1{print $7}' read_me.dat`; fi
cat > $dir/tests/BCCnb/gamma_in << __
units metal
boundary p p p
atom_style atomic
box tilt large
read_data $dir/tests/BCCnb/read_gamma.dat
mass 1 $mass1
mass 2 $mass2
pair_style $pairstyle
pair_coeff * * $potPath $elem1 $elem2
set group all type 2
thermo 1
thermo_style custom pe ke etotal
timestep 0.01
variable EATOM equal pe
run 0
print "EINIT: \${EATOM}"
fix frz all setforce 0.0 0.0 NULL
min_style fire
minimize 1e-8 1e-8 5000 50000
run 0
print "EATOM: \${EATOM}"
__
$lammps < $dir/tests/BCCnb/gamma_in > $dir/tests/BCCnb/gamma_out
if [ "$i" == "0" ]; then e0=`grep 'EATOM:' $dir/tests/BCCnb/gamma_out | awk '{print $2}'`; mevpaa=0;
			eu0=`grep 'EINIT:' $dir/tests/BCCnb/gamma_out | awk '{print $2}'`; mevini=0;
else
mevpaa=`grep 'EATOM:' $dir/tests/BCCnb/gamma_out | awk -v a=$area -v e0=$e0 '{print 1000*($2-e0)/a}'`
mevini=`grep 'EINIT:' $dir/tests/BCCnb/gamma_out | awk -v a=$area -v e0=$eu0 '{print 1000*($2-e0)/a}'`
fi

echo $deli $mevpaa >> $RILE
echo $deli $mevini >> $FILE

done
done
done
echo ""
else

cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 1

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/BCCnb/gamma.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

for PRESS in $(seq 0 25 0); do
P=`echo "10000*$PRESS" | bc -l`

for jj in $(seq 1 2); do
case $jj in
	1) XLAT="111"; ZLAT="110"; FILE=$dir/tests/BCCnb/"$PRESS"SF_110_111.dat;;		# {110}<111>
	2) XLAT="111"; ZLAT="112"; FILE=$dir/tests/BCCnb/"$PRESS"SF_112_111.dat;;		# {112}<111>
esac
printf "$jj "

echo "# DEL, GAMMA" > $FILE
rm -f $dir/tests/BCCnb/gamma.conf; touch $dir/tests/BCCnb/gamma.conf
for d in $(seq 0 $gampts); do

	DEL=`echo "scale=6; $d/$gampts" | bc -l`
	python $gammaTilt BCC ${BCCnbPLAT["$PRESS"]} $ZLAT $XLAT $DEL conf >> $dir/tests/BCCnb/gamma.conf

done

$meamz -p $dir/tests/meamz_params > $dir/tests/BCCnb/meamzilla.out

e0=`awk 'NR==3{print $6}' data.energy`
area=`awk 'NR==2{print $3}' $dir/tests/BCCnb/gamma.conf`
tail -n +3 data.energy | awk -v npts=$gampts -v e0=$e0 -v a=$area '{print $1/npts, 1000*$4*($6-e0)/a}' >> $FILE

rm data.*

done
done

bccSFsymp=`grep '^0.5 ' $dir/tests/BCCnb/0SF_110_111.dat | awk '{print $2}'`
bccSFasyp=`grep '^0.5 ' $dir/tests/BCCnb/0SF_112_111.dat | awk '{print $2}'`
echo ""
fi
fi

#---------------------------- elastic constants for bcc Nb -----------------------------------
echo -e "${prp}Elastic constants:${non}"

echo "# pressure, c11, c12, c44" > $dir/tests/BCCnb/C_VS_P.dat
for PRESS in $PRESSEQ; do 
echo "$PRESS GPa..."
for jj in $(seq 1 7); do

printf "\t $jj "

if [ $lammps_elcon_flag == "TRUE" ]; then
P=`echo $PRESS*10000 | bc -l`

case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac

cat > $dir/tests/BCCnb/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms

variable boxa equal ${BCCnbPLAT["$PRESS"]}
variable boxb equal ${BCCnbPLAT["$PRESS"]}
variable boxc equal ${BCCnbPLAT["$PRESS"]}
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice bcc ${BCCnbPLAT["$PRESS"]} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region box prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box 2 box
create_atoms ${idx["Nb"]} box
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/BCCnb/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)
variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/BCCnb/elcon.lin loop 
!!

$lammps < $dir/tests/BCCnb/elcon.lin > $dir/tests/BCCnb/elcon.out
sed -i '1d' $dir/tests/BCCnb/stresses$jj.dat

else	# compute stress-strain curves with meamzilla
cat > $dir/tests/meamz_params <<@@
ngroups 1

optstyle powell
num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/BCCnb/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
strain=`echo "$ecstr/100" | bc -l`
rm -f tmp
rm -f $dir/tests/BCCnb/elcon$jj.conf; touch $dir/tests/BCCnb/elcon$jj.conf
for i in $(seq -$DELPT $DELPT); do

	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del BCC ${BCCnbPLAT["$PRESS"]} conf > tmp
	head -1 tmp | awk '{print $1, $2, 2}' > tmp2
	head -8 tmp | tail -7 >> tmp2
	tail -2 tmp | awk '{print 1, $2, $3, $4, $5, $6, $7}' >> tmp2
	cat tmp2 >> $dir/tests/BCCnb/elcon$jj.conf
	rm -f tmp tmp2 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/BCCnb/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=pc*$5; getline	
					syy=pc*$5; getline	
					szz=pc*$5; getline	
					sxy=pc*$5; getline	
					syz=pc*$5; getline	
					szx=pc*$5;

					print del, sxx, syy, szz, syz, szx, sxy

				     }' data.stress >> $dir/tests/BCCnb/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""
rm -f $dir/tests/BCCnb/EC_fits.dat; touch $dir/tests/BCCnb/EC_fits.dat


# first row fits
awk '{print $1, $2}' $dir/tests/BCCnb/stresses1.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCnb/stresses1.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCnb/stresses1.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/BCCnb/stresses2.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCnb/stresses2.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCnb/stresses2.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/BCCnb/stresses3.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCnb/stresses3.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCnb/stresses3.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/BCCnb/stresses4.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/BCCnb/stresses4.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/BCCnb/stresses4.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/BCCnb/stresses5.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/BCCnb/stresses6.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/BCCnb/stresses7.dat > tmp; python $fit tmp >> $dir/tests/BCCnb/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/BCCnb/EC_fits.dat` )
BCCnbC11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
BCCnbC12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
BCCnbC13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
BCCnbC22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
BCCnbC23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
BCCnbC33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
BCCnbC44i=${coup[12]}
BCCnbC55i=${coup[13]}
BCCnbC66i=${coup[14]}

BCCnbC11p=`python -c "print ($BCCnbC11i+$BCCnbC22i+$BCCnbC33i)/3"`
BCCnbC12p=`python -c "print ($BCCnbC12i+$BCCnbC23i+$BCCnbC13i)/3"`
BCCnbC44p=`python -c "print ($BCCnbC44i+$BCCnbC55i+$BCCnbC66i)/3"`

echo $PRESS $BCCnbC11p $BCCnbC12p $BCCnbC44p >> $dir/tests/BCCnb/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	BCCnbC11=`printf '%3.f' $BCCnbC11p`
	BCCnbC12=`printf '%3.f' $BCCnbC12p`
	BCCnbC44=`printf '%3.f' $BCCnbC44p`

	C11BCCnbd=`printf '%3.f' $C11BCCnbd`
	C12BCCnbd=`printf '%3.f' $C12BCCnbd`
	C44BCCnbd=`printf '%3.f' $C44BCCnbd`
fi
done
#--------------------------------------------------------------------------------------------
#----------------------------- bct deformation for BCCnb -----------------------------------------
echo "BCT deformation..."
cat > $dir/tests/BCCnb/bct.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCCnb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmin equal 0.80 
variable	dmax equal 1.60 
variable	jmax equal $nevpt

variable boxa equal $BCCnbeqp
variable boxb equal $BCCnbeqp
variable boxc equal $BCCnbeqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		bcc $BCCnbeqp
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable E0 equal \${Epa}
variable DE equal v_Epa-\${E0}
variable coa equal lz/\${boxa}

fix P all print 1 "\${coa} \${DE}" file $dir/tests/BCCnb/bct.dat screen no title "# V/atom | Energy"

variable zd equal \${boxc}*\${dmin}
variable zf equal \${boxc}*\${dmax}
change_box all z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/BCCnb/bct.in loop
!!
$lammps < $dir/tests/BCCnb/bct.in > $dir/tests/BCCnb/bct.out

## PHONONS FOR BCC nb
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/BCCnb/OPT.POSCAR << -
BCC Nb
$BCCnbeqp
-0.5000	0.50000	0.50000
0.50000	-0.5000	0.50000
0.50000	0.50000	-0.5000
Nb
1
Direct
0.000000000000	0.000000000000	0.000000000000
-

cat > $dir/tests/BCCnb/__init__.py << -
-

cp $readbccphon $dir/tests/BCCnb/

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/BCCnb/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.emt import EMT
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/BCCnb/OPT.POSCAR')

atoms.set_calculator(calc)

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

points = ibz_points['bcc']
G = points['Gamma']
H = points['H']
N = points['N']
P = points['P']

path = [G, H, P, G, N]
point_names = ['\$\Gamma\$', '\$H\$', '\$P\$', '\$\Gamma\$', '\$N\$']
dirs = ['\$[00\\\xi]\$', '\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi0]\$']

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)
dft = np.loadtxt("$dftdat/Nb/bcc/phonons.dat")
b = 2*np.pi/3.309

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)

import read_bcc_phonons as RBP
exper = RBP.read_bcc_phonons("/n/jww-1/ehemann.2/testingScripts/EXP_DATA/Nb/bcc/phonons.dat", Q/np.pi)

# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    omega_nd= dft[:,n+1]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_thisd= np.max(omega_nd)
    min_thisd= np.min(omega_nd)
    max_band = np.max([max_band, max_this, max_thisd])
    min_band = np.min([min_band, min_this, min_thisd])
    plt.plot(q, omega_n, 'b-', lw=2)
    plt.plot(b*dft[:,0], omega_nd, color='gray', linestyle='--', lw=2)

plt.errorbar(np.pi*exper[:,0], conv*exper[:,1]/1000, xerr=exper[:,2], yerr=conv*exper[:,3]/1000, fmt='.', color='black')

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('bcc Nb')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='lightblue', edgecolor='b', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/bccNb_phonons.png')
ph.clean
!

python $dir/tests/BCCnb/phonons.py > $dir/tests/BCCnb/phonon_log
rm -f *.pckl
fi
#--------------------------------------------------------------------------------------------
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HCPnb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t HCP-Nb... ${non}"

cat > $dir/tests/HCPnb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid
variable boxa equal $hcpNblat
variable boxb equal $hcpNblat*sqrt(3)/2
variable boxc equal $hcpNblat*$hcpNbcoa
variable boxxy equal $hcpNblat*(-0.5)

lattice		custom $hcpNblat a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $hcpNbcoa &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
		
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0 scalexy yes scalexz yes scaleyz yes
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/HCPnb/eq.in > $dir/tests/HCPnb/eq.out

HCPnbeqp=`tail -1 $dir/tests/HCPnb/eq.out | awk '{print $1}'`
HCPnbcoap=`tail -1 $dir/tests/HCPnb/eq.out | awk '{print $2}'`
HCPnbvop=`tail -1 $dir/tests/HCPnb/eq.out | awk '{print $3}'`
HCPnbpote=`tail -1 $dir/tests/HCPnb/eq.out | awk '{print $4}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for HCPnb -----------------------------------------

cat > $dir/tests/HCPnb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR HCPnb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $HCPnbeqp
variable boxb equal $HCPnbeqp*sqrt(3)/2
variable boxc equal $HCPnbeqp*$HCPnbcoap
variable boxxy equal $HCPnbeqp*(-0.5)

lattice		custom $HCPnbeqp a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $HCPnbcoap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
		
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

dump D all xyz 1 $dir/tests/HCPnb/evsv.xyz
fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/HCPnb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/HCPnb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd  equal  \${boxa}*(1-v_dmax/2)
variable xf  equal  \${boxa}*(1+v_dmax/2)
variable yd  equal  \${boxb}*(1-v_dmax/2)
variable yf  equal  \${boxb}*(1+v_dmax/2)
variable zd  equal  \${boxc}*(1-v_dmax/2)
variable zf  equal  \${boxc}*(1+v_dmax/2)
variable xyd equal -0.5*\${xd}
variable xyf equal -0.5*\${xf}
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/HCPnb/evsv.in loop
!!
$lammps < $dir/tests/HCPnb/evsv.in > $dir/tests/HCPnb/evsv.out

sed -i '1d' $dir/tests/HCPnb/evsv.dat
line=`python $BMFit $dir/tests/HCPnb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/HCPnb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	HCPnbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	HCPnbbulkp=0
fi

sed -i '1d' $dir/tests/HCPnb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/HCPnb/pvsv.dat | sort -k2,2 -g | head -1`
	HCPnbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${HCPnbPLAT[0]} $HCPnbeqp"
		#HCPnbeqp=`echo $line | awk '{print $3}'`
		#HCPnbvop=`echo $line | awk '{print $1}'`
		#HCPnbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
HCPnbPLAT["0"]=$HCPnbeqp
awk -v Voo=$HCPnbvop '{print $1/Voo, $2, $3}' $dir/tests/HCPnb/pvsv.dat > tmp; mv tmp $dir/tests/HCPnb/pvsv.dat


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FCCnb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t FCC-Nb... ${non}"

cat > $dir/tests/FCCnb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

# define box variables
variable boxa equal $fccNblat
variable boxb equal $fccNblat
variable boxc equal $fccNblat
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		fcc $fccNblat

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/FCCnb/eq.in > $dir/tests/FCCnb/eq.out

FCCnbeqp=`tail -1 $dir/tests/FCCnb/eq.out | awk '{print $1}'`
FCCnbcoap=`tail -1 $dir/tests/FCCnb/eq.out | awk '{print $2}'`
FCCnbvop=`tail -1 $dir/tests/FCCnb/eq.out | awk '{print $3}'`
FCCnbpote=`tail -1 $dir/tests/FCCnb/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for FCCnb -----------------------------------------

cat > $dir/tests/FCCnb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR FCCnb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $FCCnbeqp
variable boxb equal $FCCnbeqp
variable boxc equal $FCCnbeqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		fcc $FCCnbeqp
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/FCCnb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/FCCnb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/FCCnb/evsv.in loop
!!
$lammps < $dir/tests/FCCnb/evsv.in > $dir/tests/FCCnb/evsv.out

sed -i '1d' $dir/tests/FCCnb/evsv.dat
line=`python $BMFit $dir/tests/FCCnb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/FCCnb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	FCCnbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	FCCnbbulkp=0
fi

sed -i '1d' $dir/tests/FCCnb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/FCCnb/pvsv.dat | sort -k2,2 -g | head -1`
	FCCnbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${FCCnbPLAT[0]} $FCCnbeqp"
		#FCCnbeqp=`echo $line | awk '{print $3}'`
		#FCCnbvop=`echo $line | awk '{print $1}'`
		#FCCnbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
FCCnbPLAT["0"]=$FCCnbeqp
awk -v Voo=$FCCnbvop '{print $1/Voo, $2, $3}' $dir/tests/FCCnb/pvsv.dat > tmp; mv tmp $dir/tests/FCCnb/pvsv.dat

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< OMGnb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t omega Nb... ${non}"

cat > $dir/tests/OMGnb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid
variable boxa equal $omgNblat
variable boxb equal $omgNblat*sqrt(3)/2
variable boxc equal $omgNblat*$omgNbcoa
variable boxxy equal $omgNblat*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $omgNblat a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $omgNbcoa &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large

create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0 scalexy yes scalexz yes scaleyz yes
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/OMGnb/eq.in > $dir/tests/OMGnb/eq.out

OMGnbeqp=`tail -1 $dir/tests/OMGnb/eq.out | awk '{print $1}'`
OMGnbcoap=`tail -1 $dir/tests/OMGnb/eq.out | awk '{print $2}'`
OMGnbvop=`tail -1 $dir/tests/OMGnb/eq.out | awk '{print $3}'`
OMGnbpote=`tail -1 $dir/tests/OMGnb/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for OMGnb -----------------------------------------

cat > $dir/tests/OMGnb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR OMGnb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

#initialization variables
variable boxa equal $OMGnbeqp
variable boxb equal $OMGnbeqp*sqrt(3)/2
variable boxc equal $OMGnbeqp*$OMGnbcoap
variable boxxy equal $OMGnbeqp*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $OMGnbeqp a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $OMGnbcoap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/OMGnb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/OMGnb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd  equal  \${boxa}*(1-v_dmax/2)
variable xf  equal  \${boxa}*(1+v_dmax/2)
variable yd  equal  \${boxb}*(1-v_dmax/2)
variable yf  equal  \${boxb}*(1+v_dmax/2)
variable zd  equal  \${boxc}*(1-v_dmax/2)
variable zf  equal  \${boxc}*(1+v_dmax/2)
variable xyd equal -0.5*\${xd}
variable xyf equal -0.5*\${xf}
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/OMGnb/evsv.in loop
!!
$lammps < $dir/tests/OMGnb/evsv.in > $dir/tests/OMGnb/evsv.out

sed -i '1d' $dir/tests/OMGnb/evsv.dat
line=`python $BMFit $dir/tests/OMGnb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/OMGnb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	OMGnbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	OMGnbbulkp=0
fi

sed -i '1d' $dir/tests/OMGnb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/OMGnb/pvsv.dat | sort -k2,2 -g | head -1`
	OMGnbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${OMGnbPLAT[0]} $OMGnbeqp"
		#OMGnbeqp=`echo $line | awk '{print $3}'`
		#OMGnbvop=`echo $line | awk '{print $1}'`
		#OMGnbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
OMGnbPLAT["0"]=$OMGnbeqp
awk -v Voo=$OMGnbvop '{print $1/Voo, $2, $3}' $dir/tests/OMGnb/pvsv.dat > tmp; mv tmp $dir/tests/OMGnb/pvsv.dat

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< A15 Nb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t A15 Nb ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/A15nb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM A15nb LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

variable SIZE equal $A15nblat

lattice		custom $A15nblat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Nb"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} & 
		basis 5 ${idx["Nb"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Nb"]} basis 8 ${idx["Nb"]}
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax iso 0.0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

run 0

print '\${a} \${vat} \${spe}'
!!

$lammps < $dir/tests/A15nb/eq.in > $dir/tests/A15nb/eq.out

A15nbeqp=`tail -1 $dir/tests/A15nb/eq.out | awk '{print $1}'`
A15nbvop=`tail -1 $dir/tests/A15nb/eq.out | awk '{print $2}'`
A15nbpote=`tail -1 $dir/tests/A15nb/eq.out | awk '{print $3}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for A15nb ------------------------------------------

cat > $dir/tests/A15nb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR A15nb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable SIZE equal $A15nbeqp

lattice		custom $A15nbeqp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Nb"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} & 
		basis 5 ${idx["Nb"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Nb"]} basis 8 ${idx["Nb"]}

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/A15nb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/A15nb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/A15nb/evsv.in loop
!!
$lammps < $dir/tests/A15nb/evsv.in > $dir/tests/A15nb/evsv.out

sed -i '1d' $dir/tests/A15nb/evsv.dat
line=`python $BMFit $dir/tests/A15nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15nb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15nbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15nbbulkp=0
fi

else # A15nb lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/A15nb/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

A15nbmmz=`echo "$A15nblat" | bc -l`
echo $A15nbmmz
rm -f $dir/tests/A15nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$A15nbmmz
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.5,	0.5,	0.5,	0,	0,	0
	print 1,	0.5,	0,	0.25,	0,	0,	0
	print 1,	0.5,	0,	0.75,	0,	0,	0
	print 1,	0,	0.25,	0.5,	0,	0,	0
	print 1,	0,	0.75,	0.5,	0,	0,	0
	print 1,	0.25,	0.5,	0,	0,	0,	0
	print 1,	0.75,	0.5,	0,	0,	0,	0
" > $dir/tests/A15nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15nb/pvsv.dat
awk -v a0=$A15nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol;
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15nb/pvsv.dat 

echo "#v, a" > $dir/tests/A15nb/evsv.dat
awk -v a0=$A15nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/8, $6;
				      }' data.energy >> $dir/tests/A15nb/evsv.dat 


sed -i '1d' $dir/tests/A15nb/evsv.dat
line=`python $BMFit $dir/tests/A15nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15nb/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15nbbulkp=`echo $line | awk '{print $1}'`
	A15nbvop=`echo $line | awk '{print $2}'`
	A15nbpote=`echo $line | awk '{print $3}'`
	A15nbeqp=`python -c "print (8.*$A15nbvop)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15nbbulkp=0
	A15nbvop=0
	A15nbpote=0
	A15nbeqp=0
fi

echo -e "${prp}E-V curve...${non}"

A15nbmmz=`echo "$A15nbeqp" | bc -l`
rm -f $dir/tests/A15nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $A15nbmmz*math.pow(a,1./3)
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.5,	0.5,	0.5,	0,	0,	0
	print 1,	0.5,	0,	0.25,	0,	0,	0
	print 1,	0.5,	0,	0.75,	0,	0,	0
	print 1,	0,	0.25,	0.5,	0,	0,	0
	print 1,	0,	0.75,	0.5,	0,	0,	0
	print 1,	0.25,	0.5,	0,	0,	0,	0
	print 1,	0.75,	0.5,	0,	0,	0,	0

" > $dir/tests/A15nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15nb/pvsv.dat
awk -v a0=$A15nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15nb/pvsv.dat 

echo "#v, a" > $dir/tests/A15nb/evsv.dat
awk -v a0=$A15nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/8.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/A15nb/evsv.dat 
fi # lammps_evsv_flag A15nb

sed -i '1d' $dir/tests/A15nb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/A15nb/pvsv.dat | sort -k2,2 -g | head -1`
	A15nbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${A15nbPLAT[0]} $A15nbeqp"
		#A15nbeqp=`echo $line | awk '{print $3}'`
		#A15nbvop=`echo $line | awk '{print $1}'`
		#A15nbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
A15nbPLAT["0"]=$A15nbeqp

awk -v Voo=$A15nbvop '{print $1/Voo, $2, $3}' $dir/tests/A15nb/pvsv.dat > tmp; mv tmp $dir/tests/A15nb/pvsv.dat



#################################################################
#################################################################

#################################################################
#	BCC TESTS
#################################################################
echo -e "${grn}FINISHED SINGLE ELEMENTS, BEGINNING ALLOY TESTS... ${non}"
echo -e "${grn}BEGINNING BCC TESTS... ${non}"

if [ "$SS_FLAG" == "TRUE" ]; then

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t Solid solution Ti-Nb... ${non}"
SSaeqp=3.25
SSbeqp=3.00
SSvop=17.1641

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curves...${non}"
seed=$RANDOM
dimalpha=12
dimbeta=24
# loop over niobium concentrations
for CONC in 25 33 50 66 75; do
	echo "$CONC % Niobium ..."
	case $CONC in
		25)		conc=$CONC;    tag="7525";;
		33)		conc=33.33333; tag="6733";;
		50)		conc=$CONC;    tag="5050";;
		66)		conc=66.66667; tag="3367";;
		75)		conc=$CONC;    tag="2575";;
	esac

for struct in bcc hcp; do

	case $struct in
		hcp) sr3=`echo "sqrt(3)" | bc -l`
		     latcom="lattice custom $SSaeqp a1 1.0 0.0 0.0 a2 0.0 $sr3 0.0 a3 0.0 0.0 $HCPticoap &
				basis 0.0 0.0 0.0 basis 0.5 0.5 0.0 basis 0.0 0.6667 0.5 basis 0.5 0.16667 0.5"
		     boxa=$dimalpha*$SSaeqp; boxb=$dimalpha*$sr3*$SSaeqp; boxc=$dimalpha*$HCPticoap*$SSaeqp;;
		bcc) sr2=`echo "sqrt(2)" | bc -l`
		     latcom="lattice custom $SSbeqp a1 1.0 0.0 0.0 a2 0.0 $sr2 0.0 a3 0.0 0.0 $sr2 &
				basis 0.0 0.0 0.0 basis 0.5 0.5 0.0 basis 0.5 0.0 0.5 basis 0.0 0.5 0.5"
		     boxa=$dimbeta*$SSbeqp; boxb=$dimbeta*$sr2*$SSbeqp; boxc=$dimbeta*$sr2*$SSbeqp;;
	esac
bigtag=`python -c "print '$struct' + '_' + '$tag'"`
#----------------------------- energy volume for SS -----------------------------------------

cat > $dir/tests/SS/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR SS LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable	seed equal $seed
variable	conc equal $conc/100

variable boxa equal $boxa 
variable boxb equal $boxb
variable boxc equal $boxc

$latcom	
	
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box
set group all type/fraction 2 \${conc} \${seed}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 0
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix REL all press/berendsen aniso 0.0 0.0 1000.0
run 2000
unfix REL
#velocity all zero linear
change_box all triclinic
#thermo 1
#fix REL all press/berendsen aniso 0.0 0.0 1000.0
#fix REL all box/relax tri 0.0 fixedpoint 0 0 0
#min_style cg
#min_modify line forcezero
#minimize 0.0 1e-6 5000 50000
#unfix REL 
#run 0

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/SS/evsv_$bigtag.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/SS/pvsv_$bigtag.dat screen no title "# V/atom | pressure | lattice constant"

variable tmp equal lx
variable aa equal \${tmp}
variable tmp equal ly
variable bb equal \${tmp}
variable tmp equal lz
variable cc equal \${tmp}
variable tmp equal xy
variable xyo equal \${tmp}/\${aa}
variable tmp equal xz
variable xzo equal \${tmp}/\${aa}
variable tmp equal yz
variable yzo equal \${tmp}/\${bb}

variable xd equal \${aa}*(1-v_dmax/2)
variable xf equal \${aa}*(1+v_dmax/2)
variable yd equal \${bb}*(1-v_dmax/2)
variable yf equal \${bb}*(1+v_dmax/2)
variable zd equal \${cc}*(1-v_dmax/2)
variable zf equal \${cc}*(1+v_dmax/2)
variable xyd equal \${xyo}*\${xd} 
variable xyf equal \${xyo}*\${xf} 
variable xzd equal \${xzo}*\${xd} 
variable xzf equal \${xzo}*\${xf} 
variable yzd equal \${yzo}*\${yd} 
variable yzf equal \${yzo}*\${yf} 

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} xz final \${xzf} yz final \${yzf} units box

reset_timestep 0
run 0

run \${jmax}
!!
$lammps < $dir/tests/SS/evsv.in > $dir/tests/SS/evsv.out

sed -i '1d' $dir/tests/SS/evsv_$bigtag.dat
line=`python $BMFit $dir/tests/SS/evsv_$bigtag.dat`
#line=`echo "data=Import[\"$dir/tests/SS/evsv_$bigtag.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	SSbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	SSbulkp=0
fi

SSBULK["$CONC"]=$SSbulkp
sed -i '1d' $dir/tests/SS/pvsv_$bigtag.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/SS/pvsv_$bigtag.dat | sort -k2,2 -g | head -1`
	SSPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${SSPLAT[0]}"
		#SSeqp=`echo $line | awk '{print $3}'`
		#SSvop=`echo $line | awk '{print $1}'`
		#SSpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
#SSPLAT["0"]=$SSeqp
awk -v Voo=$SSvop '{print $1/Voo, $2, $3}' $dir/tests/SS/pvsv_$bigtag.dat > tmp; mv tmp $dir/tests/SS/pvsv_$bigtag.dat
done
done
#---------------------------------------------------------------

fi

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< B2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t B2-TiNb... ${non}"

cat > $dir/tests/B2/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

# define box variables
variable boxa equal $B2lat
variable boxb equal $B2lat
variable boxc equal $B2lat
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		custom $B2lat &
		a1 1.00 0.00 0.00 a2 0.00 1.00 0.00 a3 0.00 0.00 1.00 &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.50

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/B2/eq.in > $dir/tests/B2/eq.out

B2eqp=`tail -1 $dir/tests/B2/eq.out | awk '{print $1}'`
B2coap=`tail -1 $dir/tests/B2/eq.out | awk '{print $2}'`
B2vop=`tail -1 $dir/tests/B2/eq.out | awk '{print $3}'`
B2pote=`tail -1 $dir/tests/B2/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for B2 -----------------------------------------

cat > $dir/tests/B2/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR B2 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $B2eqp
variable boxb equal $B2eqp
variable boxc equal $B2eqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $B2eqp &
		a1 1.00 0.00 0.00 a2 0.00 1.00 0.00 a3 0.00 0.00 1.00 &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.50
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/B2/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/B2/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${xf} z final 0 \${xf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/B2/evsv.in loop
!!
$lammps < $dir/tests/B2/evsv.in > $dir/tests/B2/evsv.out

sed -i '1d' $dir/tests/B2/evsv.dat
line=`python $BMFit $dir/tests/B2/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/B2/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	B2bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	B2bulkp=0
fi

sed -i '1d' $dir/tests/B2/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/B2/pvsv.dat | sort -k2,2 -g | head -1`
	B2PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${B2PLAT[0]} $B2eqp"
		#B2eqp=`echo $line | awk '{print $3}'`
		#B2vop=`echo $line | awk '{print $1}'`
		#B2pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
B2PLAT["0"]=$B2eqp
awk -v Voo=$B2vop '{print $1/Voo, $2, $3}' $dir/tests/B2/pvsv.dat > tmp; mv tmp $dir/tests/B2/pvsv.dat

#---------------------------- elastic constants for b2-----------------------------------
echo -e "${prp}Elastic constants:${non}"

echo "# pressure, c11, c12, c44" > $dir/tests/B2/C_VS_P.dat
for PRESS in 0; do 
echo "$PRESS GPa..."
for jj in $(seq 1 7); do

printf "\t $jj "

if [ $lammps_elcon_flag == "TRUE" ]; then
P=`echo $PRESS*10000 | bc -l`

case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac

cat > $dir/tests/B2/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms

variable boxa equal ${B2PLAT["$PRESS"]}
variable boxb equal ${B2PLAT["$PRESS"]}
variable boxc equal ${B2PLAT["$PRESS"]}
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom ${B2PLAT["$PRESS"]} &
		a1 1.00 0.00 0.00 a2 0.00 1.00 0.00 a3 0.00 0.00 1.00 &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.50 &
		orient x 1 0 0 orient y 0 1 0 orient z 0 0 1

region box prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box 2 box
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]}
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/B2/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)
variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/B2/elcon.lin loop 
!!

$lammps < $dir/tests/B2/elcon.lin > $dir/tests/B2/elcon.out
sed -i '1d' $dir/tests/B2/stresses$jj.dat

else	# compute stress-strain curves with meamzilla
cat > $dir/tests/meamz_params <<@@
ngroups 1

optstyle powell
num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/B2/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
strain=`echo "$ecstr/100" | bc -l`
rm -f tmp
rm -f $dir/tests/B2/elcon$jj.conf; touch $dir/tests/B2/elcon$jj.conf
for i in $(seq -$DELPT $DELPT); do

	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del BCC ${B2PLAT["$PRESS"]} conf > tmp
	head -1 tmp | awk '{print $1, $2, 2}' > tmp2
	head -8 tmp | tail -7 >> tmp2
	tail -2 tmp | awk '{print 1, $2, $3, $4, $5, $6, $7}' >> tmp2
	cat tmp2 >> $dir/tests/B2/elcon$jj.conf
	rm -f tmp tmp2 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/B2/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=pc*$5; getline	
					syy=pc*$5; getline	
					szz=pc*$5; getline	
					sxy=pc*$5; getline	
					syz=pc*$5; getline	
					szx=pc*$5;

					print del, sxx, syy, szz, syz, szx, sxy

				     }' data.stress >> $dir/tests/B2/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""
rm -f $dir/tests/B2/EC_fits.dat; touch $dir/tests/B2/EC_fits.dat


# first row fits
awk '{print $1, $2}' $dir/tests/B2/stresses1.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/B2/stresses1.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/B2/stresses1.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/B2/stresses2.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/B2/stresses2.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/B2/stresses2.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/B2/stresses3.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/B2/stresses3.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/B2/stresses3.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/B2/stresses4.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/B2/stresses4.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/B2/stresses4.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/B2/stresses5.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/B2/stresses6.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/B2/stresses7.dat > tmp; python $fit tmp >> $dir/tests/B2/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/B2/EC_fits.dat` )
B2C11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
B2C12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
B2C13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
B2C22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
B2C23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
B2C33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
B2C44i=${coup[12]}
B2C55i=${coup[13]}
B2C66i=${coup[14]}

B2C11p=`python -c "print ($B2C11i+$B2C22i+$B2C33i)/3"`
B2C12p=`python -c "print ($B2C12i+$B2C23i+$B2C13i)/3"`
B2C44p=`python -c "print ($B2C44i+$B2C55i+$B2C66i)/3"`

echo $PRESS $B2C11p $B2C12p $B2C44p >> $dir/tests/B2/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	B2C11=`printf '%3.f' $B2C11p`
	B2C12=`printf '%3.f' $B2C12p`
	B2C44=`printf '%3.f' $B2C44p`

	C11B2d=`printf '%3.f' $C11B2d`
	C12B2d=`printf '%3.f' $C12B2d`
	C44B2d=`printf '%3.f' $C44B2d`
fi
done
#--------------------------------------------------------------------------------------------

#----------------------------- bct deformation for B2 -----------------------------------------
echo "BCT deformation..."
cat > $dir/tests/B2/bct.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR B2 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmin equal 0.80 
variable	dmax equal 1.60 
variable	jmax equal $nevpt

variable boxa equal $B2eqp
variable boxb equal $B2eqp
variable boxc equal $B2eqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $B2eqp &
		a1 1.00 0.00 0.00 a2 0.00 1.00 0.00 a3 0.00 0.00 1.00 &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.50
		
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable E0 equal \${Epa}
variable DE equal v_Epa-\${E0}
variable coa equal lz/\${boxa}

fix P all print 1 "\${coa} \${DE}" file $dir/tests/B2/bct.dat screen no title "# V/atom | Energy"

variable zd equal \${boxc}*\${dmin}
variable zf equal \${boxc}*\${dmax}
change_box all z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/B2/bct.in loop
!!
$lammps < $dir/tests/B2/bct.in > $dir/tests/B2/bct.out

#----------------------------- point defects B2 -----------------------------------------
echo "point defects..."
declare -A B2_pds
declare -A B2_pdd
declare -A B2_pderr
HB2=`python -c "print (2*$B2pote-$HCPtipote-$BCCnbpote)/2"`
relrad=`echo "1.01*sqrt(2)" | bc -l`
for pd in 'Ti-vacancy' 'Nb-vacancy' 'Ti-antisite' 'Nb-antisite'; do
echo -e "\t $pd..."
case $pd in
	'Ti-vacancy')	pos="0.0 0.0 0.0"; cmd="delete_atoms region DEFECT";;
	'Nb-vacancy')	pos="0.5 0.5 0.5"; cmd="delete_atoms region DEFECT";;
	'Ti-antisite')	pos="0.0 0.0 0.0"; cmd="set region DEFECT type 2"  ;; # antisite on titanium sublattice
	'Nb-antisite')	pos="0.5 0.5 0.5"; cmd="set region DEFECT type 1"  ;; # antisite on niobium  sublattice
esac

cat > $dir/tests/B2/point_defects.in << !!
units metal
atom_style atomic

#define box variables
variable SIZE equal 5
variable lat  equal $B2eqp
variable boxa equal \${SIZE}*\${lat}
variable boxb equal \${SIZE}*\${lat}
variable boxc equal \${SIZE}*\${lat}

lattice		custom \${lat} &
		a1 1.00 0.00 0.00 a2 0.00 1.00 0.00 a3 0.00 0.00 1.00 &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.50

region box block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box 2 box
create_atoms ${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]}

region DEFECT sphere $pos 0.1 units lattice
$cmd

group TI type 1
group NB type 2

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

timestep 0.001
thermo 10
thermo_style custom step pe ke etotal press

variable peat equal pe
region FRZME sphere $pos $relrad side out
group FREEZE region FRZME
fix frz FREEZE setforce 0.0 0.0 0.0
min_style fire
minimize 1e-8 0.0 1000 10000

variable numTi equal count(TI)
variable numNb equal count(NB)

print "\${peat} \${numTi} \${numNb}"
!!

$lammps < $dir/tests/B2/point_defects.in > $dir/tests/B2/point_defects_"$pd".out
line=`tail -1 $dir/tests/B2/point_defects_"$pd".out`
e=`echo $line | awk '{print $1}'`
declare -a ats=( `echo $line | awk '{print $2, $3}'` )
nats=`echo "${ats[0]}+${ats[1]}" | bc`
if [ "$pd" == "Ti-vacancy" ]; then
#	#ats[0]=$((${ats[0]}+1))
	nats=$(($nats+1))
elif [ "$pd" == "Nb-vacancy" ]; then
#	#ats[1]=$((${ats[1]}+1))
	nats=$(($nats+1))
fi

B2_pds["$pd"]=`python -c "print $e-${ats[0]}*$HCPtipote-${ats[1]}*$BCCnbpote-$nats*$HB2"`
B2_pdd["$pd"]=`grep "$pd" $dftdat/Ti-Nb/TiNb/B2/point_defects.dat | awk '{print $2}'`
B2_pderr["$pd"]=`python -c "print '%.1f' % (100*(${B2_pds["$pd"]}-${B2_pdd["$pd"]})/${B2_pdd["$pd"]})"`

done

for complex in exchange divacancy triple-Ti triple-Nb inter-Ti inter-Nb; do

	case "$complex" in
		"exchange")
			e=`python -c "print ${B2_pds['Ti-antisite']}+${B2_pds['Nb-antisite']}"` ;;
		"divacancy")
			e=`python -c "print ${B2_pds['Ti-vacancy']}+${B2_pds['Nb-vacancy']}"` ;;
		"triple-Ti")
			e=`python -c "print 2*${B2_pds['Ti-vacancy']}+${B2_pds['Nb-antisite']}"` ;;
		"triple-Nb")
			e=`python -c "print 2*${B2_pds['Nb-vacancy']}+${B2_pds['Ti-antisite']}"` ;;
		"inter-Ti")
			e=`python -c "print 2*${B2_pds['Nb-vacancy']}-${B2_pds['Nb-antisite']}"` ;;
		"inter-Nb")
			e=`python -c "print -2*${B2_pds['Ti-vacancy']}+${B2_pds['Ti-antisite']}"` ;;
	esac
	B2_pds["$complex"]=$e	
	B2_pdd["$complex"]=`grep "$complex" $dftdat/Ti-Nb/TiNb/B2/point_defects.dat | awk '{print $2}'`
	B2_pderr["$complex"]=`python -c "print '%.1f' % (100*(${B2_pds["$complex"]}-${B2_pdd["$complex"]})/${B2_pdd["$complex"]})"`
done

## PHONONS FOR B2 ti
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/B2/OPT.POSCAR << -
B2-Tinb
$B2eqp
1.00	0.00	0.00
0.00	1.00	0.00
0.00	0.00	1.00
Ti Nb
1 1
Direct
0.000000000000	0.000000000000	0.000000000000
0.500000000000	0.500000000000	0.500000000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/B2/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/B2/OPT.POSCAR')

atoms.set_calculator(calc)

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

G = [0, 0, 0]
X = [1./2, 0, 0]
M = [1./2, 1./2, 0]
R = [1./2, 1./2, 1./2]


point_names = ['\$R\$', '\$\Gamma\$', '\$M\$', '\$X\$', '\$\Gamma\$']
path = [R, G, M, X, G]
dirs = ['\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi0]\$', '\$[0\\\bar{\\\xi}0]\$', '\$[\\\xi00]\$']


# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_band = np.max([max_band, max_this])
    min_band = np.min([min_band, min_this])
    plt.plot(q, omega_n, 'y-', lw=2)

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('B2 TiNb')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='yellow', edgecolor='y', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/B2_phonons.png')
ph.clean
!

python $dir/tests/B2/phonons.py > $dir/tests/B2/phonon_log
rm -f *.pckl
fi
#--------------------------------------------------------------------------------------------

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< {110}-layer bcc >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t {110}-layered bcc TiNb... ${non}"

cat > $dir/tests/bcc110/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

# define box variables
variable boxa equal $bcc110lat
variable boxb equal sqrt(2)*$bcc110lat
variable boxc equal sqrt(2)*$bcc110lat
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0
variable SQRT2 equal sqrt(2)

#define simulation region and bcc grid
lattice		custom $bcc110lat &
		a1 1.00 0.00 0.00 a2 0.00 \${SQRT2} 0.00 a3 0.00 0.00 \${SQRT2} &
		basis 0.00 0.00 0.00 basis 0.50 0.00 0.50 &
		basis 0.00 0.50 0.50 basis 0.50 0.50 0.00

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/bcc110/eq.in > $dir/tests/bcc110/eq.out

bcc110eqp=`tail -1 $dir/tests/bcc110/eq.out | awk '{print $1}'`
bcc110coap=`tail -1 $dir/tests/bcc110/eq.out | awk '{print $2}'`
bcc110vop=`tail -1 $dir/tests/bcc110/eq.out | awk '{print $3}'`
bcc110pote=`tail -1 $dir/tests/bcc110/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for bcc110 -----------------------------------------

cat > $dir/tests/bcc110/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR bcc110 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $bcc110eqp
variable boxb equal sqrt(2)*$bcc110eqp
variable boxc equal sqrt(2)*$bcc110eqp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0
variable SQRT2 equal sqrt(2)

lattice		custom $bcc110eqp &
		a1 1.00 0.00 0.00 a2 0.00 \${SQRT2} 0.00 a3 0.00 0.00 \${SQRT2} &
		basis 0.00 0.00 0.00 basis 0.50 0.00 0.50 &
		basis 0.00 0.50 0.50 basis 0.50 0.50 0.00
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/bcc110/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/bcc110/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/bcc110/evsv.in loop
!!
$lammps < $dir/tests/bcc110/evsv.in > $dir/tests/bcc110/evsv.out

sed -i '1d' $dir/tests/bcc110/evsv.dat
line=`python $BMFit $dir/tests/bcc110/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/bcc110/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	bcc110bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	bcc110bulkp=0
fi

sed -i '1d' $dir/tests/bcc110/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/bcc110/pvsv.dat | sort -k2,2 -g | head -1`
	bcc110PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${bcc110PLAT[0]} $bcc110eqp"
		#bcc110eqp=`echo $line | awk '{print $3}'`
		#bcc110vop=`echo $line | awk '{print $1}'`
		#bcc110pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
bcc110PLAT["0"]=$bcc110eqp
awk -v Voo=$bcc110vop '{print $1/Voo, $2, $3}' $dir/tests/bcc110/pvsv.dat > tmp; mv tmp $dir/tests/bcc110/pvsv.dat

#################################################################
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< A3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t A3-TiNb... ${non}"

cat > $dir/tests/A3/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

# define box variables
variable boxa equal $A3lat
variable a2y  equal sqrt(3)/2
variable boxb equal \${a2y}*$A3lat
variable a3z  equal $A3coa
variable boxc equal \${a3z}*$A3lat
variable a2x  equal 0.5
variable boxxy equal \${a2x}*$A3lat
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		custom $A3lat &
		a1 1.00 0.00 0.00 a2 \${a2x} \${a2y} 0.00 a3 0.00 0.00 \${a3z} &
		basis 0.00 0.00 0.00 basis 0.33333 0.33333 0.50

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/A3/eq.in > $dir/tests/A3/eq.out

A3eqp=`tail -1 $dir/tests/A3/eq.out | awk '{print $1}'`
A3coap=`tail -1 $dir/tests/A3/eq.out | awk '{print $2}'`
A3vop=`tail -1 $dir/tests/A3/eq.out | awk '{print $3}'`
A3pote=`tail -1 $dir/tests/A3/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for A3 -----------------------------------------

cat > $dir/tests/A3/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR A3 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

# define box variables
variable boxa equal $A3eqp
variable a2y  equal sqrt(3)/2
variable boxb equal \${a2y}*$A3eqp
variable a3z  equal $A3coap
variable boxc equal \${a3z}*$A3eqp
variable a2x  equal 0.5
variable boxxy equal \${a2x}*$A3eqp
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		custom $A3eqp &
		a1 1.00 0.00 0.00 a2 \${a2x} \${a2y} 0.00 a3 0.00 0.00 \${a3z} &
		basis 0.00 0.00 0.00 basis 0.33333 0.33333 0.50

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/A3/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/A3/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
variable xyd equal 0.5*\${xd}
variable xyf equal 0.5*\${xf}

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/A3/evsv.in loop
!!
$lammps < $dir/tests/A3/evsv.in > $dir/tests/A3/evsv.out

sed -i '1d' $dir/tests/A3/evsv.dat
line=`python $BMFit $dir/tests/A3/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A3/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A3bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A3bulkp=0
fi

sed -i '1d' $dir/tests/A3/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/A3/pvsv.dat | sort -k2,2 -g | head -1`
	A3PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${A3PLAT[0]} $A3eqp"
		#A3eqp=`echo $line | awk '{print $3}'`
		#A3vop=`echo $line | awk '{print $1}'`
		#A3pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
A3PLAT["0"]=$A3eqp
awk -v Voo=$A3vop '{print $1/Voo, $2, $3}' $dir/tests/A3/pvsv.dat > tmp; mv tmp $dir/tests/A3/pvsv.dat

## PHONONS FOR A3
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/A3/OPT.POSCAR << -
A3 tinb
$A3eqp
1.000	0.000		0.000
-0.50	0.86602540378 	0.00
0.000	0.00		$A3coap
Ti Nb
1 1
Direct
0.000000000000	0.000000000000	0.000000000000
0.333333333000	0.666666667000	0.500000000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/A3/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.emt import EMT
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/A3/OPT.POSCAR')

atoms.set_calculator(calc)

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

## High-symmetry points in the Brillouin zone
G = [0, 0, 0] 
H = [1./2, 0, 1./2]
K = [1./3, 1./3, 0]
M = [1./2, 0, 0]
A = [0, 0, 1./2]
L = [1./3, 1./3, 1./2]

point_names = ['\$\Gamma\$', '\$K\$', '\$M\$', '\$\Gamma\$', '\$A\$']
dirs = ['\$[\\\xi\\\xi0]\$', '' ,'\$[\\\xi00]\$', '\$[00\\\xi]\$']
path = [G, K, M, G, A]

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_band = np.max([max_band, max_this])
    min_band = np.min([min_band, min_this])
    plt.plot(q, omega_n, 'y-', lw=2)

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('A3 TiNb')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='yellow', edgecolor='y', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/A3_phonons.png')
ph.clean
!

python $dir/tests/A3/phonons.py > $dir/tests/A3/phonon_log
rm -f *.pckl
fi


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< D03 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t D03 ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

D03lmp=`echo "2*$D03lat" | bc -l`

cat > $dir/tests/D03/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM BCC-D03 LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and grid

variable boxa equal $D03lmp
variable boxb equal $D03lmp
variable boxc equal $D03lmp

lattice		custom $D03lmp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Ti"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Nb"]} basis 15 ${idx["Nb"]} basis 16 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx/2
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/D03/eq.in > $dir/tests/D03/eq.out

D03eqp=`tail -1 $dir/tests/D03/eq.out | awk '{print $1}'`
D03vop=`tail -1 $dir/tests/D03/eq.out | awk '{print $3}'`
D03pote=`tail -1 $dir/tests/D03/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for D03 -----------------------------------------

D03lmp=`echo "2*$D03eqp" | bc -l`

cat > $dir/tests/D03/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCC-D03 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $D03lmp
variable boxb equal $D03lmp
variable boxc equal $D03lmp

lattice		custom $D03lmp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Ti"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Nb"]} basis 15 ${idx["Nb"]} basis 16 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx/2
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/D03/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/D03/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/D03/evsv.in loop
!!
$lammps < $dir/tests/D03/evsv.in > $dir/tests/D03/evsv.out

sed -i '1d' $dir/tests/D03/evsv.dat
line=`python $BMFit $dir/tests/D03/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/D03/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	D03bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	D03bulkp=0
fi


else # D03lammps flag

## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/D03/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

D03mmz=`echo "2*$D03lat" | bc -l`
rm -f $dir/tests/D03/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)*$D03mmz*$D03mmz*$D03mmz
	a = math.pow(a,1./3)
	
	print '#N', 16, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, a/2, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, a/2, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, 0.0, a/2, 0.0, 0.0, 0.0
	print 0, a/2, a/2, a/2, 0.0, 0.0, 0.0
	print 0, a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 1, a/2, a/2, 0.0, 0.0, 0.0, 0.0
	print 1, a/2, 0.0, a/2, 0.0, 0.0, 0.0
	print 1, 0.0, a/2, a/2, 0.0, 0.0, 0.0
" > $dir/tests/D03/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/D03/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/D03/pvsv.dat
awk -v a0=$D03mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/16, pum, a;
				     }' data.stress >> $dir/tests/D03/pvsv.dat 

echo "#v, a" > $dir/tests/D03/evsv.dat
awk -v a0=$D03mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/16, $6;
				      }' data.energy >> $dir/tests/D03/evsv.dat 


sed -i '1d' $dir/tests/D03/evsv.dat
line=`python $BMFit $dir/tests/D03/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/D03/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	D03bulkp=`echo $line | awk '{print $1}'`
	D03vop=`echo $line | awk '{print $2}'`
	D03pote=`echo $line | awk '{print $3}'`
	D03eqp=`python -c "print (2.*$D03vop)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	D03bulkp=0
	D03vop=0
	D03pote=0
	D03eqp=0
fi

echo -e "${prp}E-V curve...${non}"

D03mmz=`echo "2*$D03eqp" | bc -l`
rm -f $dir/tests/D03/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)*$D03mmz*$D03mmz*$D03mmz
	a = math.pow(a,1./3)
	
	print '#N', 16, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, a/2, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, a/2, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, 0.0, a/2, 0.0, 0.0, 0.0
	print 0, a/2, a/2, a/2, 0.0, 0.0, 0.0
	print 0, a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 1, a/2, a/2, 0.0, 0.0, 0.0, 0.0
	print 1, a/2, 0.0, a/2, 0.0, 0.0, 0.0
	print 1, 0.0, a/2, a/2, 0.0, 0.0, 0.0

" > $dir/tests/D03/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/D03/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/D03/pvsv.dat
awk -v a0=$D03mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/16, pum, a;
				     }' data.stress >> $dir/tests/D03/pvsv.dat 

echo "#v, a" > $dir/tests/D03/evsv.dat
awk -v a0=$D03mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/16.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/D03/evsv.dat 
fi # lammps_evsv_flag D03

sed -i '1d' $dir/tests/D03/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/D03/pvsv.dat | sort -k2,2 -g | head -1`
	D03PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${D03PLAT[0]} $D03eqp"
		#D03eqp=`echo $line | awk '{print $3}'`
		#D03vop=`echo $line | awk '{print $1}'`
		#D03pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
D03PLAT["0"]=$D03eqp
awk -v Voo=$D03vop '{print $1/Voo, $2, $3}' $dir/tests/D03/pvsv.dat > tmp; mv tmp $dir/tests/D03/pvsv.dat

#---------------------------- elastic constants for D03 -----------------------------------
echo -e "${prp}Elastic constants:${non}"

echo "# pressure, c11, c12, c44" > $dir/tests/D03/C_VS_P.dat
for PRESS in $PRESSEQ; do 
echo "$PRESS GPa..."
D03lmp=`echo 2*${D03PLAT["$PRESS"]} | bc -l`
for jj in $(seq 1 7); do

printf "\t $jj "

if [ $lammps_elcon_flag == "TRUE" ]; then

P=`echo $PRESS*10000 | bc -l`
case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac

cat > $dir/tests/D03/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms
variable boxa equal $D03lmp
variable boxb equal $D03lmp
variable boxc equal $D03lmp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $D03lmp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 

region box prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box

create_box	2 box
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Ti"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Nb"]} basis 15 ${idx["Nb"]} basis 16 ${idx["Nb"]}

mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/D03/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)

variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/D03/elcon.lin loop 
!!

$lammps < $dir/tests/D03/elcon.lin > $dir/tests/D03/elcon.out
sed -i '1d' $dir/tests/D03/stresses$jj.dat

else

cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/D03/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
strain=0.002
rm -f $dir/tests/D03/elcon$jj.conf
for i in $(seq -$DELPT $DELPT); do
	D03mmz=`echo 2*${D03PLAT["$PRESS"]} | bc -l`
	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del D03 $D03mmz conf >> $dir/tests/D03/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/D03/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=-pc*$5; getline	
					syy=-pc*$5; getline	
					szz=-pc*$5; getline	
					sxy=-pc*$5; getline	
					syz=-pc*$5; getline	
					szx=-pc*$5; 	
					print del, sxx, syy, szz, syz, szx, sxy
				     }' data.stress >> $dir/tests/D03/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""

rm -f $dir/tests/D03/EC_fits.dat; touch $dir/tests/D03/EC_fits.dat

# first row fits
awk '{print $1, $2}' $dir/tests/D03/stresses1.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/D03/stresses1.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/D03/stresses1.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/D03/stresses2.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/D03/stresses2.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/D03/stresses2.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/D03/stresses3.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/D03/stresses3.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/D03/stresses3.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/D03/stresses4.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/D03/stresses4.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/D03/stresses4.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/D03/stresses5.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/D03/stresses6.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/D03/stresses7.dat > tmp; python $fit tmp >> $dir/tests/D03/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/D03/EC_fits.dat` )
D03C11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
D03C12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
D03C13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
D03C22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
D03C23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
D03C33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
D03C44i=${coup[12]}
D03C55i=${coup[13]}
D03C66i=${coup[14]}

D03C11p=`python -c "print ($D03C11i+$D03C22i+$D03C33i)/3"`
D03C12p=`python -c "print ($D03C12i+$D03C23i+$D03C13i)/3"`
D03C44p=`python -c "print ($D03C44i+$D03C55i+$D03C66i)/3"`

echo $PRESS $D03C11p $D03C12p $D03C44p >> $dir/tests/D03/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	D03C11=`printf '%3.f' $D03C11p`
	D03C12=`printf '%3.f' $D03C12p`
	D03C44=`printf '%3.f' $D03C44p`

	C11D03d=`printf '%3.f' $C11D03d`
	C12D03d=`printf '%3.f' $C12D03d`
	C44D03d=`printf '%3.f' $C44D03d`
fi
done

#----------------------------- bct deformation for D03 -----------------------------------------
echo "BCT deformation..."
D03lmp=`echo "2*$D03eqp" | bc -l`
cat > $dir/tests/D03/bct.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR D03 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmin equal 0.80 
variable	dmax equal 1.60 
variable	jmax equal $nevpt

variable boxa equal $D03lmp
variable boxb equal $D03lmp
variable boxc equal $D03lmp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $D03lmp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Ti"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Nb"]} basis 15 ${idx["Nb"]} basis 16 ${idx["Nb"]}
		
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable E0 equal \${Epa}
variable DE equal v_Epa-\${E0}
variable coa equal lz/\${boxa}

fix P all print 1 "\${coa} \${DE}" file $dir/tests/D03/bct.dat screen no title "# V/atom | Energy"

variable zd equal \${boxc}*\${dmin}
variable zf equal \${boxc}*\${dmax}
change_box all z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/D03/bct.in loop
!!
$lammps < $dir/tests/D03/bct.in > $dir/tests/D03/bct.out

## PHONONS FOR D03
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

D03mmz=`echo "2*$D03lat" | bc -l`
cat > $dir/tests/D03/OPT.POSCAR << -
D03 Ti3Nb
$D03mmz
0.50	0.50	0.0
0.50	0.00	-0.50
0.00	0.50	-0.50
Ti Nb
3 1
Direct
0.250000	0.250000	0.250000
0.500000	0.500000	0.500000
0.750000	0.750000	0.750000
0.000000	0.000000	0.000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/D03/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.emt import EMT
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/D03/OPT.POSCAR')

atoms.set_calculator(calc)

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

G = [0, 0, 0]
X = [1./2, 0, 0]
M = [1./2, 1./2, 0]
R = [1./2, 1./2, 1./2]

point_names = ['', '\$\Gamma\$', '', '', '\$\Gamma\$']
path = [R, G, M, X, G]
dirs = ['\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi0]\$', '\$[0\\\bar{\\\xi}0]\$', '\$[\\\xi00]\$']

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_band = np.max([max_band, max_this])
    min_band = np.min([min_band, min_this])
    plt.plot(q, omega_n, 'g-', lw=2)

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('D0\$_{3}\$ Ti\$_{3}\$Nb')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='lightgreen', edgecolor='g', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/D03_phonons.png')
ph.clean
!

python $dir/tests/D03/phonons.py > $dir/tests/D03/phonon_log
rm -f *.pckl
fi

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< G1  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t G1 ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

G1lmp=`echo "2*$G1lat" | bc -l`

cat > $dir/tests/G1/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM BCC-G1 LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid
variable boxa equal $G1lmp
variable boxb equal $G1lmp
variable boxc equal $G1lmp

lattice		custom $G1lmp  &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Ti"]} basis 15 ${idx["Ti"]} basis 16 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx/2
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/G1/eq.in > $dir/tests/G1/eq.out

G1eqp=`tail -1 $dir/tests/G1/eq.out | awk '{print $1}'`
G1vop=`tail -1 $dir/tests/G1/eq.out | awk '{print $3}'`
G1pote=`tail -1 $dir/tests/G1/eq.out | awk '{print $4}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for G1 ------------------------------------------

G1lmp=`echo "2*$G1eqp" | bc -l`

cat > $dir/tests/G1/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCC-G1 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $G1lmp
variable boxb equal $G1lmp
variable boxc equal $G1lmp

lattice		custom $G1lmp  &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Ti"]} basis 15 ${idx["Ti"]} basis 16 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx/2
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/G1/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/G1/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/G1/evsv.in loop
!!
$lammps < $dir/tests/G1/evsv.in > $dir/tests/G1/evsv.out

sed -i '1d' $dir/tests/G1/evsv.dat
line=`python $BMFit $dir/tests/G1/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/G1/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	G1bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	G1bulkp=0
fi

sed -i '1d' $dir/tests/G1/pvsv.dat




else # G1 lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/G1/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

G1mmz=`echo "2*$G1lat" | bc -l`
rm -f $dir/tests/G1/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)*$G1mmz*$G1mmz*$G1mmz
	a = math.pow(a,1./3)
	
	print '#N', 16, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, a/2, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, a/2, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, 0.0, a/2, 0.0, 0.0, 0.0
	print 1, a/2, a/2, a/2, 0.0, 0.0, 0.0
	print 0, a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 1, 3*a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 1, a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, a/2, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, 0.0, a/2, 0.0, 0.0, 0.0
	print 0, 0.0, a/2, a/2, 0.0, 0.0, 0.0
" > $dir/tests/G1/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/G1/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/G1/pvsv.dat
awk -v a0=$G1mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/16, pum, a;
				     }' data.stress >> $dir/tests/G1/pvsv.dat 

echo "#v, a" > $dir/tests/G1/evsv.dat
awk -v a0=$G1mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/16, $6;
				      }' data.energy >> $dir/tests/G1/evsv.dat 


sed -i '1d' $dir/tests/G1/evsv.dat
line=`python $BMFit $dir/tests/D03/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/G1/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	G1bulkp=`echo $line | awk '{print $1}'`
	G1vop=`echo $line | awk '{print $2}'`
	G1pote=`echo $line | awk '{print $3}'`
	G1eqp=`python -c "print (2.*$G1vop)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	G1bulkp=0
	G1vop=0
	G1pote=0
	G1eqp=0
fi

echo -e "${prp}E-V curve...${non}"

G1mmz=`echo "2*$G1eqp" | bc -l`
rm -f $dir/tests/G1/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)*$G1mmz*$G1mmz*$G1mmz
	a = math.pow(a,1./3)
	
	print '#N', 16, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, a/2, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, a/2, 0.0, 0.0, 0.0, 0.0
	print 0, 0.0, 0.0, a/2, 0.0, 0.0, 0.0
	print 1, a/2, a/2, a/2, 0.0, 0.0, 0.0
	print 0, a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 1, 3*a/4, a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, a/4, 3*a/4, 0.0, 0.0, 0.0
	print 1, a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 0, 3*a/4, 3*a/4, 3*a/4, 0.0, 0.0, 0.0
	print 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, a/2, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, 0.0, a/2, 0.0, 0.0, 0.0
	print 0, 0.0, a/2, a/2, 0.0, 0.0, 0.0

" > $dir/tests/G1/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/G1/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/G1/pvsv.dat
awk -v a0=$G1mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/16, pum, a;
				     }' data.stress >> $dir/tests/G1/pvsv.dat 

echo "#v, a" > $dir/tests/G1/evsv.dat
awk -v a0=$G1mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/16.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/G1/evsv.dat 
fi # lammps_evsv_flag G1

awk -v Voo=$G1vop '{print $1/Voo, $2, $3}' $dir/tests/G1/pvsv.dat > tmp; mv tmp $dir/tests/G1/pvsv.dat
for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/G1/pvsv.dat | sort -k2,2 -g | head -1`
	G1PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${G1PLAT[0]} $G1eqp"
		#G1eqp=`echo $line | awk '{print $3}'`
		#G1vop=`echo $line | awk '{print $1}'`
		#G1pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
G1PLAT["0"]=$G1eqp

# ---------------------- elastic constants for G1 ---------------------------------
echo -e "${prp}Elastic constants:${non}"
echo "# pressure, c11, c12, c44" > $dir/tests/G1/C_VS_P.dat
for PRESS in $PRESSEQ; do 
echo "$PRESS GPa..."
G1lmp=`echo 2*${G1PLAT["$PRESS"]} | bc -l`
for jj in $(seq 1 7); do

printf "\t $jj "

if [ "$lammps_elcon_g1" == "TRUE" ]; then

P=`echo $PRESS*10000 | bc -l`
case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac

cat > $dir/tests/G1/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms
variable boxa equal $G1lmp
variable boxb equal $G1lmp
variable boxc equal $G1lmp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $G1lmp  &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 
		
region box prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box

#create atoms
create_box      2 box
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Ti"]} basis 15 ${idx["Ti"]} basis 16 ${idx["Ti"]}
 
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/G1/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)

variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/G1/elcon.lin loop 
!!

$lammps < $dir/tests/G1/elcon.lin > $dir/tests/G1/elcon.out
sed -i '1d' $dir/tests/G1/stresses$jj.dat

else # else stress strain-curves with meamz

cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 1

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/G1/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
strain=0.002
rm -f $dir/tests/G1/elcon$jj.conf
for i in $(seq -$DELPT $DELPT); do
	G1mmz=`echo 2*${G1PLAT["$PRESS"]} | bc -l`
	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del G1 $G1mmz conf >> $dir/tests/G1/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/G1/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=-pc*$5; getline	
					syy=-pc*$5; getline	
					szz=-pc*$5; getline	
					sxy=-pc*$5; getline	
					syz=-pc*$5; getline	
					szx=-pc*$5; 	
					print del, sxx, syy, szz, syz, szx, sxy
				     }' data.stress >> $dir/tests/G1/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""

rm -f $dir/tests/G1/EC_fits.dat; touch $dir/tests/G1/EC_fits.dat

# first row fits
awk '{print $1, $2}' $dir/tests/G1/stresses1.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/G1/stresses1.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/G1/stresses1.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/G1/stresses2.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/G1/stresses2.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/G1/stresses2.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/G1/stresses3.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/G1/stresses3.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/G1/stresses3.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/G1/stresses4.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/G1/stresses4.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/G1/stresses4.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/G1/stresses5.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/G1/stresses6.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/G1/stresses7.dat > tmp; python $fit tmp >> $dir/tests/G1/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/G1/EC_fits.dat` )
G1C11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
G1C12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
G1C13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
G1C22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
G1C23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
G1C33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
G1C44i=${coup[12]}
G1C55i=${coup[13]}
G1C66i=${coup[14]}

G1C11p=`python -c "print ($G1C11i+$G1C22i+$G1C33i)/3"`
G1C12p=`python -c "print ($G1C12i+$G1C23i+$G1C13i)/3"`
G1C44p=`python -c "print ($G1C44i+$G1C55i+$G1C66i)/3"`

echo $PRESS $G1C11p $G1C12p $G1C44p >> $dir/tests/G1/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	G1C11=`printf '%3.f' $G1C11p`
	G1C12=`printf '%3.f' $G1C12p`
	G1C44=`printf '%3.f' $G1C44p`

	C11G1d=`printf '%3.f' $C11G1d`
	C12G1d=`printf '%3.f' $C12G1d`
	C44G1d=`printf '%3.f' $C44G1d`
fi
done

#----------------------------- bct deformation for G1 -----------------------------------------
echo "BCT deformation..."
G1lmp=`echo "2*$G1eqp" | bc -l`
cat > $dir/tests/G1/bct.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR G1 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmin equal 0.80 
variable	dmax equal 1.60 
variable	jmax equal $nevpt

variable boxa equal $G1lmp
variable boxb equal $G1lmp
variable boxc equal $G1lmp
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $G1lmp  &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.50 0.00 0.00 basis 0.00 0.50 0.00 &
		basis 0.00 0.00 0.50 basis 0.50 0.50 0.50 &
		basis 0.25 0.25 0.25 basis 0.75 0.25 0.25 &
		basis 0.25 0.75 0.25 basis 0.25 0.25 0.75 &
		basis 0.75 0.75 0.25 basis 0.75 0.25 0.75 &
		basis 0.25 0.75 0.75 basis 0.75 0.75 0.75 &
		basis 0.00 0.00 0.00 basis 0.00 0.50 0.50 &
		basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]} &
		basis 5 ${idx["Ti"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} &
		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Ti"]} &
		basis 13 ${idx["Nb"]} basis 14 ${idx["Ti"]} basis 15 ${idx["Ti"]} basis 16 ${idx["Ti"]}
		
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable E0 equal \${Epa}
variable DE equal v_Epa-\${E0}
variable coa equal lz/\${boxa}

fix P all print 1 "\${coa} \${DE}" file $dir/tests/G1/bct.dat screen no title "# V/atom | Energy"

variable zd equal \${boxc}*\${dmin}
variable zf equal \${boxc}*\${dmax}
change_box all z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/G1/bct.in loop
!!
$lammps < $dir/tests/G1/bct.in > $dir/tests/G1/bct.out

## PHONONS FOR G1
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

G1mmz=`echo "2*$G1lat" | bc -l`
cat > $dir/tests/G1/OPT.POSCAR << -
G1 Ti3Nb
$G1mmz
0.25	-0.25	0.75
0.25	0.75	-0.25
-0.75	-0.25	-0.25
Ti Nb
3 1
Direct
0.500000000000	0.000000000000	0.000000000000
0.000000000000	0.500000000000	0.000000000000
0.000000000000	0.000000000000	0.500000000000
0.000000000000	0.000000000000	0.000000000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/G1/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import FIRE
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/G1/OPT.POSCAR')

atoms.set_calculator(calc)

opt = FIRE(atoms=atoms)
opt.run()

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

G = [0, 0, 0]
X = [1./2, 0, 0]
M = [1./2, 1./2, 0]
R = [1./2, 1./2, 1./2]
Rp = [1./2, -1./2, -1./2]

point_names = ['', '\$\Gamma\$', '', '', '\$\Gamma\$']
path = [R, G, M, X, G]
dirs = ['\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi0]\$', '\$[0\\\bar{\\\xi}0]\$', '\$[\\\xi00]\$']

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_band = np.max([max_band, max_this])
    min_band = np.min([min_band, min_this])
    plt.plot(q, omega_n, 'g-', lw=2)

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('G1 Ti\$_{3}\$Nb')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='lightgreen', edgecolor='g', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/G1_phonons.png')
ph.clean
!

python $dir/tests/G1/phonons.py > $dir/tests/G1/phonon_log
rm -f *.pckl
fi
#--------------------------------------------------------------------------------------------



#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SQS7525 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if [ "$SQS_FLAG" == "TRUE" ]; then

echo -e "${yel} \t SQS 75/25 ${non}"
echo -e "${prp}equilibria...${non}"

if [ "$lammps_SQS_flag" == "TRUE" ]; then

$sqsgen 75 $SQS7525lat lmp > $dir/tests/SQS7525/SQS.dat

cat > $dir/tests/SQS7525/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM BCC-SQS LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

box tilt large
read_data $dir/tests/SQS7525/SQS.dat
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz xy xz yz cellalpha cellbeta cellgamma
thermo 1
dump d all xyz 1 /n/jww-1/ehemann.2/testingScripts/binary/sqs.xyz

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms
variable	XY equal xy
variable	XZ equal xz
variable	YZ equal yz
variable	XL equal xlo
variable	XH equal xhi
variable	YL equal ylo
variable	YH equal yhi
variable	ZL equal zlo
variable	ZH equal zhi

dump 1 all xyz 1 $dir/tests/SQS7525/SQS_RELAXED.xyz
run 0

print '\${a} \${coa} \${vat} \${spe} \${XL} \${XH} \${YL} \${YH} \${ZL} \${ZH} \${XY} \${XZ} \${YZ}'
!!

$lammps < $dir/tests/SQS7525/eq.in > $dir/tests/SQS7525/eq.out

SQS7525eqp=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $1/sqrt(5)}'`
SQS7525vop=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $3}'`
SQS7525pote=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $4}'`

# now get relaxed triclinic box parameters
xlSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $5}'`
xhSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $6}'`
ylSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $7}'`
yhSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $8}'`
zlSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $9}'`
zhSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $10}'`
xySQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $11}'`
xzSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $12}'`
yzSQS7525=`tail -1 $dir/tests/SQS7525/eq.out | awk '{print $13}'`

# now write relaxed data file to be read in by e-v calculator
cat > $dir/tests/SQS7525/SQS_RELAXED.coords << -+
## comment line
16 atoms
2 atom types
$xlSQS7525 $xhSQS7525 xlo xhi
$ylSQS7525 $yhSQS7525 ylo yhi
$zlSQS7525 $zhSQS7525 zlo zhi
$xySQS7525 $xzSQS7525 $yzSQS7525 xy xz yz
Atoms

-+

tail -16 $dir/tests/SQS7525/SQS_RELAXED.xyz | awk '{print NR, $1, $2, $3, $4}' >> $dir/tests/SQS7525/SQS_RELAXED.coords

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for SQS -----------------------------------------
cat > $dir/tests/SQS7525/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCC-SQS LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

box tilt large
read_data $dir/tests/SQS7525/SQS_RELAXED.coords


mass 1 $mass1
mass 2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cella cellb cellc cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}
variable tmp equal xy
variable XY0 equal \${tmp}
variable tmp equal xz
variable XZ0 equal \${tmp}
variable tmp equal yz
variable YZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/SQS7525/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/SQS7525/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)
variable XYI equal \${XY0}*(1-v_dmax/2)
variable XZI equal \${XZ0}*(1-v_dmax/2)
variable YZI equal \${YZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)
variable XYF equal \${XY0}*(1+v_dmax/2)
variable XZF equal \${XZ0}*(1+v_dmax/2)
variable YZF equal \${YZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} xy final \${XYI} xz final \${XZI} yz final \${YZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} xy final \${XYF} xz final \${XZF} yz final \${YZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/SQS7525/evsv.in loop
!!
$lammps < $dir/tests/SQS7525/evsv.in > $dir/tests/SQS7525/evsv.out

sed -i '1d' $dir/tests/SQS7525/evsv.dat
line=`python $BMFit $dir/tests/SQS7525/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/SQS7525/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	SQS7525bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	SQS7525bulkp=0
fi
sed -i '1d' $dir/tests/SQS7525/pvsv.dat


else # SQS lammps flag


## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/SQS7525/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

SQS7525mmz=`echo "$SQS7525lat" | bc -l`
rm -f $dir/tests/SQS7525/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$SQS7525mmz
	
	print '#N', 16, 2
	print '##'
	print '#X', a, -2.005055*a, 0.012647*a
	print '#Y', 0.012647*a, -2.005055*a, a
	print '#Z', -2.00209*a, 0.0116144*a, -2.00209*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, -1.4409*a, -0.559058*a, -1.4409*a, 0.0, 0.0, 0.0
	print 0, -1.48021*a, -1.47393*a, -1.48021*a, 0.0, 0.0, 0.0
	print 0, -0.92789*a, -1.05406*a, -0.92789*a, 0.0, 0.0, 0.0
	print 0, -0.526474*a, -0.49512*a, -0.526474*a, 0.0, 0.0, 0.0
	print 0, -1.97366*a, -0.0238696*a, -1.97366*a, 0.0, 0.0, 0.0
	print 0, -1.00347*a, -2.00829*a, -1.00347*a, 0.0, 0.0, 0.0
	print 0, -0.448366*a, -1.55546*a, -0.448366*a, 0.0, 0.0, 0.0
	print 0, -0.001167*a, -0.984167*a, -0.001167*a, 0.0, 0.0, 0.0
	print 0, -1.05055*a, -2.94944*a, -1.04044*a, 0.0, 0.0, 0.0
	print 0, -0.491212*a, -2.50102*a, -0.491212*a, 0.0, 0.0, 0.0
	print 0, -0.549298*a, -3.44459*a, -0.549298*a, 0.0, 0.0, 0.0
	print 0, 0.00407678*a, -3.0102*a, 0.00407678*a, 0.0, 0.0, 0.0
	print 1, -0.0402315*a, -3.95859*a, -0.0402315*a, 0.0, 0.0, 0.0
	print 1, 0.0223491*a, -2.00314*a, 0.0223491*a, 0.0, 0.0, 0.0
	print 1, 0.477775*a, -2.47544*a, 0.477775*a, 0.0, 0.0, 0.0
	print 1, 0.512644*a, -3.5214*a, 0.512644*a, 0.0, 0.0, 0.0

" > $dir/tests/SQS7525/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/SQS7525/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/SQS7525/pvsv.dat
awk -v a0=$SQS7525mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = 7.9155*a0*a0*a0*vol
					   print vol/16, pum, a;
				     }' data.stress >> $dir/tests/SQS7525/pvsv.dat 

echo "#v, a" > $dir/tests/SQS7525/evsv.dat
awk -v a0=$SQSmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*7.9155;
					   print vol/16, $6;
				      }' data.energy >> $dir/tests/SQS7525/evsv.dat 


sed -i '1d' $dir/tests/SQS7525/evsv.dat
line=`python $BMFit $dir/tests/SQS7525/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/SQS7525/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	SQS7525bulkp=`echo $line | awk '{print $1}'`
	SQS7525vop=`echo $line | awk '{print $2}'`
	SQS7525pote=`echo $line | awk '{print $3}'`
	SQS7525eqp=`python -c "print (16.*$SQS7525vop/7.9155)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	SQS7525bulkp=0
	SQS7525vop=0
	SQS7525pote=0
	SQS7525eqp=0
fi

echo -e "${prp}E-V curve...${non}"

SQSmmz=`echo "$SQS7525eqp" | bc -l`
rm -f $dir/tests/SQS7525/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$SQSmmz
	
	print '#N', 16, 2
	print '##'
	print '#X', a, -2.005055*a, 0.012647*a
	print '#Y', 0.012647*a, -2.005055*a, a
	print '#Z', -2.00209*a, 0.0116144*a, -2.00209*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, -1.4409*a, -0.559058*a, -1.4409*a, 0.0, 0.0, 0.0
	print 0, -1.48021*a, -1.47393*a, -1.48021*a, 0.0, 0.0, 0.0
	print 0, -0.92789*a, -1.05406*a, -0.92789*a, 0.0, 0.0, 0.0
	print 0, -0.526474*a, -0.49512*a, -0.526474*a, 0.0, 0.0, 0.0
	print 0, -1.97366*a, -0.0238696*a, -1.97366*a, 0.0, 0.0, 0.0
	print 0, -1.00347*a, -2.00829*a, -1.00347*a, 0.0, 0.0, 0.0
	print 0, -0.448366*a, -1.55546*a, -0.448366*a, 0.0, 0.0, 0.0
	print 0, -0.001167*a, -0.984167*a, -0.001167*a, 0.0, 0.0, 0.0
	print 0, -1.05055*a, -2.94944*a, -1.04044*a, 0.0, 0.0, 0.0
	print 0, -0.491212*a, -2.50102*a, -0.491212*a, 0.0, 0.0, 0.0
	print 0, -0.549298*a, -3.44459*a, -0.549298*a, 0.0, 0.0, 0.0
	print 0, 0.00407678*a, -3.0102*a, 0.00407678*a, 0.0, 0.0, 0.0
	print 1, -0.0402315*a, -3.95859*a, -0.0402315*a, 0.0, 0.0, 0.0
	print 1, 0.0223491*a, -2.00314*a, 0.0223491*a, 0.0, 0.0, 0.0
	print 1, 0.477775*a, -2.47544*a, 0.477775*a, 0.0, 0.0, 0.0
	print 1, 0.512644*a, -3.5214*a, 0.512644*a, 0.0, 0.0, 0.0

" > $dir/tests/SQS7525/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/SQS7525/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/SQS7525/pvsv.dat
awk -v a0=$SQSmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = vol*a0*a0*a0*7.9155;
					   print vol/16, pum, a;
				     }' data.stress >> $dir/tests/SQS7525/pvsv.dat 

echo "#v, a" > $dir/tests/SQS7525/evsv.dat
awk -v a0=$SQS7525mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*7.9155/16.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/SQS7525/evsv.dat 
fi # lammps_evsv_flag SQS

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/SQS7525/pvsv.dat | sort -k2,2 -g | head -1`
	SQS7525PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${SQS7525PLAT[0]} $SQS7525eqp"
	#	SQS7525eqp=`echo $line | awk '{print $3}'`
	#	SQS7525vop=`echo $line | awk '{print $1}'`
		#SQS7525pote=`echo $line | awk '{print 1000*$4}'`
	fi
done

awk -v Voo=$SQS7525vop '{print $1/Voo, $2, $3}' $dir/tests/SQS7525/pvsv.dat > tmp; mv tmp $dir/tests/SQS7525/pvsv.dat
SQS7525PLAT["0"]=$SQS7525eqp


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SQS5050 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

echo -e "${yel} \t SQS 50/50 ${non}"
echo -e "${prp}equilibria...${non}"

if [ "$lammps_SQS_flag" == "TRUE" ]; then

$sqsgen 75 $SQS5050lat lmp > $dir/tests/SQS5050/SQS.dat

cat > $dir/tests/SQS5050/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM BCC-SQS LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

box tilt large
read_data $dir/tests/SQS5050/SQS.dat
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz xy xz yz cellalpha cellbeta cellgamma
thermo 1
dump d all xyz 1 /n/jww-1/ehemann.2/testingScripts/binary/sqs.xyz

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms
variable	XY equal xy
variable	XZ equal xz
variable	YZ equal yz
variable	XL equal xlo
variable	XH equal xhi
variable	YL equal ylo
variable	YH equal yhi
variable	ZL equal zlo
variable	ZH equal zhi

dump 1 all xyz 1 $dir/tests/SQS5050/SQS_RELAXED.xyz
run 0

print '\${a} \${coa} \${vat} \${spe} \${XL} \${XH} \${YL} \${YH} \${ZL} \${ZH} \${XY} \${XZ} \${YZ}'
!!

$lammps < $dir/tests/SQS5050/eq.in > $dir/tests/SQS5050/eq.out

SQS5050eqp=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $1}'`
SQS5050vop=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $3}'`
SQS5050pote=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $4}'`

# now get relaxed triclinic box parameters
xlSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $5}'`
xhSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $6}'`
ylSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $7}'`
yhSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $8}'`
zlSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $9}'`
zhSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $10}'`
xySQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $11}'`
xzSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $12}'`
yzSQS5050=`tail -1 $dir/tests/SQS5050/eq.out | awk '{print $13}'`

# now write relaxed data file to be read in by e-v calculator
cat > $dir/tests/SQS5050/SQS_RELAXED.coords << -+
## comment line
16 atoms
2 atom types
$xlSQS5050 $xhSQS5050 xlo xhi
$ylSQS5050 $yhSQS5050 ylo yhi
$zlSQS5050 $zhSQS5050 zlo zhi
$xySQS5050 $xzSQS5050 $yzSQS5050 xy xz yz
Atoms

-+

tail -16 $dir/tests/SQS5050/SQS_RELAXED.xyz | awk '{print NR, $1, $2, $3, $4}' >> $dir/tests/SQS5050/SQS_RELAXED.coords

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for SQS -----------------------------------------
cat > $dir/tests/SQS5050/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCC-SQS LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

box tilt large
read_data $dir/tests/SQS5050/SQS_RELAXED.coords


mass 1 $mass1
mass 2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cella cellb cellc cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}
variable tmp equal xy
variable XY0 equal \${tmp}
variable tmp equal xz
variable XZ0 equal \${tmp}
variable tmp equal yz
variable YZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/SQS5050/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/SQS5050/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)
variable XYI equal \${XY0}*(1-v_dmax/2)
variable XZI equal \${XZ0}*(1-v_dmax/2)
variable YZI equal \${YZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)
variable XYF equal \${XY0}*(1+v_dmax/2)
variable XZF equal \${XZ0}*(1+v_dmax/2)
variable YZF equal \${YZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} xy final \${XYI} xz final \${XZI} yz final \${YZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} xy final \${XYF} xz final \${XZF} yz final \${YZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/SQS5050/evsv.in loop
!!
$lammps < $dir/tests/SQS5050/evsv.in > $dir/tests/SQS5050/evsv.out

sed -i '1d' $dir/tests/SQS5050/evsv.dat
line=`python $BMFit $dir/tests/SQS5050/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/SQS5050/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	SQS5050bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	SQS5050bulkp=0
fi
sed -i '1d' $dir/tests/SQS5050/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/SQS5050/pvsv.dat | sort -k2,2 -g | head -1`
	SQS5050PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${SQS5050PLAT[0]} $SQS5050eqp"
	#	SQS5050eqp=`echo $line | awk '{print $3}'`
	#	SQS5050vop=`echo $line | awk '{print $1}'`
		#SQS5050pote=`echo $line | awk '{print 1000*$4}'`
	fi
done

awk -v Voo=$SQS5050vop '{print $1/Voo, $2, $3}' $dir/tests/SQS5050/pvsv.dat > tmp; mv tmp $dir/tests/SQS5050/pvsv.dat
SQS5050PLAT["0"]=$SQS5050eqp

fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SQS2575 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

echo -e "${yel} \t SQS 75/25 ${non}"
echo -e "${prp}equilibria...${non}"

if [ "$lammps_SQS_flag" == "TRUE" ]; then

$sqsgen 75 $SQS2575lat lmp > $dir/tests/SQS2575/SQS.dat

cat > $dir/tests/SQS2575/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM BCC-SQS LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

box tilt large
read_data $dir/tests/SQS2575/SQS.dat
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz xy xz yz cellalpha cellbeta cellgamma
thermo 1
dump d all xyz 1 /n/jww-1/ehemann.2/testingScripts/binary/sqs.xyz

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
min_style	cg
minimize	$BETA_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms
variable	XY equal xy
variable	XZ equal xz
variable	YZ equal yz
variable	XL equal xlo
variable	XH equal xhi
variable	YL equal ylo
variable	YH equal yhi
variable	ZL equal zlo
variable	ZH equal zhi

dump 1 all xyz 1 $dir/tests/SQS2575/SQS_RELAXED.xyz
run 0

print '\${a} \${coa} \${vat} \${spe} \${XL} \${XH} \${YL} \${YH} \${ZL} \${ZH} \${XY} \${XZ} \${YZ}'
!!

$lammps < $dir/tests/SQS2575/eq.in > $dir/tests/SQS2575/eq.out

SQS2575eqp=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $1}'`
SQS2575vop=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $3}'`
SQS2575pote=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $4}'`

# now get relaxed triclinic box parameters
xlSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $5}'`
xhSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $6}'`
ylSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $7}'`
yhSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $8}'`
zlSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $9}'`
zhSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $10}'`
xySQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $11}'`
xzSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $12}'`
yzSQS2575=`tail -1 $dir/tests/SQS2575/eq.out | awk '{print $13}'`

# now write relaxed data file to be read in by e-v calculator
cat > $dir/tests/SQS2575/SQS_RELAXED.coords << -+
## comment line
16 atoms
2 atom types
$xlSQS2575 $xhSQS2575 xlo xhi
$ylSQS2575 $yhSQS2575 ylo yhi
$zlSQS2575 $zhSQS2575 zlo zhi
$xySQS2575 $xzSQS2575 $yzSQS2575 xy xz yz
Atoms

-+

tail -16 $dir/tests/SQS2575/SQS_RELAXED.xyz | awk '{print NR, $1, $2, $3, $4}' >> $dir/tests/SQS2575/SQS_RELAXED.coords

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for SQS -----------------------------------------
cat > $dir/tests/SQS2575/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR BCC-SQS LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

box tilt large
read_data $dir/tests/SQS2575/SQS_RELAXED.coords


mass 1 $mass1
mass 2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cella cellb cellc cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}
variable tmp equal xy
variable XY0 equal \${tmp}
variable tmp equal xz
variable XZ0 equal \${tmp}
variable tmp equal yz
variable YZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/SQS2575/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/SQS2575/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)
variable XYI equal \${XY0}*(1-v_dmax/2)
variable XZI equal \${XZ0}*(1-v_dmax/2)
variable YZI equal \${YZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)
variable XYF equal \${XY0}*(1+v_dmax/2)
variable XZF equal \${XZ0}*(1+v_dmax/2)
variable YZF equal \${YZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} xy final \${XYI} xz final \${XZI} yz final \${YZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} xy final \${XYF} xz final \${XZF} yz final \${YZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/SQS2575/evsv.in loop
!!
$lammps < $dir/tests/SQS2575/evsv.in > $dir/tests/SQS2575/evsv.out

sed -i '1d' $dir/tests/SQS2575/evsv.dat
line=`python $BMFit $dir/tests/SQS2575/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/SQS2575/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	SQS2575bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	SQS2575bulkp=0
fi
sed -i '1d' $dir/tests/SQS2575/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/SQS2575/pvsv.dat | sort -k2,2 -g | head -1`
	SQS2575PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${SQS2575PLAT[0]} $SQS2575eqp"
	#	SQS2575eqp=`echo $line | awk '{print $3}'`
	#	SQS2575vop=`echo $line | awk '{print $1}'`
		#SQS2575pote=`echo $line | awk '{print 1000*$4}'`
	fi
done

awk -v Voo=$SQS2575vop '{print $1/Voo, $2, $3}' $dir/tests/SQS2575/pvsv.dat > tmp; mv tmp $dir/tests/SQS2575/pvsv.dat
SQS2575PLAT["0"]=$SQS2575eqp

fi
fi # SQS flag

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< L60 Ti3Nb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t L60 Ti3Nb ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/L60Ti3Nb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM L60Ti3Nb LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

variable SIZE equal $L60Ti3Nblat
variable YIZE equal 1.403836667*$L60Ti3Nblat

lattice		custom $L60Ti3Nblat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.403836667 0.0 &
		a3 0.0 0.0 1.403836667 &
		basis 0.0	0.5	0.50 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.00 &
		basis 0.5	0.0	0.50
		
region		mybox block 0 \${SIZE} 0 \${YIZE} 0 \${YIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Nb"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} 
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax iso 0.0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

run 0

print '\${a} \${vat} \${spe}'
!!

$lammps < $dir/tests/L60Ti3Nb/eq.in > $dir/tests/L60Ti3Nb/eq.out

L60Ti3Nbeqp=`tail -1 $dir/tests/L60Ti3Nb/eq.out | awk '{print $1}'`
L60Ti3Nbvop=`tail -1 $dir/tests/L60Ti3Nb/eq.out | awk '{print $2}'`
L60Ti3Nbpote=`tail -1 $dir/tests/L60Ti3Nb/eq.out | awk '{print $3}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for L60Ti3Nb ------------------------------------------

cat > $dir/tests/L60Ti3Nb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR L60Ti3Nb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable SIZE equal $L60Ti3Nblat
variable YIZE equal 1.403836667*$L60Ti3Nblat

lattice		custom $L60Ti3Nblat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.403836667 0.0 &
		a3 0.0 0.0 1.403836667 &
		basis 0.0	0.5	0.50 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.00 &
		basis 0.5	0.0	0.50
		
region		mybox block 0 \${SIZE} 0 \${YIZE} 0 \${YIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Nb"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} 

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/L60Ti3Nb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/L60Ti3Nb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/L60Ti3Nb/evsv.in loop
!!
$lammps < $dir/tests/L60Ti3Nb/evsv.in > $dir/tests/L60Ti3Nb/evsv.out

sed -i '1d' $dir/tests/L60Ti3Nb/evsv.dat
line=`python $BMFit $dir/tests/L60Ti3Nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/L60Ti3Nb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	L60Ti3Nbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	L60Ti3Nbbulkp=0
fi

else # L60Ti3Nb lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/L60Ti3Nb/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

L60Ti3Nbmmz=`echo "$L60Ti3Nblat" | bc -l`
echo $L60Ti3Nbmmz
rm -f $dir/tests/L60Ti3Nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$L60Ti3Nbmmz
	
	print '#N', 4, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, 1.403836667*a, 0
	print '#Z', 0, 0, 1.403836667*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0.5,	0.5,	0,	0,	0
	print 0,	0,	0,	0,	0,	0,	0
	print 0,	0.5,	0.5,	0,	0,	0,	0
	print 0,	0.5,	0.0,	0.5,	0,	0,	0
" > $dir/tests/L60Ti3Nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/L60Ti3Nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/L60Ti3Nb/pvsv.dat
awk -v a0=$L60Ti3Nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol;
					   print vol/4, pum, a;
				     }' data.stress >> $dir/tests/L60Ti3Nb/pvsv.dat 

echo "#v, a" > $dir/tests/L60Ti3Nb/evsv.dat
awk -v a0=$L60Ti3Nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/4, $6;
				      }' data.energy >> $dir/tests/L60Ti3Nb/evsv.dat 


sed -i '1d' $dir/tests/L60Ti3Nb/evsv.dat
line=`python $BMFit $dir/tests/L60Ti3Nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/L60Ti3Nb/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	L60Ti3Nbbulkp=`echo $line | awk '{print $1}'`
	L60Ti3Nbvop=`echo $line | awk '{print $2}'`
	L60Ti3Nbpote=`echo $line | awk '{print $3}'`
	L60Ti3Nbeqp=`python -c "print (4.*$L60Ti3Nbvop/1.97075739)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	L60Ti3Nbbulkp=0
	L60Ti3Nbvop=0
	L60Ti3Nbpote=0
	L60Ti3Nbeqp=0
fi

echo -e "${prp}E-V curve...${non}"

L60Ti3Nbmmz=`echo "$L60Ti3Nbeqp" | bc -l`
rm -f $dir/tests/L60Ti3Nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $L60Ti3Nbmmz*math.pow(a,1./3)
	
	print '#N', 4, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, 1.403836667*a, 0
	print '#Z', 0, 0, 1.403836667*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0.5,	0.5,	0,	0,	0
	print 0,	0,	0,	0,	0,	0,	0
	print 0,	0.5,	0.5,	0,	0,	0,	0
	print 0,	0.5,	0.0,	0.5,	0,	0,	0

" > $dir/tests/L60Ti3Nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/L60Ti3Nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/L60Ti3Nb/pvsv.dat
awk -v a0=$L60Ti3Nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/4, pum, a;
				     }' data.stress >> $dir/tests/L60Ti3Nb/pvsv.dat 

echo "#v, a" > $dir/tests/L60Ti3Nb/evsv.dat
awk -v a0=$L60Ti3Nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/8.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/L60Ti3Nb/evsv.dat 
fi # lammps_evsv_flag L60Ti3Nb

sed -i '1d' $dir/tests/L60Ti3Nb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/L60Ti3Nb/pvsv.dat | sort -k2,2 -g | head -1`
	L60Ti3NbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${L60Ti3NbPLAT[0]} $L60Ti3Nbeqp"
		#L60Ti3Nbeqp=`echo $line | awk '{print $3}'`
		#L60Ti3Nbvop=`echo $line | awk '{print $1}'`
		#L60Ti3Nbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
L60Ti3NbPLAT["0"]=$L60Ti3Nbeqp


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< L60 TiNb3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t L60 TiNb3 ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/L60TiNb3/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM L60TiNb3 LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

variable SIZE equal $L60TiNb3lat
variable YIZE equal 1.403836667*$L60TiNb3lat

lattice		custom $L60TiNb3lat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.403836667 0.0 &
		a3 0.0 0.0 1.403836667 &
		basis 0.0	0.5	0.50 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.00 &
		basis 0.5	0.0	0.50
		
region		mybox block 0 \${SIZE} 0 \${YIZE} 0 \${YIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} 
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax iso 0.0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

run 0

print '\${a} \${vat} \${spe}'
!!

$lammps < $dir/tests/L60TiNb3/eq.in > $dir/tests/L60TiNb3/eq.out

L60TiNb3eqp=`tail -1 $dir/tests/L60TiNb3/eq.out | awk '{print $1}'`
L60TiNb3vop=`tail -1 $dir/tests/L60TiNb3/eq.out | awk '{print $2}'`
L60TiNb3pote=`tail -1 $dir/tests/L60TiNb3/eq.out | awk '{print $3}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for L60TiNb3 ------------------------------------------

cat > $dir/tests/L60TiNb3/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR L60TiNb3 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable SIZE equal $L60TiNb3lat
variable YIZE equal 1.403836667*$L60TiNb3lat

lattice		custom $L60TiNb3lat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.403836667 0.0 &
		a3 0.0 0.0 1.403836667 &
		basis 0.0	0.5	0.50 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.00 &
		basis 0.5	0.0	0.50
		
region		mybox block 0 \${SIZE} 0 \${YIZE} 0 \${YIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} 

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/L60TiNb3/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/L60TiNb3/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/L60TiNb3/evsv.in loop
!!
$lammps < $dir/tests/L60TiNb3/evsv.in > $dir/tests/L60TiNb3/evsv.out

sed -i '1d' $dir/tests/L60TiNb3/evsv.dat
line=`python $BMFit $dir/tests/L60TiNb3/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/L60TiNb3/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	L60TiNb3bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	L60TiNb3bulkp=0
fi

else # L60TiNb3 lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/L60TiNb3/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

L60TiNb3mmz=`echo "$L60TiNb3lat" | bc -l`
echo $L60TiNb3mmz
rm -f $dir/tests/L60TiNb3/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$L60TiNb3mmz
	
	print '#N', 4, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, 1.403836667*a, 0
	print '#Z', 0, 0, 1.403836667*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0,	0,	0.5,	0.5,	0,	0,	0
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.5,	0.5,	0,	0,	0,	0
	print 1,	0.5,	0.0,	0.5,	0,	0,	0
" > $dir/tests/L60TiNb3/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/L60TiNb3/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/L60TiNb3/pvsv.dat
awk -v a0=$L60TiNb3mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol;
					   print vol/4, pum, a;
				     }' data.stress >> $dir/tests/L60TiNb3/pvsv.dat 

echo "#v, a" > $dir/tests/L60TiNb3/evsv.dat
awk -v a0=$L60TiNb3mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/4, $6;
				      }' data.energy >> $dir/tests/L60TiNb3/evsv.dat 


sed -i '1d' $dir/tests/L60TiNb3/evsv.dat
line=`python $BMFit $dir/tests/L60TiNb3/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/L60TiNb3/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	L60TiNb3bulkp=`echo $line | awk '{print $1}'`
	L60TiNb3vop=`echo $line | awk '{print $2}'`
	L60TiNb3pote=`echo $line | awk '{print $3}'`
	L60TiNb3eqp=`python -c "print (4.*$L60TiNb3vop/1.97075739)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	L60TiNb3bulkp=0
	L60TiNb3vop=0
	L60TiNb3pote=0
	L60TiNb3eqp=0
fi

echo -e "${prp}E-V curve...${non}"

L60TiNb3mmz=`echo "$L60TiNb3eqp" | bc -l`
rm -f $dir/tests/L60TiNb3/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $L60TiNb3mmz*math.pow(a,1./3)
	
	print '#N', 4, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, 1.403836667*a, 0
	print '#Z', 0, 0, 1.403836667*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0,	0,	0.5,	0.5,	0,	0,	0
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.5,	0.5,	0,	0,	0,	0
	print 1,	0.5,	0.0,	0.5,	0,	0,	0

" > $dir/tests/L60TiNb3/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/L60TiNb3/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/L60TiNb3/pvsv.dat
awk -v a0=$L60TiNb3mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/4, pum, a;
				     }' data.stress >> $dir/tests/L60TiNb3/pvsv.dat 

echo "#v, a" > $dir/tests/L60TiNb3/evsv.dat
awk -v a0=$L60TiNb3mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/8.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/L60TiNb3/evsv.dat 
fi # lammps_evsv_flag L60TiNb3

sed -i '1d' $dir/tests/L60TiNb3/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/L60TiNb3/pvsv.dat | sort -k2,2 -g | head -1`
	L60TiNb3PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${L60TiNb3PLAT[0]} $L60TiNb3eqp"
		#L60TiNb3eqp=`echo $line | awk '{print $3}'`
		#L60TiNb3vop=`echo $line | awk '{print $1}'`
		#L60TiNb3pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
L60TiNb3PLAT["0"]=$L60TiNb3eqp


######################################################################
#	METASTABLE PHASES
######################################################################
echo -e "${grn}BEGINNING METASTABLE TESTS... ${non}"

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< APP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t alpha'' ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/APP/eq.in << !!
###############################################
#	CALCULATES EQUILIBRIUM APP LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

variable boxa equal $APPlat
variable boxb equal $APPlat*$APPboa
variable boxc equal $APPlat*$APPcoa

lattice		custom $APPlat  &
		a1 1.0 0.0 0.0 &
		a2 0.0 $APPboa 0.0 &
		a3 0.0 0.0000000 $APPcoa &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.00 &
		basis 0.50 0.10 0.50 basis 0.00 0.60 0.50 
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]} 
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 fixedpoint 0 0 0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx
variable	boa equal ly/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${boa} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/APP/eq.in > $dir/tests/APP/eq.out

APPeqp=`tail -1 $dir/tests/APP/eq.out | awk '{print $1}'`
APPboap=`tail -1 $dir/tests/APP/eq.out | awk '{print $2}'`
APPcoap=`tail -1 $dir/tests/APP/eq.out | awk '{print $3}'`
APPvop=`tail -1 $dir/tests/APP/eq.out | awk '{print $4}'`
APPpote=`tail -1 $dir/tests/APP/eq.out | awk '{print $5}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for APP ------------------------------------------

cat > $dir/tests/APP/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR APP LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $APPeqp
variable boxb equal $APPeqp*$APPboap
variable boxc equal $APPeqp*$APPcoap

lattice		custom $APPlat  &
		a1 1.0 0.0 0.0 &
		a2 0.0 $APPboap 0.0 &
		a3 0.0 0.0000000 $APPcoap &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.00 &
		basis 0.50 0.10 0.50 basis 0.00 0.60 0.50 
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]} 
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/APP/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/APP/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/APP/evsv.in loop
!!
$lammps < $dir/tests/APP/evsv.in > $dir/tests/APP/evsv.out

sed -i '1d' $dir/tests/APP/evsv.dat
line=`python $BMFit $dir/tests/APP/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/APP/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	APPbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	APPbulkp=0
fi

else # APP lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/APP/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

APPmmz=`echo "$APPlat" | bc -l`
rm -f $dir/tests/APP/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$APPmmz
	
	print '#N', 4, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, 1.4281437*a, 0
	print '#Z', 0, 0, 1.320359*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, 0.0, -0.00029*a, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, 0.726343*a, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, 0.143208*a, 0.66018*a, 0.0, 0.0, 0.0
	print 1, 0.0, 0.884535*a, 0.66018*a, 0.0, 0.0, 0.0
" > $dir/tests/APP/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/APP/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/APP/pvsv.dat
awk -v a0=$APPmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*1.88566*vol 
					   print vol/4, pum, a;
				     }' data.stress >> $dir/tests/APP/pvsv.dat 

echo "#v, a" > $dir/tests/APP/evsv.dat
awk -v a0=$APPmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*1.88566;
					   print vol/4, $6;
				      }' data.energy >> $dir/tests/APP/evsv.dat 


sed -i '1d' $dir/tests/APP/evsv.dat
line=`python $BMFit $dir/tests/APP/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/APP/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	APPbulkp=`echo $line | awk '{print $1}'`
	APPvop=`echo $line | awk '{print $2}'`
	APPpote=`echo $line | awk '{print $3}'`
	APPeqp=`python -c "print (4.*$APPvop/1.88566)**(1./3)"`

else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	APPbulkp=0
	APPvop=0
	APPpote=0
	APPeqp=0
fi
echo -e "${prp}E-V curve...${non}"

APPmmz=`echo "$APPeqp" | bc -l`
rm -f $dir/tests/APP/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$APPmmz
	
	print '#N', 4, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, 1.4281437*a, 0
	print '#Z', 0, 0, 1.320359*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, 0.0, -0.00029*a, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, 0.726343*a, 0.0, 0.0, 0.0, 0.0
	print 0, a/2, 0.143208*a, 0.66018*a, 0.0, 0.0, 0.0
	print 1, 0.0, 0.884535*a, 0.66018*a, 0.0, 0.0, 0.0

" > $dir/tests/APP/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/APP/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/APP/pvsv.dat
awk -v a0=$APPmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*1.88566*vol;
					   print vol/4, pum, a;
				     }' data.stress >> $dir/tests/APP/pvsv.dat 

echo "#v, a" > $dir/tests/APP/evsv.dat
awk -v a0=$APPmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*1.88566/4;
					   print vol, $6;
				     }' data.energy >> $dir/tests/APP/evsv.dat 
fi # lammps_evsv_flag APP

sed -i '1d' $dir/tests/APP/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/APP/pvsv.dat | sort -k2,2 -g | head -1`
	APPPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${APPPLAT[0]} $APPeqp"
		#APPeqp=`echo $line | awk '{print $3}'`
		#APPvop=`echo $line | awk '{print $1}'`
		#APPpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
APPPLAT["0"]=$APPeqp

awk -v Voo=$APPvop '{print $1/Voo, $2, $3}' $dir/tests/APP/pvsv.dat > tmp; mv tmp $dir/tests/APP/pvsv.dat

#-----------------------------------------------------------------------------------
# ---------------------- pressure dependence of lattice constants for app ----------
echo -e "${prp}Pressure dependence of lattice constants...${non}"

rm -f $dir/tests/APP/abcvp.dat; echo "#a, b/a, c/a" > $dir/tests/APP/abcvp.dat
for PRESS in $PRESSLAT; do
PRES=`echo "10000*$PRESS" | bc -l`
cat > $dir/tests/APP/abcvp.lin << __
## LAMMPS script
units metal
boundary p p p
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal ${APPPLAT["$PRESS"]}
variable boxb equal ${APPPLAT["$PRESS"]}*$APPboap
variable boxc equal ${APPPLAT["$PRESS"]}*$APPcoap

lattice		custom ${APPPLAT["$PRESS"]} &
		a1 1.0 0.0 0.0 &
		a2 0.0 $APPboap 0.0 &
		a3 0.0 0.0000000 $APPcoap &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.00 &
		basis 0.50 0.10 0.50 basis 0.00 0.60 0.50 
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]} 
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#variable bat equal $APPboap*0.5
#region frz1 sphere 0.0 0.0 0.0 0.1 units box
#region frz2 sphere 0.5 \${bat} 0.0 0.1 units box
#region frz union 2 frz1 frz2
#group frz region frz
#
#fix FREEZE freeze frz

thermo 1000
thermo_style custom step temp press pe ke etotal lx ly lz
timestep 0.001

variable a equal lx
variable boa equal ly/lx
variable coa equal lz/lx

fix P all press/berendsen aniso $PRES $PRES 1000.0
run 10000

dump D all custom 1 $dir/tests/APP/abcvp.coords xs ys zs

print "$PRESS \${a} \${boa} \${coa}"
__

$lammps < $dir/tests/APP/abcvp.lin > $dir/tests/APP/abcvp.out
tail -1 $dir/tests/APP/abcvp.out >> $dir/tests/APP/abcvp.dat

done

# ---------------------- elastic constants for APP ---------------------------------
echo -e "${prp}Elastic constants:${non}"
echo "# pressure, c11, c12, c44" > $dir/tests/APP/C_VS_P.dat
for PRESS in $PRESSEQ; do 
echo "$PRESS GPa..."
for jj in $(seq 1 7); do

printf "\t $jj "

if [ "$lammps_elcon_APP" == "TRUE" ]; then

P=`echo $PRESS*10000 | bc -l`
case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac

cat > $dir/tests/APP/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

# lattice and atoms
variable boxa equal ${APPPLAT["$PRESS"]}
variable boxb equal ${APPPLAT["$PRESS"]}*$APPboap
variable boxc equal ${APPPLAT["$PRESS"]}*$APPcoap
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom ${APPPLAT["$PRESS"]}  &
		a1 1.0 0.0 0.0 &
		a2 0.0 $APPboap 0.0 &
		a3 0.0 0.0000000 $APPcoap &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.00 &
		basis 0.50 0.10 0.50 basis 0.00 0.60 0.50 
		
region		box prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 box

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Nb"]}
 
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
#fix rel all box/relax iso $P vmax 0.001 fixedpoint 0 0 0
#min_style cg
#minimize 1e-20 1e-20 100000 100000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
timestep 0.001
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/APP/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)

variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/APP/elcon.lin loop 
!!

$lammps < $dir/tests/APP/elcon.lin > $dir/tests/APP/elcon.out
sed -i '1d' $dir/tests/APP/stresses$jj.dat

else # else stress strain-curves with meamz

cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/APP/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
strain=0.002
rm -f $dir/tests/APP/elcon$jj.conf
for i in $(seq -$DELPT $DELPT); do
	APPmmz=`echo ${APPPLAT["$PRESS"]} | bc -l`
	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del APP $APPmmz conf >> $dir/tests/APP/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/APP/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=-pc*$5; getline	
					syy=-pc*$5; getline	
					szz=-pc*$5; getline	
					sxy=-pc*$5; getline	
					syz=-pc*$5; getline	
					szx=-pc*$5; 	
					print del, sxx, syy, szz, syz, szx, sxy
				     }' data.stress >> $dir/tests/APP/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""

rm -f $dir/tests/APP/EC_fits.dat; touch $dir/tests/APP/EC_fits.dat

# first row fits
awk '{print $1, $2}' $dir/tests/APP/stresses1.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/APP/stresses1.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/APP/stresses1.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/APP/stresses2.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/APP/stresses2.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/APP/stresses2.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/APP/stresses3.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/APP/stresses3.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/APP/stresses3.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/APP/stresses4.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/APP/stresses4.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/APP/stresses4.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/APP/stresses5.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/APP/stresses6.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/APP/stresses7.dat > tmp; python $fit tmp >> $dir/tests/APP/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/APP/EC_fits.dat` )
APPC11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
APPC12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
APPC13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
APPC22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
APPC23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
APPC33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
APPC44i=${coup[12]}
APPC55i=${coup[13]}
APPC66i=${coup[14]}

#APPC11p=`python -c "print ($APPC11i+$APPC22i+$APPC33i)/3"`
#APPC12p=`python -c "print ($APPC12i+$APPC23i+$APPC13i)/3"`
#APPC44p=`python -c "print ($APPC44i+$APPC55i+$APPC66i)/3"`

echo $PRESS $APPC11i $APPC12i $APPC13i $APPC22i $APPC23i $APPC33i $APPC44i $APPC55i $APPC66i >> $dir/tests/APP/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	APPC11=`printf '%3.f' $APPC11i`
	APPC12=`printf '%3.f' $APPC12i`
	APPC13=`printf '%3.f' $APPC13i`
	APPC22=`printf '%3.f' $APPC22i`
	APPC23=`printf '%3.f' $APPC23i`
	APPC33=`printf '%3.f' $APPC33i`
	APPC44=`printf '%3.f' $APPC44i`
	APPC55=`printf '%3.f' $APPC55i`
	APPC66=`printf '%3.f' $APPC66i`
	
	C11APPd=`printf '%3.f' $C11APPd`
	C12APPd=`printf '%3.f' $C12APPd`
	C13APPd=`printf '%3.f' $C13APPd`
	C22APPd=`printf '%3.f' $C22APPd`
	C23APPd=`printf '%3.f' $C23APPd`
	C33APPd=`printf '%3.f' $C33APPd`
	C44APPd=`printf '%3.f' $C44APPd`
	C55APPd=`printf '%3.f' $C55APPd`
	C66APPd=`printf '%3.f' $C66APPd`
fi
done

# PHONONS FOR APP
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/APP/OPT.POSCAR << -
APP Ti3Nb
$APPeqp
1.000000000000	0.000000000000	0.000000000000
0.000000000000	$APPboap	0.000000000000
0.000000000000	0.000000000000	$APPcoap
Ti Nb
3 1
Direct
0.000000000000	0.000000000000	0.000000000000
0.500000000000	0.500000000000	0.000000000000
0.500000000000	0.100000000000	0.500000000000
0.000000000000	0.600000000000	0.500000000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/APP/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS
from ase.constraints import StrainFilter
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/APP/OPT.POSCAR')

atoms.set_calculator(calc)

opt = BFGS(atoms=atoms)
sf = StrainFilter(atoms)
opt.run()

# Phonon calculator
N = 10
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

G = [0, 0, 0]
X = [1./2, 0, 0]
S = [1./2, 1./2, 0]
R = [1./2, 1./2, 1./2]

point_names = ['\$R\$', '\$\Gamma\$', '\$S\$', '\$X\$']
path = [R, G, S, X]
dirs = ['\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi0]\$', '\$[0\\\bar{\\\xi}0]\$']

# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv
dft = np.loadtxt("$dftdat/Ti-Nb/APP/phonons.dat")
b = 2*np.pi/3.34

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0

omega_nd= dft[:,1]
max_thisd= np.max(omega_nd)
min_thisd= np.min(omega_nd)

for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_band = np.max([max_band, max_this, max_thisd])
    min_band = np.min([min_band, min_this, min_thisd])
    plt.plot(q, omega_n, 'g-', lw=2)

plt.scatter(b*dft[:,0], omega_nd, c='gray', s=0.1, marker='o')

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('\$\\\alpha\'\'\$ Ti\$_{3}\$Nb')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='lightgreen', edgecolor='g', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/APP_phonons.png')
ph.clean
!

python $dir/tests/APP/phonons.py > $dir/tests/APP/phonon_log
rm -f *.pckl
fi



#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< AP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t alpha' ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/AP/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM AP LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

##define simulation region and bcc grid
#lattice		custom $APlat  &
#		a1 1.0 0.0 0.0 &
#		a2 0.4994779 0.70283364 0.0 &
#		a3 0.0 0.0000000 1.2923942497 &
#		basis 0.00 0.00 0.00 basis 0.50 0.00 0.00 &
#		basis 0.00 0.50 0.00 basis 0.50 0.50 0.00 &
#		basis 0.00 0.50 0.50 basis 0.50 0.50 0.50 &
#		basis 0.00 0.00 0.50 basis 0.50 0.00 0.50 &
#		basis 0.1666667 0.8333333 0.25 basis 0.8333333 0.1666667 0.25 &
#		basis 0.8333333 0.8333333 0.75 basis 0.1666667 0.1666667 0.75 &
#		basis 0.1666667 0.1666667 0.25 basis 0.8333333 0.8333333 0.25 &
#		basis 0.8333333 0.1666667 0.75 basis 0.1666667 0.8333333 0.75
#		
#region		mybox block 0 1 0 1 0 1
#create_box	2 mybox
#
##create atoms
#create_atoms 	${idx["Nb"]} box &
#		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} & 
#		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} & 
#		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Ti"]} basis 12 ${idx["Ti"]} & 
#		basis 13 ${idx["Nb"]} basis 14 ${idx["Nb"]} basis 15 ${idx["Nb"]} basis 16 ${idx["Nb"]}
# 

box tilt large
read_data $dftdat/$elem1-$elem2/AP/AP.COORDS

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
#fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xyz fixedpoint 0 0 0
#min_style	cg
#minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx
variable	boa equal ly/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms
variable	XY equal xy
variable	XZ equal xz
variable	YZ equal yz
variable	XH equal xhi
variable	XL equal xlo
variable	YH equal yhi
variable	YL equal ylo
variable	ZH equal zhi
variable	ZL equal zlo

dump 1 all xyz 1 $dir/tests/AP/AP_RELAXED.xyz
run 0

print '\${a} \${coa} \${vat} \${spe} \${XL} \${XH} \${YL} \${YH} \${ZL} \${ZH} \${XY} \${XZ} \${YZ}'
!!

$lammps < $dir/tests/AP/eq.in > $dir/tests/AP/eq.out

APeqp=`tail -1 $dir/tests/AP/eq.out | awk '{print $1}'`
APcoap=`tail -1 $dir/tests/AP/eq.out | awk '{print $2}'`
APvop=`tail -1 $dir/tests/AP/eq.out | awk '{print $3}'`
APpote=`tail -1 $dir/tests/AP/eq.out | awk '{print $4}'`
echo "ALPHA PRIME POTENTIAL ENERGY = $APpote"
# now get relaxed triclinic box parameters
xlAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $5}'`
xhAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $6}'`
ylAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $7}'`
yhAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $8}'`
zlAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $9}'`
zhAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $10}'`
xyAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $11}'`
xzAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $12}'`
yzAP=`tail -1 $dir/tests/AP/eq.out | awk '{print $13}'`

# now write relaxed data file to be read in by e-v calculator
cat > $dir/tests/AP/AP_RELAXED.coords << -+
## comment line
16 atoms
2 atom types
$xlAP $xhAP xlo xhi
$ylAP $yhAP ylo yhi
$zlAP $zhAP zlo zhi
$xyAP $xzAP $yzAP xy xz yz
Atoms

-+

tail -16 $dir/tests/AP/AP_RELAXED.xyz | awk '{print NR, $1, $2, $3, $4}' >> $dir/tests/AP/AP_RELAXED.coords
##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for AP ------------------------------------------
#if [ "$ALPHA_PRIME_FLAG" == "TRUE" ]; then
if [ "$lammps_evsv_flag" == "TRUE" ] && [ "$lammps_AP_flag" == "TRUE" ]; then

cat > $dir/tests/AP/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR AP LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

##define simulation region and bcc grid
#lattice		custom $APeqp  &
#		a1 1.0 0.0 0.0 &
#		a2 0.4994779 0.70283364 0.0 &
#		a3 0.0 0.0000000 1.2923942497 &
#		basis 0.00 0.00 0.00 basis 0.50 0.00 0.00 &
#		basis 0.00 0.50 0.00 basis 0.50 0.50 0.00 &
#		basis 0.00 0.50 0.50 basis 0.50 0.50 0.50 &
#		basis 0.00 0.00 0.50 basis 0.50 0.00 0.50 &
#		basis 0.1666667 0.8333333 0.25 basis 0.8333333 0.1666667 0.25 &
#		basis 0.8333333 0.8333333 0.75 basis 0.1666667 0.1666667 0.75 &
#		basis 0.1666667 0.1666667 0.25 basis 0.8333333 0.8333333 0.25 &
#		basis 0.8333333 0.1666667 0.75 basis 0.1666667 0.8333333 0.75
#		
#region		mybox block 0 1 0 1 0 1
#create_box	2 mybox
#
##create atoms
#create_atoms 	${idx["Nb"]} box &
#		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} & 
#		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} & 
#		basis 9 ${idx["Ti"]} basis 10 ${idx["Ti"]} basis 11 ${idx["Ti"]} basis 12 ${idx["Ti"]} & 
#		basis 13 ${idx["Nb"]} basis 14 ${idx["Nb"]} basis 15 ${idx["Nb"]} basis 16 ${idx["Nb"]}
# 
#group		grp region mybox

box tilt large
read_data $dir/tests/AP/AP_RELAXED.coords

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}
variable tmp equal xy
variable XY0 equal \${tmp}
variable tmp equal xz
variable XZ0 equal \${tmp}
variable tmp equal yz
variable YZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/AP/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/AP/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)
variable XYI equal \${XY0}*(1-v_dmax/2)
variable XZI equal \${XZ0}*(1-v_dmax/2)
variable YZI equal \${YZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)
variable XYF equal \${XY0}*(1+v_dmax/2)
variable XZF equal \${XZ0}*(1+v_dmax/2)
variable YZF equal \${YZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} xy final \${XYI} xz final \${XZI} yz final \${YZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} xy final \${XYF} xz final \${XZF} yz final \${YZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/AP/evsv.in loop
!!
$lammps < $dir/tests/AP/evsv.in > $dir/tests/AP/evsv.out

sed -i '1d' $dir/tests/AP/evsv.dat
line=`python $BMFit $dir/tests/APP/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/AP/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	APbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	APbulkp=0
fi



else # AP lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/AP/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

APmmz=`echo "$APlat" | bc -l`
rm -f $dir/tests/AP/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$APmmz
	
	print '#N', 16, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0.499478*a, 0.702834*a, 0
	print '#Z', 0, 0, 1.29239*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.5*a, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.249739*a, 0.351418*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.749739*a, 0.351418*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.249739*a, 0.351418*a, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.749739*a, 0.351418*a, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.0,   0.0, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.5*a, 0.0, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.499651*a, 0.468555*a, 0.329099*a, 0.0, 0.0, 0.0
	print 0, 0.749913*a, 0.117139*a, 0.323099*a, 0.0, 0.0, 0.0
	print 0, 0.999651*a, 0.468555*a, 0.969296*a, 0.0, 0.0, 0.0
	print 0, 0.249913*a, 0.117139*a, 0.969296*a, 0.0, 0.0, 0.0
	print 1, 0.249913*a, 0.117139*a, 0.323099*a, 0.0, 0.0, 0.0
	print 1, 0.999651*a, 0.468555*a, 0.323099*a, 0.0, 0.0, 0.0
	print 1, 0.749913*a, 0.117139*a, 0.969296*a, 0.0, 0.0, 0.0
	print 1, 0.499651*a, 0.468555*a, 0.969296*a, 0.0, 0.0, 0.0 

" > $dir/tests/AP/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/AP/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/AP/pvsv.dat
awk -v a0=$APmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*0.908336*vol
					   print vol/16, pum, a;
				     }' data.stress >> $dir/tests/AP/pvsv.dat 

echo "#v, a" > $dir/tests/AP/evsv.dat
awk -v a0=$APmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*0.908336;
					   print vol/16, $6;
				      }' data.energy >> $dir/tests/AP/evsv.dat 


sed -i '1d' $dir/tests/AP/evsv.dat
line=`python $BMFit $dir/tests/APP/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/AP/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	APbulkp=`echo $line | awk '{print $1}'`
	#APvop=`echo $line | awk '{print $2}'`
	#APpote=`echo $line | awk '{print $3}'`
	#APeqp=`python -c "print (16.*$APvop/0.908336)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	APbulkp=0
	#APvop=0
	#APpote=0
	#APeqp=0
fi

echo -e "${prp}E-V curve...${non}"

APmmz=`echo "$APeqp" | bc -l`
rm -f $dir/tests/AP/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$APmmz
	
	print '#N', 16, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0.499478*a, 0.702834*a, 0
	print '#Z', 0, 0, 1.29239*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.5*a, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 0.249739*a, 0.351418*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.749739*a, 0.351418*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.249739*a, 0.351418*a, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.749739*a, 0.351418*a, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.0,   0.0, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.5*a, 0.0, 0.646197*a, 0.0, 0.0, 0.0
	print 0, 0.499651*a, 0.468555*a, 0.329099*a, 0.0, 0.0, 0.0
	print 0, 0.749913*a, 0.117139*a, 0.323099*a, 0.0, 0.0, 0.0
	print 0, 0.999651*a, 0.468555*a, 0.969296*a, 0.0, 0.0, 0.0
	print 0, 0.249913*a, 0.117139*a, 0.969296*a, 0.0, 0.0, 0.0
	print 1, 0.249913*a, 0.117139*a, 0.323099*a, 0.0, 0.0, 0.0
	print 1, 0.999651*a, 0.468555*a, 0.323099*a, 0.0, 0.0, 0.0
	print 1, 0.749913*a, 0.117139*a, 0.969296*a, 0.0, 0.0, 0.0
	print 1, 0.499651*a, 0.468555*a, 0.969296*a, 0.0, 0.0, 0.0 

" > $dir/tests/AP/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/AP/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/AP/pvsv.dat
awk -v a0=$APmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*0.908336*vol
					   print vol, pum, a;
				     }' data.stress >> $dir/tests/AP/pvsv.dat 

echo "#v, a" > $dir/tests/AP/evsv.dat
awk -v a0=$APmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*0.908336/16.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/AP/evsv.dat 
fi # lammps_evsv_flag AP

sed -i '1d' $dir/tests/AP/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/AP/pvsv.dat | sort -k2,2 -g | head -1`
	APPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${APPLAT[0]} $APeqp"
		APeqp=`echo $line | awk '{print $3}'`
		APvop=`echo $line | awk '{print $1}'`
		#APpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
APPLAT["0"]=$APeqp

awk -v Voo=$APvop '{print $1/Voo, $2, $3}' $dir/tests/AP/pvsv.dat > tmp; mv tmp $dir/tests/AP/pvsv.dat
fi


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< OMEGA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t omega ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/omega/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM omega LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

##define simulation region and bcc grid
variable sqrt32 equal sqrt(3)/2
lattice		custom $omglat &
		a1 1.0 0.0 0.0 &
		a2 0.5 \${sqrt32} 0.0 &
		a3 0.0 0.0 0.3 &
		basis 0.5 0.5 0.0 basis 0.0 0.5 0.0 basis 0.5 0.0 0.0 basis 0.8333333 0.3333333 0.5 &
		basis 0.3333333 0.8333333 0.5 basis 0.8333333 0.8333333 0.5 &
		basis 0.1666667 0.1666667 0.5 basis 0.6666667 0.1666667 0.5 &
		basis 0.1666667 0.6666667 0.5 basis 0.0 0.0 0.0 &
		basis 0.3333333 0.3333333 0.5 basis 0.6666667 0.6666667 0.5



# now create box and define atom basis
# a1 lattice constnats
variable 	aOMGAs equal $omglat
variable	bOMGAs equal $omglat*sqrt(3)/2
variable	cOMGAs equal 0.3*$omglat
variable 	xyOMGAs equal $omglat/2

box tilt large
region		mybox prism 0 \${aOMGAs} 0 \${bOMGAs} 0 \${cOMGAs} \${xyOMGAs} 0.0 0.0 units box
create_box	2 mybox

create_atoms ${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} &
		basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} basis 9 ${idx["Ti"]} basis 10 ${idx["Nb"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Nb"]}

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx
variable	boa equal ly/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms
variable	XY equal xy
variable	XZ equal xz
variable	YZ equal yz
variable	XH equal xhi
variable	XL equal xlo
variable	YH equal yhi
variable	YL equal ylo
variable	ZH equal zhi
variable	ZL equal zlo

dump 1 all xyz 1 $dir/tests/omega/omega_RELAXED.xyz
run 0

print '\${a} \${coa} \${vat} \${spe} \${XL} \${XH} \${YL} \${YH} \${ZL} \${ZH} \${XY} \${XZ} \${YZ}'
!!

$lammps < $dir/tests/omega/eq.in > $dir/tests/omega/eq.out

omgeqp=`tail -1 $dir/tests/omega/eq.out | awk '{print $1}'`
omgcoap=`tail -1 $dir/tests/omega/eq.out | awk '{print $2}'`
omgvop=`tail -1 $dir/tests/omega/eq.out | awk '{print $3}'`
omgpote=`tail -1 $dir/tests/omega/eq.out | awk '{print $4}'`

# now get relaxed triclinic box parameters
xlomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $5}'`
xhomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $6}'`
ylomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $7}'`
yhomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $8}'`
zlomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $9}'`
zhomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $10}'`
xyomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $11}'`
xzomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $12}'`
yzomega=`tail -1 $dir/tests/omega/eq.out | awk '{print $13}'`

# now write relaxed data file to be read in by e-v calculator
cat > $dir/tests/omega/omega_RELAXED.coords << -+
## comment line
12 atoms
2 atom types
$xlomega $xhomega xlo xhi
$ylomega $yhomega ylo yhi
$zlomega $zhomega zlo zhi
$xyomega $xzomega $yzomega xy xz yz
Atoms

-+

tail -12 $dir/tests/omega/omega_RELAXED.xyz | awk '{print NR, $1, $2, $3, $4}' >> $dir/tests/omega/omega_RELAXED.coords

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for omega ------------------------------------------

cat > $dir/tests/omega/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR omega LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

##define simulation region and bcc grid
variable sqrt32 equal sqrt(3)/2
lattice		custom $omgeqp &
		a1 1.0 0.0 0.0 &
		a2 0.5 \${sqrt32} 0.0 &
		a3 0.0 0.0 $omgcoap &
		basis 0.5 0.5 0.0 basis 0.0 0.5 0.0 basis 0.5 0.0 0.0 basis 0.8333333 0.3333333 0.5 &
		basis 0.3333333 0.8333333 0.5 basis 0.8333333 0.8333333 0.5 &
		basis 0.1666667 0.1666667 0.5 basis 0.6666667 0.1666667 0.5 &
		basis 0.1666667 0.6666667 0.5 basis 0.0 0.0 0.0 &
		basis 0.3333333 0.3333333 0.5 basis 0.6666667 0.6666667 0.5



# now create box and define atom basis
# a1 lattice constnats
variable 	aOMGAs equal $omgeqp
variable	bOMGAs equal $omgeqp*sqrt(3)/2
variable	cOMGAs equal $omgcoap*$omgeqp
variable 	xyOMGAs equal $omgeqp/2

box tilt large
region		mybox prism 0 \${aOMGAs} 0 \${bOMGAs} 0 \${cOMGAs} \${xyOMGAs} 0.0 0.0 units box
create_box	2 mybox

create_atoms ${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} &
		basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} basis 9 ${idx["Ti"]} basis 10 ${idx["Nb"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Nb"]}

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}
variable tmp equal xy
variable XY0 equal \${tmp}
variable tmp equal xz
variable XZ0 equal \${tmp}
variable tmp equal yz
variable YZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/omega/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/omega/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)
variable XYI equal \${XY0}*(1-v_dmax/2)
variable XZI equal \${XZ0}*(1-v_dmax/2)
variable YZI equal \${YZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)
variable XYF equal \${XY0}*(1+v_dmax/2)
variable XZF equal \${XZ0}*(1+v_dmax/2)
variable YZF equal \${YZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} xy final \${XYI} xz final \${XZI} yz final \${YZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} xy final \${XYF} xz final \${XZF} yz final \${YZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/omega/evsv.in loop
!!
$lammps < $dir/tests/omega/evsv.in > $dir/tests/omega/evsv.out

sed -i '1d' $dir/tests/omega/evsv.dat
line=`python $BMFit $dir/tests/omega/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/omega/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	omgbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	omgbulkp=0
fi

else # omega lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/omega/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

omgmmz=`echo "$omglat" | bc -l`
echo $omgmmz
rm -f $dir/tests/omega/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$omgmmz
	
	print '#N', 12, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0.5*a, 0.8660254*a, 0
	print '#Z', 0, 0, 0.3*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, 0.75*a, 0.433013*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.25*a, 0.433013*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.50*a, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 1.00*a, 0.288674*a, 0.15*a, 0.0, 0.0, 0.0
	print 0, 0.75*a, 0.721688*a, 0.15*a, 0.0, 0.0, 0.0	
	print 0, 1.25*a, 0.721688*a, 0.15*a, 0.0, 0.0, 0.0	
	print 0, 0.25*a, 0.144337*a, 0.15*a, 0.0, 0.0, 0.0	
	print 0, 0.75*a, 0.144337*a, 0.15*a, 0.0, 0.0, 0.0
	print 0, 0.50*a, 0.577351*a, 0.15*a, 0.0, 0.0, 0.0
	print 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 1, 0.50*a, 0.288672*a, 0.15*a, 0.0, 0.0, 0.0
	print 1, 1.00*a, 0.577351*a, 0.15*a, 0.0, 0.0, 0.0	

" > $dir/tests/omega/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/omega/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/omega/pvsv.dat
awk -v a0=$omgmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*0.259808*vol;
					   print vol/12, pum, a;
				     }' data.stress >> $dir/tests/omega/pvsv.dat 

echo "#v, a" > $dir/tests/omega/evsv.dat
awk -v a0=$omgmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*0.259808;
					   print vol/12, $6;
				      }' data.energy >> $dir/tests/omega/evsv.dat 


sed -i '1d' $dir/tests/omega/evsv.dat
line=`python $BMFit $dir/tests/omega/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/omega/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	omgbulkp=`echo $line | awk '{print $1}'`
	omgvop=`echo $line | awk '{print $2}'`
	omgpote=`echo $line | awk '{print $3}'`
	omgeqp=`python -c "print (12.*$omgvop/0.259808)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	omgbulkp=0
	omgvop=0
	omgpote=0
	omgeqp=0
fi

echo -e "${prp}E-V curve...${non}"

omgmmz=`echo "$omgeqp" | bc -l`
rm -f $dir/tests/omega/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $omgmmz*math.pow(a,1./3)
	
	print '#N', 12, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0.5*a, 0.8660254*a, 0
	print '#Z', 0, 0, 0.3*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0, 0.75*a, 0.433013*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.25*a, 0.433013*a, 0.0, 0.0, 0.0, 0.0
	print 0, 0.50*a, 0.0, 0.0, 0.0, 0.0, 0.0
	print 0, 1.00*a, 0.288674*a, 0.15*a, 0.0, 0.0, 0.0
	print 0, 0.75*a, 0.721688*a, 0.15*a, 0.0, 0.0, 0.0	
	print 0, 1.25*a, 0.721688*a, 0.15*a, 0.0, 0.0, 0.0	
	print 0, 0.25*a, 0.144337*a, 0.15*a, 0.0, 0.0, 0.0	
	print 0, 0.75*a, 0.144337*a, 0.15*a, 0.0, 0.0, 0.0
	print 0, 0.50*a, 0.577351*a, 0.15*a, 0.0, 0.0, 0.0
	print 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print 1, 0.50*a, 0.288672*a, 0.15*a, 0.0, 0.0, 0.0
	print 1, 1.00*a, 0.577351*a, 0.15*a, 0.0, 0.0, 0.0	

" > $dir/tests/omega/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/omega/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/omega/pvsv.dat
awk -v a0=$omgmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*0.259808*vol
					   print vol/12, pum, a;
				     }' data.stress >> $dir/tests/omega/pvsv.dat 

echo "#v, a" > $dir/tests/omega/evsv.dat
awk -v a0=$omgmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*0.259808/12.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/omega/evsv.dat 
fi # lammps_evsv_flag omega

sed -i '1d' $dir/tests/omega/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/omega/pvsv.dat | sort -k2,2 -g | head -1`
	omgPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${omgPLAT[0]} $omgeqp"
		#omgeqp=`echo $line | awk '{print $3}'`
		#omgvop=`echo $line | awk '{print $1}'`
		#omegapote=`echo $line | awk '{print 1000*$4}'`
	fi
done
omgPLAT["0"]=$omgeqp

awk -v Voo=$omgvop '{print $1/Voo, $2, $3}' $dir/tests/omega/pvsv.dat > tmp; mv tmp $dir/tests/omega/pvsv.dat

# ---------------------- elastic constants for omega ---------------------------------
echo -e "${prp}Elastic constants:${non}"
echo "# pressure, c11, c12, c44" > $dir/tests/omega/C_VS_P.dat
for PRESS in $PRESSEQ; do 
echo "$PRESS GPa..."
for jj in $(seq 1 7); do

printf "\t $jj "

if [ "$lammps_elcon_flag" == "TRUE" ]; then

P=`echo $PRESS*10000 | bc -l`
case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac
#case $jj in
#	1) XD="1+v_d"; YD="sqrt(3)*(1+v_d)/2"; ZD="1+v_d"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	2) XD="1+v_d"; YD="sqrt(3)*(1-v_d)/2"; ZD="1/(1-v_d2)"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	3) XD="1/(1-v_d2)"; YD="sqrt(3)*(1+v_d)/2"; ZD="1-v_d"; XYD="-0.5/(1-v_d2)"; XZD="0"; YZD="0" ;;
#	4) XD="1-v_d"; YD="sqrt(3)/(2-2*v_d2)"; ZD="1+v_d"; XYD="-0.5*(1-v_d)"; XZD="0"; YZD="0" ;;
#	5) XD="1/(1-(v_d2)/4)"; YD="sqrt(3)/2"; ZD="1"; XYD="-0.5/(1-(v_d2)/4)"; XZD="0"; YZD="v_d" ;;
#	6) XD="1"; YD="sqrt(3)/((1-(v_d2)/4)*2)"; ZD="1"; XYD="-0.5"; XZD="v_d"; YZD="0" ;;
#	7) XD="1"; YD="sqrt(3)/2"; ZD="1/(1-(v_d2)/4)"; XYD="-0.5+v_d"; XZD="0"; YZD="0" ;;
#esac

cat > $dir/tests/omega/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

variable boxa equal ${omgPLAT["$PRESS"]}
variable boxb equal ${omgPLAT["$PRESS"]}*sqrt(3)/2
variable boxc equal ${omgPLAT["$PRESS"]}*$omgcoap
variable boxxy equal ${omgPLAT["$PRESS"]}*(0.5)
variable boxxz equal 0
variable boxyz equal 0


##define simulation region and omega grid
lattice		custom ${omgPLAT["$PRESS"]}  &
		a1 1.0 0.0 0.0 &
		a2 0.50 0.86602540 0.0 &
		a3 0.0 0.0 $omgcoap &
		basis 0.50 0.50 0.00 basis 0.00 0.50 0.00 &
		basis 0.50 0.00 0.00 basis 0.833333 0.333333 0.50 &
		basis 0.333333 0.833333 0.50 basis 0.833333 0.833333 0.50 &
		basis 0.1666667 0.1666667 0.50 basis 0.6666667 0.1666667 0.50 &
		basis 0.1666667 0.6666667 0.55 basis 0.0000000 0.0000000 0.00 &
		basis 0.3333333 0.3333333 0.50 basis 0.6666667 0.6666667 0.50
		
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} & 
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} & 
		basis 9 ${idx["Ti"]} basis 10 ${idx["Nb"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Nb"]} 
 

dump D all xyz 1 elcon$jj.xyz
#read_data $dir/tests/omega/omega_RELAXED.coords
 
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
timestep 0.001
#fix rel all press/berendsen iso $P $P 1.0
#run 1000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

variable tmp equal lx
variable LX equal \${tmp}
variable tmp equal ly
variable LY equal \${tmp}
variable tmp equal lz
variable LZ equal \${tmp}
variable XY equal xy
variable lat equal \${tmp}
variable YZ equal yz
variable lat equal \${tmp}
variable XZ equal xz
variable lat equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/omega/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)

variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/omega/elcon.lin loop 
!!

$lammps < $dir/tests/omega/elcon.lin > $dir/tests/omega/elcon.out
sed -i '1d' $dir/tests/omega/stresses$jj.dat

else # else stress strain-curves with meamz

cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/omega/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
strain=0.002
rm -f $dir/tests/omega/elcon$jj.conf
for i in $(seq -$DELPT $DELPT); do
	omgmmz=`echo ${omgPLAT["$PRESS"]} | bc -l`
	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del omega $omgmmz conf >> $dir/tests/omega/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/omega/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=-pc*$5; getline	
					syy=-pc*$5; getline	
					szz=-pc*$5; getline	
					sxy=-pc*$5; getline	
					syz=-pc*$5; getline	
					szx=-pc*$5; 	
					print del, sxx, syy, szz, syz, szx, sxy
				     }' data.stress >> $dir/tests/omega/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""

rm -f $dir/tests/omega/EC_fits.dat; touch $dir/tests/omega/EC_fits.dat

# first row fits
awk '{print $1, $2}' $dir/tests/omega/stresses1.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/omega/stresses1.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/omega/stresses1.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/omega/stresses2.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/omega/stresses2.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/omega/stresses2.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/omega/stresses3.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/omega/stresses3.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/omega/stresses3.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/omega/stresses4.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/omega/stresses4.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/omega/stresses4.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/omega/stresses5.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/omega/stresses6.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/omega/stresses7.dat > tmp; python $fit tmp >> $dir/tests/omega/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/omega/EC_fits.dat` )
omgC11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
omgC12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
omgC13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
omgC22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
omgC23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
omgC33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
omgC44i=${coup[12]}
omgC55i=${coup[13]}
omgC66i=${coup[14]}

omgC11p=`python -c "print ($omgC11i+$omgC22i)/2"`
omgC13p=`python -c "print ($omgC13i+$omgC23i)/2"`
omgC44p=`python -c "print ($omgC44i+$omgC55i)/2"`

echo $PRESS $omgC11p $omgC12i $omgC13p $omgC33i $omgC44p >> $dir/tests/omega/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	omgC11=`printf '%3.f' $omgC11p`
	omgC12=`printf '%3.f' $omgC12i`
	omgC13=`printf '%3.f' $omgC13p`
	omgC33=`printf '%3.f' $omgC33i`
	omgC44=`printf '%3.f' $omgC44p`
	omgC66=`printf '%3.f' $omgC66i`

	C11omgd=`printf '%3.f' $C11omgd`
	C12omgd=`printf '%3.f' $C12omgd`
	C13omgd=`printf '%3.f' $C13omgd`
	C33omgd=`printf '%3.f' $C33omgd`
	C44omgd=`printf '%3.f' $C44omgd`
	C66omgd=`printf '%3.f' $C66omgd`
fi
done


# PHONONS FOR omega
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/omega/OPT.POSCAR << -
omega-Ti3Nb 
$omgeqp
1.000000000000	0.000000000000	0.000000000000
0.500000000000	0.866025403780	0.000000000000
0.000000000000	0.000000000000	$omgcoap
Ti Nb
9 3
Direct
0.500000000000	0.500000000000	0.000000000000
0.000000000000	0.500000000000	0.000000000000
0.500000000000	0.000000000000	0.000000000000
0.833333333333	0.333333333333	0.500000000000
0.333333333333	0.833333333333	0.500000000000
0.833333333333	0.833333333333	0.500000000000
0.166666666667	0.166666666667	0.500000000000
0.666666666667	0.166666666667	0.500000000000
0.166666666667	0.666666666667	0.500000000000
0.000000000000	0.000000000000	0.000000000000
0.333333333333	0.333333333333	0.500000000000
0.666666666667	0.666666666667	0.500000000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/omega/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS
from ase.units import _hbar, _e
from ase.io import read
from ase.constraints import StrainFilter
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/omega/OPT.POSCAR')

atoms.set_calculator(calc)

opt = BFGS(atoms=atoms)
sf = StrainFilter(atoms)
opt.run()

# Phonon calculator
N = 4
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

points = ibz_points['hexagonal']
G = points['Gamma']
H = points['H']
K = points['K']
M = points['M']
A = points['A']

point_names = ['\$\Gamma\$', '\$K\$', '\$M\$', '\$\Gamma\$', '\$A\$']
dirs = ['\$[\\\xi\\\xi0]\$', '' ,'\$[\\\xi00]\$', '\$[00\\\xi]\$']
path = [G, K, M, G, A]


# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_band = np.max([max_band, max_this])
    min_band = np.min([min_band, min_this])
    plt.plot(q, omega_n, 'g-', lw=2)

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('\$\\omega\$ Ti\$_{3}\$Nb')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='lightgreen', edgecolor='g', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/omega_phonons.png')
ph.clean
!

python $dir/tests/omega/phonons.py > $dir/tests/omega/phonon_log
rm -f *.pckl
fi

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< D019 Ti3Nb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t D019 Ti3Nb ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/D019Ti3Nb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM D019Ti3Nb LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

variable LAT equal $D019Ti3Nblat
variable COA equal $D019Ti3Nbcoa
variable XY  equal 0.5*$D019Ti3Nblat
variable Y   equal sqrt(3)/2
variable LY  equal $D019Ti3Nblat*\${Y}
variable LZ  equal $D019Ti3Nblat*\${COA}

lattice custom \${LAT} &
	a1 1.00		0.00		0.00 &
	a2 0.50		\${Y}		0.00 &
	a3 0.00		0.00		\${COA} &
	basis 0.0 0.0 0.0 basis 0.66667 0.66667 0.5 &
	basis 0.16667 0.16667 0.50 basis 0.50 0.0 0.0 &
	basis 0.0 0.50 0.0 basis 0.16667 0.66667 0.50 &
	basis 0.66667 0.16667 0.50 basis 0.50 0.50 0.0

box tilt large
region mybox prism 0 \${LAT} 0 \${LY} 0 \${LZ} \${XY} 0 0 units box
create_box 2 mybox
create_atoms ${idx["Nb"]} box basis 1 ${idx["Nb"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} &
		   basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax aniso 0.0 couple xy scalexy yes
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	a equal lx
variable	coa equal lz/lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

run 0

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/D019Ti3Nb/eq.in > $dir/tests/D019Ti3Nb/eq.out

D019Ti3Nbeqp=`tail -1 $dir/tests/D019Ti3Nb/eq.out | awk '{print $1}'`
D019Ti3Nbcoap=`tail -1 $dir/tests/D019Ti3Nb/eq.out | awk '{print $2}'`
D019Ti3Nbvop=`tail -1 $dir/tests/D019Ti3Nb/eq.out | awk '{print $3}'`
D019Ti3Nbpote=`tail -1 $dir/tests/D019Ti3Nb/eq.out | awk '{print $4}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for D019Ti3Nb ------------------------------------------

cat > $dir/tests/D019Ti3Nb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR D019Ti3Nb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable LAT equal $D019Ti3Nbeqp
variable COA equal $D019Ti3Nbcoap
variable XY  equal 0.5*$D019Ti3Nbeqp
variable Y   equal sqrt(3)/2
variable LY  equal $D019Ti3Nbeqp*\${Y}
variable LZ  equal $D019Ti3Nbeqp*\${COA}

lattice custom \${LAT} &
	a1 1.00		0.00		0.00 &
	a2 0.50		\${Y}		0.00 &
	a3 0.00		0.00		\${COA} &
	basis 0.0 0.0 0.0 basis 0.66667 0.66667 0.5 &
	basis 0.16667 0.16667 0.50 basis 0.50 0.0 0.0 &
	basis 0.0 0.50 0.0 basis 0.16667 0.66667 0.50 &
	basis 0.66667 0.16667 0.50 basis 0.50 0.50 0.0

box tilt large
region mybox prism 0 \${LAT} 0 \${LY} 0 \${LZ} \${XY} 0 0 units box
create_box 2 mybox
create_atoms ${idx["Nb"]} box basis 1 ${idx["Nb"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} &
		   basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]}

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz xy pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}
variable tmp equal xy
variable XY0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/D019Ti3Nb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/D019Ti3Nb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)
variable XYI equal 0.5*\${LXI}

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)
variable XYF equal 0.5*\${LXF}

variable xys equal xy/lx

thermo_style custom step pe etotal press vol lx ly lz xy v_xys pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} xy final \${XYI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} xy final \${XYF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/D019Ti3Nb/evsv.in loop
!!
$lammps < $dir/tests/D019Ti3Nb/evsv.in > $dir/tests/D019Ti3Nb/evsv.out

sed -i '1d' $dir/tests/D019Ti3Nb/evsv.dat
line=`python $BMFit $dir/tests/D019Ti3Nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/D019Ti3Nb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	D019Ti3Nbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	D019Ti3Nbbulkp=0
fi

else # D019Ti3Nb lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/D019Ti3Nb/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

D019Ti3Nbmmz=`echo "$D019Ti3Nblat" | bc -l`
echo $D019Ti3Nbmmz
rm -f $dir/tests/D019Ti3Nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$D019Ti3Nbmmz
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0.50, 0.866025, 0
	print '#Z', 0, 0, $D019coa*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.66667,0.66667,0.5,	0,	0,	0
	print 0,	0.16667,0.16667,0.5,	0,	0,	0
	print 0,	0.5,	0,	0.0,	0,	0,	0
	print 0,	0,	0.5,	0.0,	0,	0,	0
	print 0,	0.16667,0.66667,0.5,	0,	0,	0
	print 0,	0.66667,0.16667,0.5,	0,	0,	0
	print 0,	0.5,	0.5,	0,	0,	0,	0
" > $dir/tests/D019Ti3Nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/D019Ti3Nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/D019Ti3Nb/pvsv.dat
awk -v a0=$D019Ti3Nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol*0.866025;
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/D019Ti3Nb/pvsv.dat 

echo "#v, a" > $dir/tests/D019Ti3Nb/evsv.dat
awk -v a0=$D019Ti3Nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*0.866025;
					   print vol/8, $6;
				      }' data.energy >> $dir/tests/D019Ti3Nb/evsv.dat 


sed -i '1d' $dir/tests/D019Ti3Nb/evsv.dat
line=`python $BMFit $dir/tests/D019Ti3Nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/D019Ti3Nb/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	D019Ti3Nbbulkp=`echo $line | awk '{print $1}'`
	D019Ti3Nbvop=`echo $line | awk '{print $2}'`
	D019Ti3Nbpote=`echo $line | awk '{print $3}'`
	D019Ti3Nbeqp=`python -c "print (8.*$D019Ti3Nbvop/0.866025)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	D019Ti3Nbbulkp=0
	D019Ti3Nbvop=0
	D019Ti3Nbpote=0
	D019Ti3Nbeqp=0
fi

echo -e "${prp}E-V curve...${non}"

D019Ti3Nbmmz=`echo "$D019Ti3Nbeqp" | bc -l`
rm -f $dir/tests/D019Ti3Nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $D019Ti3Nbmmz*math.pow(a,1./3)
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0.50, 0.866025, 0
	print '#Z', 0, 0, $D019coap*a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.66667,0.66667,0.5,	0,	0,	0
	print 0,	0.16667,0.16667,0.5,	0,	0,	0
	print 0,	0.5,	0,	0.0,	0,	0,	0
	print 0,	0,	0.5,	0.0,	0,	0,	0
	print 0,	0.16667,0.66667,0.5,	0,	0,	0
	print 0,	0.66667,0.16667,0.5,	0,	0,	0
	print 0,	0.5,	0.5,	0,	0,	0,	0

" > $dir/tests/D019Ti3Nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/D019Ti3Nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/D019Ti3Nb/pvsv.dat
awk -v a0=$D019Ti3Nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*0.866025*vol
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/D019Ti3Nb/pvsv.dat 

echo "#v, a" > $dir/tests/D019Ti3Nb/evsv.dat
awk -v a0=$D019Ti3Nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0*0.866025/8.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/D019Ti3Nb/evsv.dat 
fi # lammps_evsv_flag D019Ti3Nb

sed -i '1d' $dir/tests/D019Ti3Nb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/D019Ti3Nb/pvsv.dat | sort -k2,2 -g | head -1`
	D019Ti3NbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${D019Ti3NbPLAT[0]} $D019Ti3Nbeqp"
		#D019Ti3Nbeqp=`echo $line | awk '{print $3}'`
		#D019Ti3Nbvop=`echo $line | awk '{print $1}'`
		#D019Ti3Nbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
D019Ti3NbPLAT["0"]=$D019Ti3Nbeqp

awk -v Voo=$D019Ti3Nbvop '{print $1/Voo, $2, $3}' $dir/tests/D019Ti3Nb/pvsv.dat > tmp; mv tmp $dir/tests/D019Ti3Nb/pvsv.dat


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< A15 Ti3Nb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t A15 Ti3Nb ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/A15Ti3Nb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM A15Ti3Nb LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

variable SIZE equal $A15Ti3Nblat

lattice		custom $A15Ti3Nblat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Nb"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} & 
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]}
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax iso 0.0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

run 0

print '\${a} \${vat} \${spe}'
!!

$lammps < $dir/tests/A15Ti3Nb/eq.in > $dir/tests/A15Ti3Nb/eq.out

A15Ti3Nbeqp=`tail -1 $dir/tests/A15Ti3Nb/eq.out | awk '{print $1}'`
A15Ti3Nbvop=`tail -1 $dir/tests/A15Ti3Nb/eq.out | awk '{print $2}'`
A15Ti3Nbpote=`tail -1 $dir/tests/A15Ti3Nb/eq.out | awk '{print $3}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for A15Ti3Nb ------------------------------------------

cat > $dir/tests/A15Ti3Nb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR A15Ti3Nb LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable SIZE equal $A15Ti3Nbeqp

lattice		custom $A15Ti3Nbeqp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Nb"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} & 
		basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]}

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/A15Ti3Nb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/A15Ti3Nb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/A15Ti3Nb/evsv.in loop
!!
$lammps < $dir/tests/A15Ti3Nb/evsv.in > $dir/tests/A15Ti3Nb/evsv.out

sed -i '1d' $dir/tests/A15Ti3Nb/evsv.dat
line=`python $BMFit $dir/tests/A15Ti3Nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15Ti3Nb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15Ti3Nbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15Ti3Nbbulkp=0
fi

else # A15Ti3Nb lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/A15Ti3Nb/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

A15Ti3Nbmmz=`echo "$A15Ti3Nblat" | bc -l`
echo $A15Ti3Nbmmz
rm -f $dir/tests/A15Ti3Nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$A15Ti3Nbmmz
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.5,	0.5,	0.5,	0,	0,	0
	print 0,	0.5,	0,	0.25,	0,	0,	0
	print 0,	0.5,	0,	0.75,	0,	0,	0
	print 0,	0,	0.25,	0.5,	0,	0,	0
	print 0,	0,	0.75,	0.5,	0,	0,	0
	print 0,	0.25,	0.5,	0,	0,	0,	0
	print 0,	0.75,	0.5,	0,	0,	0,	0
" > $dir/tests/A15Ti3Nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15Ti3Nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15Ti3Nb/pvsv.dat
awk -v a0=$A15Ti3Nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol;
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15Ti3Nb/pvsv.dat 

echo "#v, a" > $dir/tests/A15Ti3Nb/evsv.dat
awk -v a0=$A15Ti3Nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/8, $6;
				      }' data.energy >> $dir/tests/A15Ti3Nb/evsv.dat 


sed -i '1d' $dir/tests/A15Ti3Nb/evsv.dat
line=`python $BMFit $dir/tests/A15Ti3Nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15Ti3Nb/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15Ti3Nbbulkp=`echo $line | awk '{print $1}'`
	A15Ti3Nbvop=`echo $line | awk '{print $2}'`
	A15Ti3Nbpote=`echo $line | awk '{print $3}'`
	A15Ti3Nbeqp=`python -c "print (8.*$A15Ti3Nbvop)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15Ti3Nbbulkp=0
	A15Ti3Nbvop=0
	A15Ti3Nbpote=0
	A15Ti3Nbeqp=0
fi

echo -e "${prp}E-V curve...${non}"

A15Ti3Nbmmz=`echo "$A15Ti3Nbeqp" | bc -l`
rm -f $dir/tests/A15Ti3Nb/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $A15Ti3Nbmmz*math.pow(a,1./3)
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 1,	0,	0,	0,	0,	0,	0
	print 1,	0.5,	0.5,	0.5,	0,	0,	0
	print 0,	0.5,	0,	0.25,	0,	0,	0
	print 0,	0.5,	0,	0.75,	0,	0,	0
	print 0,	0,	0.25,	0.5,	0,	0,	0
	print 0,	0,	0.75,	0.5,	0,	0,	0
	print 0,	0.25,	0.5,	0,	0,	0,	0
	print 0,	0.75,	0.5,	0,	0,	0,	0

" > $dir/tests/A15Ti3Nb/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15Ti3Nb/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15Ti3Nb/pvsv.dat
awk -v a0=$A15Ti3Nbmmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15Ti3Nb/pvsv.dat 

echo "#v, a" > $dir/tests/A15Ti3Nb/evsv.dat
awk -v a0=$A15Ti3Nbmmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/8.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/A15Ti3Nb/evsv.dat 
fi # lammps_evsv_flag A15Ti3Nb

sed -i '1d' $dir/tests/A15Ti3Nb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/A15Ti3Nb/pvsv.dat | sort -k2,2 -g | head -1`
	A15Ti3NbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${A15Ti3NbPLAT[0]} $A15Ti3Nbeqp"
		#A15Ti3Nbeqp=`echo $line | awk '{print $3}'`
		#A15Ti3Nbvop=`echo $line | awk '{print $1}'`
		#A15Ti3Nbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
A15Ti3NbPLAT["0"]=$A15Ti3Nbeqp

awk -v Voo=$A15Ti3Nbvop '{print $1/Voo, $2, $3}' $dir/tests/A15Ti3Nb/pvsv.dat > tmp; mv tmp $dir/tests/A15Ti3Nb/pvsv.dat

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< L10 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t FCT L10 TiNb... ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/L10/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

# define box variables
variable boxa equal $L10lat
variable boxb equal $L10lat
variable boxc equal $L10lat*$L10coa
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

#define simulation region and bcc grid
lattice		custom $L10lat &
		a1 1.00 0.00 0.00 a2 0.00 1.00 0.00 a3 0.00 0.00 $L10coa &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 basis 0.00 0.50 0.50

region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} 
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/L10/eq.in > $dir/tests/L10/eq.out

L10eqp=`tail -1 $dir/tests/L10/eq.out | awk '{print $1}'`
L10coap=`tail -1 $dir/tests/L10/eq.out | awk '{print $2}'`
L10vop=`tail -1 $dir/tests/L10/eq.out | awk '{print $3}'`
L10pote=`tail -1 $dir/tests/L10/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for L10 -----------------------------------------

cat > $dir/tests/L10/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR L10 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable boxa equal $L10eqp
variable boxb equal $L10eqp
variable boxc equal $L10eqp*$L10coap
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $L10lat &
		a1 1.00 0.00 0.00 a2 0.00 1.00 0.00 a3 0.00 0.00 $L10coap &
		basis 0.00 0.00 0.00 basis 0.50 0.50 0.00 basis 0.50 0.00 0.50 basis 0.00 0.50 0.50
		
region		mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} \${boxxz} \${boxyz} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/L10/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/L10/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd equal \${boxa}*(1-v_dmax/2)
variable xf equal \${boxa}*(1+v_dmax/2)
variable yd equal \${boxb}*(1-v_dmax/2)
variable yf equal \${boxb}*(1+v_dmax/2)
variable zd equal \${boxc}*(1-v_dmax/2)
variable zf equal \${boxc}*(1+v_dmax/2)
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/L10/evsv.in loop
!!
$lammps < $dir/tests/L10/evsv.in > $dir/tests/L10/evsv.out

sed -i '1d' $dir/tests/L10/evsv.dat
line=`python $BMFit $dir/tests/L10/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/L10/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	L10bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	L10bulkp=0
fi

sed -i '1d' $dir/tests/L10/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/L10/pvsv.dat | sort -k2,2 -g | head -1`
	L10PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${L10PLAT[0]} $L10eqp"
		#L10eqp=`echo $line | awk '{print $3}'`
		#L10vop=`echo $line | awk '{print $1}'`
		#L10pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
L10PLAT["0"]=$L10eqp
awk -v Voo=$L10vop '{print $1/Voo, $2, $3}' $dir/tests/L10/pvsv.dat > tmp; mv tmp $dir/tests/L10/pvsv.dat
fi


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< A15 TiNb3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t A15 TiNb3 ${non}"
echo -e "${prp}equilibria...${non}"
if [ "$lammps_evsv_flag" == "TRUE" ]; then

cat > $dir/tests/A15TiNb3/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM A15TiNb3 LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid

variable SIZE equal $A15TiNb3lat

lattice		custom $A15TiNb3lat &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} & 
		basis 5 ${idx["Nb"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Nb"]} basis 8 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz cellalpha cellbeta cellgamma
thermo 1

# minimize
fix 1 all box/relax iso 0.0
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

run 0

print '\${a} \${vat} \${spe}'
!!

$lammps < $dir/tests/A15TiNb3/eq.in > $dir/tests/A15TiNb3/eq.out

A15TiNb3eqp=`tail -1 $dir/tests/A15TiNb3/eq.out | awk '{print $1}'`
A15TiNb3vop=`tail -1 $dir/tests/A15TiNb3/eq.out | awk '{print $2}'`
A15TiNb3pote=`tail -1 $dir/tests/A15TiNb3/eq.out | awk '{print $3}'`

##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for A15TiNb3 ------------------------------------------

cat > $dir/tests/A15TiNb3/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR A15TiNb3 LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable SIZE equal $A15TiNb3eqp

lattice		custom $A15TiNb3eqp &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		
region		mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} & 
		basis 5 ${idx["Nb"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Nb"]} basis 8 ${idx["Nb"]}

mass		1 $mass1 
mass		2 $mass2


#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz cellalpha cellbeta cellgamma
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000
variable tmp equal lx
variable LX0 equal \${tmp}
variable tmp equal ly
variable LY0 equal \${tmp}
variable tmp equal lz
variable LZ0 equal \${tmp}

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/A15TiNb3/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/A15TiNb3/pvsv.dat screen no title "# V/atom | pressure | lattice constant"
reset_timestep 0

variable LXI equal \${LX0}*(1-v_dmax/2)
variable LYI equal \${LY0}*(1-v_dmax/2)
variable LZI equal \${LZ0}*(1-v_dmax/2)

variable LXF equal \${LX0}*(1+v_dmax/2)
variable LYF equal \${LY0}*(1+v_dmax/2)
variable LZF equal \${LZ0}*(1+v_dmax/2)

change_box all x final 0 \${LXI} y final 0 \${LYI} z final 0 \${LZI} remap units box

fix def all deform 1 x final 0 \${LXF} y final 0 \${LYF} z final 0 \${LZF} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/A15TiNb3/evsv.in loop
!!
$lammps < $dir/tests/A15TiNb3/evsv.in > $dir/tests/A15TiNb3/evsv.out

sed -i '1d' $dir/tests/A15TiNb3/evsv.dat
line=`python $BMFit $dir/tests/A15TiNb3/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15TiNb3/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15TiNb3bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15TiNb3bulkp=0
fi

else # A15TiNb3 lammps flag
## calculate ev curve with fitting code!
numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/A15TiNb3/evsv.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

A15TiNb3mmz=`echo "$A15TiNb3lat" | bc -l`
echo $A15TiNb3mmz
rm -f $dir/tests/A15TiNb3/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = math.pow(a,1./3)*$A15TiNb3mmz
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0,	0,	0,	0,	0,	0,	0
	print 0,	0.5,	0.5,	0.5,	0,	0,	0
	print 1,	0.5,	0,	0.25,	0,	0,	0
	print 1,	0.5,	0,	0.75,	0,	0,	0
	print 1,	0,	0.25,	0.5,	0,	0,	0
	print 1,	0,	0.75,	0.5,	0,	0,	0
	print 1,	0.25,	0.5,	0,	0,	0,	0
	print 1,	0.75,	0.5,	0,	0,	0,	0
" > $dir/tests/A15TiNb3/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15TiNb3/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15TiNb3/pvsv.dat
awk -v a0=$A15TiNb3mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol;
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15TiNb3/pvsv.dat 

echo "#v, a" > $dir/tests/A15TiNb3/evsv.dat
awk -v a0=$A15TiNb3mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0;
					   print vol/8, $6;
				      }' data.energy >> $dir/tests/A15TiNb3/evsv.dat 


sed -i '1d' $dir/tests/A15TiNb3/evsv.dat
line=`python $BMFit $dir/tests/A15TiNb3/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/A15TiNb3/evsv.dat\"];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
#echo $line
var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	A15TiNb3bulkp=`echo $line | awk '{print $1}'`
	A15TiNb3vop=`echo $line | awk '{print $2}'`
	A15TiNb3pote=`echo $line | awk '{print $3}'`
	A15TiNb3eqp=`python -c "print (8.*$A15TiNb3vop)**(1./3)"`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	A15TiNb3bulkp=0
	A15TiNb3vop=0
	A15TiNb3pote=0
	A15TiNb3eqp=0
fi

echo -e "${prp}E-V curve...${non}"

A15TiNb3mmz=`echo "$A15TiNb3eqp" | bc -l`
rm -f $dir/tests/A15TiNb3/evsv.conf; 
python -c "import math
for i in range(0, $numpts+1):
	a = (0.80 + (float(i)/$numpts)*0.4)
	a = $A15TiNb3mmz*math.pow(a,1./3)
	
	print '#N', 8, 2
	print '##'
	print '#X', a, 0, 0
	print '#Y', 0, a, 0
	print '#Z', 0, 0, a
	print '#E', 0.0
	print '#S', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	print '#F'
	print 0,	0,	0,	0,	0,	0,	0
	print 0,	0.5,	0.5,	0.5,	0,	0,	0
	print 1,	0.5,	0,	0.25,	0,	0,	0
	print 1,	0.5,	0,	0.75,	0,	0,	0
	print 1,	0,	0.25,	0.5,	0,	0,	0
	print 1,	0,	0.75,	0.5,	0,	0,	0
	print 1,	0.25,	0.5,	0,	0,	0,	0
	print 1,	0.75,	0.5,	0,	0,	0,	0

" > $dir/tests/A15TiNb3/evsv.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/A15TiNb3/meamz_evsv.out

echo "#v/v0, P, a" > $dir/tests/A15TiNb3/pvsv.dat
awk -v a0=$A15TiNb3mmz -v np=$numpts -v pc=$PCONV 'NR>1{
					   pum=0;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; pum+=$5;
					   getline; getline;
					   pum=pc*pum/3;

					   vol = (0.80 + ($1/np)*0.4);
					   a = a0*(vol)^(1./3);
					   vol = a0*a0*a0*vol
					   print vol/8, pum, a;
				     }' data.stress >> $dir/tests/A15TiNb3/pvsv.dat 

echo "#v, a" > $dir/tests/A15TiNb3/evsv.dat
awk -v a0=$A15TiNb3mmz -v np=$numpts 'NR>1{
					   getline;
					   vol = (0.80 + ($1/np)*0.4)*a0*a0*a0/8.;
					   print vol, $6;
				     }' data.energy >> $dir/tests/A15TiNb3/evsv.dat 
fi # lammps_evsv_flag A15TiNb3

sed -i '1d' $dir/tests/A15TiNb3/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/A15TiNb3/pvsv.dat | sort -k2,2 -g | head -1`
	A15TiNb3PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${A15TiNb3PLAT[0]} $A15TiNb3eqp"
		#A15TiNb3eqp=`echo $line | awk '{print $3}'`
		#A15TiNb3vop=`echo $line | awk '{print $1}'`
		#A15TiNb3pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
A15TiNb3PLAT["0"]=$A15TiNb3eqp

awk -v Voo=$A15TiNb3vop '{print $1/Voo, $2, $3}' $dir/tests/A15TiNb3/pvsv.dat > tmp; mv tmp $dir/tests/A15TiNb3/pvsv.dat

# ---------------------- elastic constants for A15 TiNb3----------------------------
echo -e "${prp}Elastic constants:${non}"
echo "# pressure, c11, c12, c44" > $dir/tests/A15TiNb3/C_VS_P.dat
for PRESS in 0; do 
echo "$PRESS GPa..."
for jj in $(seq 1 7); do

printf "\t $jj "

if [ "$lammps_elcon_flag" == "TRUE" ]; then

P=`echo $PRESS*10000 | bc -l`
case $jj in
	1) e1="(1+v_d)";	    e2="(1+v_d)";		e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	2) e1="(1+v_d)";	    e2="(1-v_d)";		e3="(1+v_d2/(1-v_d2))";	e4=0;	  e5=0;	    e6=0	;;
	3) e1="(1+v_d2/(1-v_d2))";  e2="(1+v_d)";		e3="(1-v_d)";		e4=0;	  e5=0;	    e6=0	;;
	4) e1="(1-v_d)";	    e2="(1+v_d2/(1-v_d2))";	e3="(1+v_d)";		e4=0;	  e5=0;	    e6=0	;;
	5) e1="(1+v_d2/(4-v_d2))";  e2="1";			e3="1";			e4="v_d"; e5=0;	    e6=0	;;
	6) e1="1";		    e2="(1+v_d2/(4-v_d2))";	e3="1";			e4=0;	  e5="v_d"; e6=0	;;
	7) e1="1";		    e2="1";			e3="(1+v_d2/(4-v_d2))";	e4=0;	  e5=0;	    e6="v_d"	;;
esac
#case $jj in
#	1) XD="1+v_d"; YD="sqrt(3)*(1+v_d)/2"; ZD="1+v_d"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	2) XD="1+v_d"; YD="sqrt(3)*(1-v_d)/2"; ZD="1/(1-v_d2)"; XYD="-0.5*(1+v_d)"; XZD="0"; YZD="0" ;;
#	3) XD="1/(1-v_d2)"; YD="sqrt(3)*(1+v_d)/2"; ZD="1-v_d"; XYD="-0.5/(1-v_d2)"; XZD="0"; YZD="0" ;;
#	4) XD="1-v_d"; YD="sqrt(3)/(2-2*v_d2)"; ZD="1+v_d"; XYD="-0.5*(1-v_d)"; XZD="0"; YZD="0" ;;
#	5) XD="1/(1-(v_d2)/4)"; YD="sqrt(3)/2"; ZD="1"; XYD="-0.5/(1-(v_d2)/4)"; XZD="0"; YZD="v_d" ;;
#	6) XD="1"; YD="sqrt(3)/((1-(v_d2)/4)*2)"; ZD="1"; XYD="-0.5"; XZD="v_d"; YZD="0" ;;
#	7) XD="1"; YD="sqrt(3)/2"; ZD="1/(1-(v_d2)/4)"; XYD="-0.5+v_d"; XZD="0"; YZD="0" ;;
#esac

cat > $dir/tests/A15TiNb3/elcon.lin << !!
###############################################################
# for use in script looping over the seven strains of Trinkle #
###############################################################

units metal
atom_style atomic

variable boxa equal ${A15TiNb3PLAT["$PRESS"]}
variable boxb equal ${A15TiNb3PLAT["$PRESS"]}
variable boxc equal ${A15TiNb3PLAT["$PRESS"]}
variable boxxy equal 0
variable boxxz equal 0
variable boxyz equal 0

lattice		custom ${A15TiNb3PLAT["$PRESS"]} &
		a1 1.0 0.0 0.0 &
		a2 0.0 1.0 0.0 &
		a3 0.0 0.0 1.0 &
		basis 0.0	0.0	0.00 &
		basis 0.5	0.5	0.50 &
		basis 0.5	0.0	0.25 &
		basis 0.5	0.0	0.75 &
		basis 0.0	0.25	0.50 &
		basis 0.0	0.75	0.50 &
		basis 0.25	0.50	0.00 & 
		basis 0.75	0.50	0.00
		

region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} 0 0 0 units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} & 
		basis 5 ${idx["Nb"]} basis 6 ${idx["Nb"]} basis 7 ${idx["Nb"]} basis 8 ${idx["Nb"]}
 

dump D all xyz 1 elcon$jj.xyz
 
mass 1 $mass1
mass 2 $mass2

# variables for loop
variable dmax equal $ecstr/100	# strain percent (max is half of this)
variable jmax equal 100		# number of steps
variable conv equal 160.217656  # GPa per eV/A^3

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# computes
compute strs all stress/atom NULL
compute sigxx all reduce sum c_strs[1]
compute sigyy all reduce sum c_strs[2]
compute sigzz all reduce sum c_strs[3]
compute sigxy all reduce sum c_strs[4]
compute sigxz all reduce sum c_strs[5]
compute sigyz all reduce sum c_strs[6]

#variable tmp equal lx
#variable lxi equal \${tmp}
#
timestep 0.001
#fix rel all press/berendsen iso $P $P 1.0
#run 1000
#unfix rel
#
#variable tmp equal lx/v_lxi
#variable sca equal \${tmp}
timestep 0.01
fix rel all box/relax iso $P fixedpoint 0 0 0
min_style hftn
min_modify line forcezero
minimize 0.0 1e-4 100000 100000
unfix rel

# initialization
thermo 1
thermo_style custom pe vol c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz
run 0
variable tmp equal pe
variable e0 equal \${tmp}
variable ndef equal 10
variable tmp equal c_sigxx
variable Sxx0 equal \${tmp}

variable tmp equal c_sigyy
variable Syy0 equal \${tmp}

variable tmp equal c_sigzz
variable Szz0 equal \${tmp}

variable tmp equal c_sigxy
variable Sxy0 equal \${tmp}

variable tmp equal c_sigxz
variable Sxz0 equal \${tmp}

variable tmp equal c_sigyz
variable Syz0 equal \${tmp}

variable tmp equal lx
variable LX equal \${tmp}
variable tmp equal ly
variable LY equal \${tmp}
variable tmp equal lz
variable LZ equal \${tmp}
variable XY equal xy
variable lat equal \${tmp}
variable YZ equal yz
variable lat equal \${tmp}
variable XZ equal xz
variable lat equal \${tmp}

# variables
variable Sxx equal (c_sigxx-\${Sxx0})/vol/10000
variable Syy equal (c_sigyy-\${Syy0})/vol/10000
variable Szz equal (c_sigzz-\${Szz0})/vol/10000
variable Sxy equal (c_sigxy-\${Sxy0})/vol/10000
variable Sxz equal (c_sigxz-\${Sxz0})/vol/10000
variable Syz equal (c_sigyz-\${Syz0})/vol/10000

variable tmp equal lx
variable lx0 equal \${tmp}
variable tmp equal ly
variable ly0 equal \${tmp}
variable tmp equal lz
variable lz0 equal \${tmp}
variable tmp equal xy
variable xy0 equal \${tmp}
variable tmp equal xz
variable xz0 equal \${tmp}
variable tmp equal yz
variable yz0 equal \${tmp}

# new thermo
thermo 10
thermo_style custom step pe vol lx ly lz c_sigxx c_sigyy c_sigzz c_sigxy c_sigxz c_sigyz v_Sxx v_Syy v_Szz v_Sxy v_Sxz v_Syz

# fixes
fix P all print 1 "\${d} \${Sxx} \${Syy} \${Szz} \${Syz} \${Sxz} \${Sxy}" file $dir/tests/A15TiNb3/stresses$jj.dat	# printed in voigt notation 1->2->3->4->5->6

# loop:
variable j loop 0 \${jmax}
label loop
variable d equal v_dmax*((v_j)/(v_jmax)-1/2)
variable d2 equal (v_d*v_d)

variable xd  equal ($e1*\${lx0})
variable yd  equal ($e2*\${ly0}+$e6*\${xy0})
variable zd  equal ($e3*\${lz0}+$e5*\${xz0}+$e4*\${yz0})
variable xyd equal ($e1*\${xy0}+$e6*\${ly0})
variable xzd equal ($e1*\${xz0}+$e6*\${yz0}+$e5*\${lz0}) 
variable yzd equal ($e6*\${xz0}+$e2*\${yz0}+$e4*\${lz0})

change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} xz final \${xzd} yz final \${yzd} remap units box

#min_style cg
#minimize 0.0 1e-10 100 1000

run 1

next j
jump $dir/tests/A15TiNb3/elcon.lin loop 
!!

$lammps < $dir/tests/A15TiNb3/elcon.lin > $dir/tests/A15TiNb3/elcon.out
sed -i '1d' $dir/tests/A15TiNb3/stresses$jj.dat

else # else stress strain-curves with meamz

cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/A15TiNb3/elcon$jj.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

DELPT=5
strain=0.002
rm -f $dir/tests/A15TiNb3/elcon$jj.conf
for i in $(seq -$DELPT $DELPT); do
	omgmmz=`echo ${omgPLAT["$PRESS"]} | bc -l`
	del=`echo "$strain*($i/$DELPT)"	| bc -l`
	python $ecgen $jj $del A15TiNb3 $omgmmz conf >> $dir/tests/A15TiNb3/elcon$jj.conf 

done

$meamz -p $dir/tests/meamz_params > $dir/tests/omega/meamz_elcon.out
awk -v e0=$strain -v np=$DELPT -v pc=$PCONV 'NR>2{
					
					del=(($1-np)/np)*e0
					sxx=-pc*$5; getline	
					syy=-pc*$5; getline	
					szz=-pc*$5; getline	
					sxy=-pc*$5; getline	
					syz=-pc*$5; getline	
					szx=-pc*$5; 	
					print del, sxx, syy, szz, syz, szx, sxy
				     }' data.stress >> $dir/tests/A15TiNb3/stresses$jj.dat 

fi	# lammps elcon flag

done
echo ""

rm -f $dir/tests/A15TiNb3/EC_fits.dat; touch $dir/tests/A15TiNb3/EC_fits.dat

# first row fits
awk '{print $1, $2}' $dir/tests/A15TiNb3/stresses1.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/A15TiNb3/stresses1.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/A15TiNb3/stresses1.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;

# second row fits
awk '{print $1, $2}' $dir/tests/A15TiNb3/stresses2.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/A15TiNb3/stresses2.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/A15TiNb3/stresses2.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;

# third row fits
awk '{print $1, $2}' $dir/tests/A15TiNb3/stresses3.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/A15TiNb3/stresses3.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/A15TiNb3/stresses3.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;

# fourth row fits
awk '{print $1, $2}' $dir/tests/A15TiNb3/stresses4.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $3}' $dir/tests/A15TiNb3/stresses4.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;
awk '{print $1, $4}' $dir/tests/A15TiNb3/stresses4.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;

# fifth row fit
awk '{print $1, $5}' $dir/tests/A15TiNb3/stresses5.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;

# sixth row fit
awk '{print $1, $6}' $dir/tests/A15TiNb3/stresses6.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;

# seventh row fit
awk '{print $1, $7}' $dir/tests/A15TiNb3/stresses7.dat > tmp; python $fit tmp >> $dir/tests/A15TiNb3/EC_fits.dat;

# now decouple!
declare -a coup=( `cat $dir/tests/A15TiNb3/EC_fits.dat` )
A15TiNb3C11i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]}-3*${coup[9]})/3"`
A15TiNb3C12i=`python -c "print (${coup[2]}+2*${coup[5]}+3*${coup[6]}+ ${coup[8]})/3"`
A15TiNb3C13i=`python -c "print (${coup[2]}+2*${coup[5]}+${coup[8]})/3"`
A15TiNb3C22i=`python -c "print (${coup[2]}-${coup[5]}+3*${coup[7]}+${coup[8]})/3"`
A15TiNb3C23i=`python -c "print (${coup[2]}-${coup[5]}+${coup[8]})/3"`
A15TiNb3C33i=`python -c "print (${coup[2]}-${coup[5]}-2*${coup[8]})/3"`
A15TiNb3C44i=${coup[12]}
A15TiNb3C55i=${coup[13]}
A15TiNb3C66i=${coup[14]}

A15TiNb3C11p=`python -c "print ($A15TiNb3C11i+$A15TiNb3C22i+$A15TiNb3C33i)/3"`
A15TiNb3C13p=`python -c "print ($A15TiNb3C13i+$A15TiNb3C23i+$A15TiNb3C13i)/3"`
A15TiNb3C44p=`python -c "print ($A15TiNb3C44i+$A15TiNb3C55i+$A15TiNb3C66i)/3"`

echo $PRESS $A15TiNb3C11p $A15TiNb3C12i $A15TiNb3C13p $A15TiNb3C33i $A15TiNb3C44p >> $dir/tests/A15TiNb3/C_VS_P.dat

if [ "$PRESS" == "0" ]; then
	A15TiNb3C11p=`printf '%3.f' $A15TiNb3C11p`
	A15TiNb3C12p=`printf '%3.f' $A15TiNb3C12p`
	A15TiNb3C44p=`printf '%3.f' $A15TiNb3C44p`

	C11A15TiNb3d=`printf '%3.f' $C11A15TiNb3d`
	C12A15TiNb3d=`printf '%3.f' $C12A15TiNb3d`
	C44A15TiNb3d=`printf '%3.f' $C44A15TiNb3d`
fi
done

## PHONONS FOR A15 TiNb3
if [ "$phonon_flag" == "TRUE" ]; then
echo "phonons..."

cat > $dir/tests/A15TiNb3/OPT.POSCAR << -
A15 TiNb3
$A15TiNb3eqp
1.00	0.00	0.00
0.00	1.00	0.00
0.00	0.00	1.00
Ti Nb
2 6
Direct
0.000000000000	0.000000000000	0.000000000000
0.500000000000	0.500000000000	0.500000000000
0.500000000000	0.000000000000	0.250000000000
0.500000000000	0.000000000000	0.750000000000
0.000000000000	0.250000000000	0.500000000000
0.000000000000	0.750000000000	0.500000000000
0.250000000000	0.500000000000	0.000000000000
0.750000000000	0.500000000000	0.000000000000
-

if [ "$typ" == "GMEAM" ]; then
	PS="gmeam/spline"
else
	PS="meam/alloy/spline"
fi 

cat > $dir/tests/A15TiNb3/phonons.py << !
#
# script using ASE to compute phonons
#

from ase.lattice import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons
from ase.calculators.lammpsrun import LAMMPS
from ase.units import _hbar, _e
from ase.io import read
import numpy as np

# set lammps calculator parameters ('dictionary' data type)
ps = "$PS"
pc = ["* * $dir/lammps.pt $elem1 $elem2"]
ms = ["1 $mass1", "2 $mass2"]
so=['$elem1','$elem2']

params = dict(pair_style=ps, pair_coeff=pc, mass=ms)
calc = LAMMPS(parameters=params, specorder=so)

## Setup crystal and EMT calculator
atoms = read('$dir/tests/A15TiNb3/OPT.POSCAR')

atoms.set_calculator(calc)

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.001)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True, method='standard', symmetrize=5)

G = [0, 0, 0]
X = [1./2, 0, 0]
M = [1./2, 1./2, 0]
R = [1./2, 1./2, 1./2]


point_names = ['\$R\$', '\$\Gamma\$', '\$M\$', '\$X\$', '\$\Gamma\$']
path = [R, G, M, X, G]
dirs = ['\$[\\\xi\\\xi\\\xi]\$', '\$[\\\xi\\\xi0]\$', '\$[0\\\bar{\\\xi}0]\$', '\$[\\\xi00]\$']


# Band structure in THz
conv = 241.79893	# eV to THz
path_kc, q, Q = get_bandpath(path, atoms.cell, 1000)
omega_kn = conv * ph.band_structure(path_kc, verbose=False)

# Calculate phonon DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= conv
dos_e /= conv

# directions
dirQ = np.array([])
for i in range(0,np.size(Q)-1):
	dirQ = np.append(dirQ, (Q[i+1]+Q[i])/2)


# Plot the band structure and DOS
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])

dft = np.loadtxt('$dftdat/Ti-Nb/A15TiNb3/phonons.dat')
b = 2*np.pi/5.25368042

max_band = 0
min_band = 0
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    omega_nd= dft[:, n+1]
    max_this = np.max(omega_n)
    min_this = np.min(omega_n)
    max_thisd= np.max(omega_nd)
    min_thisd= np.min(omega_nd)
    max_band = np.max([max_band, max_this, max_thisd])
    min_band = np.min([min_band, min_this, min_thisd])
    plt.plot(q, omega_n, 'k-', lw=2)
    plt.plot(b*dft[:,0], omega_nd, color='gray', linestyle='--', lw=2)

max_band *= 1.05 # max band >= 0
min_band *= 1.05 # min band <= 0
plt.title('A15 TiNb\$_{3}\$')
plt.xticks(Q, point_names, fontsize=18)
for i in range(0,np.size(dirQ)):
	plt.text(dirQ[i], min_band-0.02*max_band, dirs[i], fontsize=15, ha='center', va='top')
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(min_band, max_band)
plt.ylabel("Frequency ($\mathrm{THz}$)", fontsize=18)
plt.grid('on')
plt.axes([.771, .07, .17, .85])
plt.fill_between(np.absolute(dos_e), omega_e, y2=0, color='gray', edgecolor='k', lw=1)
plt.ylim(min_band, max_band)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("\$DOS\$", fontsize=18)
plt.savefig('$dir/tests/A15TiNb3_phonons.png')
ph.clean
!

python $dir/tests/A15TiNb3/phonons.py > $dir/tests/A15TiNb3/phonon_log
rm -f *.pckl
fi


# omega TiNb2
#<><><><><><><>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< omega TiNb2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t omega-TiNb2... ${non}"

cat > $dir/tests/omgTiNb2/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid
variable boxa equal $omgTiNb2lat
variable boxb equal $omgTiNb2lat*sqrt(3)/2
variable boxc equal $omgTiNb2lat*$omgTiNb2coa
variable boxxy equal $omgTiNb2lat*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $omgTiNb2lat a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $omgTiNb2coa &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large

create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0 scalexy yes scalexz yes scaleyz yes
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/omgTiNb2/eq.in > $dir/tests/omgTiNb2/eq.out

omgTiNb2eqp=`tail -1 $dir/tests/omgTiNb2/eq.out | awk '{print $1}'`
omgTiNb2coap=`tail -1 $dir/tests/omgTiNb2/eq.out | awk '{print $2}'`
omgTiNb2vop=`tail -1 $dir/tests/omgTiNb2/eq.out | awk '{print $3}'`
omgTiNb2pote=`tail -1 $dir/tests/omgTiNb2/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for OMGti -----------------------------------------

cat > $dir/tests/omgTiNb2/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR OMGti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

#initialization variables
variable boxa equal $omgTiNb2eqp
variable boxb equal $omgTiNb2eqp*sqrt(3)/2
variable boxc equal $omgTiNb2eqp*$omgTiNb2coap
variable boxxy equal $omgTiNb2eqp*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $omgTiNb2eqp a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $omgTiNb2coap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/omgTiNb2/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/omgTiNb2/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd  equal  \${boxa}*(1-v_dmax/2)
variable xf  equal  \${boxa}*(1+v_dmax/2)
variable yd  equal  \${boxb}*(1-v_dmax/2)
variable yf  equal  \${boxb}*(1+v_dmax/2)
variable zd  equal  \${boxc}*(1-v_dmax/2)
variable zf  equal  \${boxc}*(1+v_dmax/2)
variable xyd equal -0.5*\${xd}
variable xyf equal -0.5*\${xf}
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/omgTiNb2/evsv.in loop
!!
$lammps < $dir/tests/omgTiNb2/evsv.in > $dir/tests/omgTiNb2/evsv.out

sed -i '1d' $dir/tests/omgTiNb2/evsv.dat
line=`python $BMFit $dir/tests/omgTiNb2/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/omgTiNb2/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	omgTiNb2bulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	omgTiNb2bulkp=0
fi

sed -i '1d' $dir/tests/omgTiNb2/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/omgTiNb2/pvsv.dat | sort -k2,2 -g | head -1`
	omgTiNb2PLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${omgTiNb2PLAT[0]} $omgTiNb2eqp"
		#omgTiNb2eqp=`echo $line | awk '{print $3}'`
		#omgTiNb2vop=`echo $line | awk '{print $1}'`
		#omgTiNb2pote=`echo $line | awk '{print 1000*$4}'`
	fi
done
omgTiNb2PLAT["0"]=$omgTiNb2eqp
awk -v Voo=$omgTiNb2vop '{print $1/Voo, $2, $3}' $dir/tests/omgTiNb2/pvsv.dat > tmp; mv tmp $dir/tests/omgTiNb2/pvsv.dat

# ---------------------- pressure dependence of lattice constants for omg tinb2 ----------
echo -e "${prp}Pressure dependence of lattice constants...${non}"

rm -f $dir/tests/omgTiNb2/abcvp.dat; echo "#a, c/a, gamma" > $dir/tests/omgTiNb2/abcvp.dat
for PRESS in $PRESSLAT; do
PRES=`echo "10000*$PRESS" | bc -l`
cat > $dir/tests/omgTiNb2/abcvp.lin << __
## LAMMPS script
units metal
boundary p p p
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable a1	equal	1.00
variable a2	equal	sqrt(3)
variable a3	equal	$omgTiNb2coap
variable boxa equal $omgTiNb2eqp*\${a1}
variable boxb equal $omgTiNb2eqp*\${a2}
variable boxc equal $omgTiNb2eqp*\${a3}

lattice custom $omgTiNb2eqp &
	a1 \${a1}	0.0	0.0 &
	a2 0.0		\${a2}	0.0 &
	a3 0.0		0.0	\${a3} &
	basis 0.0	0.0	0.0 basis 0.5	0.5	0.0 &
	basis 0.0	0.33333 0.5 basis 0.5   0.83333 0.5 &
	basis 0.0	0.66667 0.5 basis 0.5   0.16667 0.5
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Nb"]} basis 4 ${idx["Nb"]} basis 5 ${idx["Nb"]} basis 6 ${idx["Nb"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#region frz1 sphere 0.0 0.0 0.0 0.1 units box
#region frz2 sphere 0.5 0.5 0.0 0.1 units box
#region frz union 2 frz1 frz2
#group frz region frz
#
#fix FREEZE freeze frz

thermo 1000
thermo_style custom step temp press pe ke etotal lx ly lz
timestep 0.001

variable a equal lx
variable coa equal lz/lx
variable Gamma equal 180-2*180*atan(ly/lx)/3.141592654

fix P all press/berendsen aniso $PRES $PRES 1000.0
run 10000

dump D all custom 1 $dir/tests/omgTiNb2/abcvp.coords xs ys zs

print "$PRESS \${a} \${coa} \${Gamma}"
__

$lammps < $dir/tests/omgTiNb2/abcvp.lin > $dir/tests/omgTiNb2/abcvp.out
tail -1 $dir/tests/omgTiNb2/abcvp.out >> $dir/tests/omgTiNb2/abcvp.dat

done





# omega Ti2Nb
# <><><><><><
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< omega Ti2Nb >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo -e "${yel} \t omega-Ti2Nb... ${non}"

cat > $dir/tests/omgTi2Nb/eq.in << !!
################################################
#	CALCULATES EQUILIBRIUM HCP-Ti LATTICE
#	PARAMETER
################################################

units		metal
atom_style	atomic

#define simulation region and bcc grid
variable boxa equal $omgTi2Nblat
variable boxb equal $omgTi2Nblat*sqrt(3)/2
variable boxc equal $omgTi2Nblat*$omgTi2Nbcoa
variable boxxy equal $omgTi2Nblat*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $omgTi2Nblat a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $omgTi2Nbcoa &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large

create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box basis 1 ${idx["Nb"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#set up thermo style
thermo_style custom step etotal pe ke vol press temp lx ly lz
thermo 1

# minimize
fix 1 all box/relax x 0.0 y 0.0 z 0.0 couple xy fixedpoint 0 0 0 scalexy yes scalexz yes scaleyz yes
min_style	cg
minimize	$META_RELAX_ETOL 0.0 10000 1000000

variable	coa equal lz/lx 
variable	a equal lx
variable	vat equal vol/atoms
variable	spe equal pe/atoms

print '\${a} \${coa} \${vat} \${spe}'
!!

$lammps < $dir/tests/omgTi2Nb/eq.in > $dir/tests/omgTi2Nb/eq.out

omgTi2Nbeqp=`tail -1 $dir/tests/omgTi2Nb/eq.out | awk '{print $1}'`
omgTi2Nbcoap=`tail -1 $dir/tests/omgTi2Nb/eq.out | awk '{print $2}'`
omgTi2Nbvop=`tail -1 $dir/tests/omgTi2Nb/eq.out | awk '{print $3}'`
omgTi2Nbpote=`tail -1 $dir/tests/omgTi2Nb/eq.out | awk '{print $4}'`


##--------------------------------------------------------------------------------------------
echo -e "${prp}E-V curve...${non}"
#----------------------------- energy volume for OMGti -----------------------------------------

cat > $dir/tests/omgTi2Nb/evsv.in << !!
################################################
#  CALCULATES ENERGY VERSUS VOLUME CURVE
# FOR OMGti LATTICE
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

#initialization variables
variable boxa equal $omgTi2Nbeqp
variable boxb equal $omgTi2Nbeqp*sqrt(3)/2
variable boxc equal $omgTi2Nbeqp*$omgTi2Nbcoap
variable boxxy equal $omgTi2Nbeqp*(-0.5)
variable boxxz equal 0
variable boxyz equal 0

lattice		custom $omgTi2Nbeqp a1 1.0 0.0 0.0 a2 -0.5 0.86602540378 0.0 a3 0.0 0.0 $omgTi2Nbcoap &
		basis 0.0 0.0 0.0 basis 0.6666666 0.3333333 0.5 basis 0.3333333 0.6666666 0.5
region mybox prism 0 \${boxa} 0 \${boxb} 0 \${boxc} \${boxxy} 0 0 units box
box tilt large
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box basis 1 ${idx["Nb"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#initilization variables
thermo_style custom step pe etotal press vol lx ly lz pxx pyy pzz pxy pxz pyz
thermo 1
timestep 0.001
run 1
variable Epa equal pe/atoms	#energy and volume PER ATOM
variable Vpa equal vol/atoms           #
variable a equal lx
variable prz equal press/10000

fix P all print 1 "\${Vpa} \${Epa}" file $dir/tests/omgTi2Nb/evsv.dat screen no title "# V/atom | Energy"
fix P2 all print 1 "\${Vpa} \${prz} \$a \${Epa}" file $dir/tests/omgTi2Nb/pvsv.dat screen no title "# V/atom | pressure | lattice constant"

variable xd  equal  \${boxa}*(1-v_dmax/2)
variable xf  equal  \${boxa}*(1+v_dmax/2)
variable yd  equal  \${boxb}*(1-v_dmax/2)
variable yf  equal  \${boxb}*(1+v_dmax/2)
variable zd  equal  \${boxc}*(1-v_dmax/2)
variable zf  equal  \${boxc}*(1+v_dmax/2)
variable xyd equal -0.5*\${xd}
variable xyf equal -0.5*\${xf}
change_box all x final 0 \${xd} y final 0 \${yd} z final 0 \${zd} xy final \${xyd} remap units box

reset_timestep 0
fix def all deform 1 x final 0 \${xf} y final 0 \${yf} z final 0 \${zf} xy final \${xyf} units box

run \${jmax}

#variable j loop 0 \${jmax}
#label loop
#min_style fire
#minimize 1e-5 0.0 100 1000
#run 1
#next j
#jump $dir/tests/omgTi2Nb/evsv.in loop
!!
$lammps < $dir/tests/omgTi2Nb/evsv.in > $dir/tests/omgTi2Nb/evsv.out

sed -i '1d' $dir/tests/omgTi2Nb/evsv.dat
line=`python $BMFit $dir/tests/omgTi2Nb/evsv.dat`
#line=`echo "data=Import[\"$dir/tests/omgTi2Nb/evsv.dat\"];
#data=Take[data,{$MDIN,$MDOU}];
#fit=FindFit[data, Eo + (9/8)*B*Vo*((Vo/V)^(2/3) - 1)^2 + (9/16)*B*Vo*(Bp-4)*((Vo/V)^(2/3) - 1)^3, {Eo,{Vo,18},Bp,{B,100}}, V];
#B=160.256*Replace[B,fit[[4]]];
#Eo=Replace[Eo,fit[[1]]];
#Vo=Replace[Vo,fit[[2]]];
#Print[B]
#Print[Vo]
#Print[Eo]
#" | math -noprompt`
##echo $line

var=`echo $line | awk '{print int($1)}'`
if [ "$var" -ne "0" ] 2>/dev/null; then
	echo -e "\t ${dbl} B-M fit successful! ${non}"
	omgTi2Nbbulkp=`echo $line | awk '{print $1}'`
else
	echo -e "\t ${red} B-M fit unsuccessful... ${non}"
	omgTi2Nbbulkp=0
fi

sed -i '1d' $dir/tests/omgTi2Nb/pvsv.dat

for pressure in 0 10 20 25 30 40 50 60 70 75 80 90 100; do
	line=`awk -v p=$pressure '{print $1, ($2-p)**2, $3, $4}' $dir/tests/omgTi2Nb/pvsv.dat | sort -k2,2 -g | head -1`
	omgTi2NbPLAT["$pressure"]=`echo $line | awk '{print $3}'`
	
	if [ "$pressure" == "0" ]; then
		echo -e "\\t ${omgTi2NbPLAT[0]} $omgTi2Nbeqp"
		#omgTi2Nbeqp=`echo $line | awk '{print $3}'`
		#omgTi2Nbvop=`echo $line | awk '{print $1}'`
		#omgTi2Nbpote=`echo $line | awk '{print 1000*$4}'`
	fi
done
omgTi2NbPLAT["0"]=$omgTi2Nbeqp
awk -v Voo=$omgTi2Nbvop '{print $1/Voo, $2, $3}' $dir/tests/omgTi2Nb/pvsv.dat > tmp; mv tmp $dir/tests/omgTi2Nb/pvsv.dat

# ---------------------- pressure dependence of lattice constants for omg tinb2 ----------
echo -e "${prp}Pressure dependence of lattice constants...${non}"

rm -f $dir/tests/omgTi2Nb/abcvp.dat; echo "#a, c/a, gamma" > $dir/tests/omgTi2Nb/abcvp.dat
for PRESS in $PRESSLAT; do
PRES=`echo "10000*$PRESS" | bc -l`
cat > $dir/tests/omgTi2Nb/abcvp.lin << __
## LAMMPS script
units metal
boundary p p p
atom_style	atomic

#initialization variables
variable	dmax equal $stpct/100
variable	jmax equal $nevpt

variable a1	equal	1.00
variable a2	equal	sqrt(3)
variable a3	equal	$omgTi2Nbcoap
variable boxa equal $omgTi2Nbeqp*\${a1}
variable boxb equal $omgTi2Nbeqp*\${a2}
variable boxc equal $omgTi2Nbeqp*\${a3}

lattice custom $omgTi2Nbeqp &
	a1 \${a1}	0.0	0.0 &
	a2 0.0		\${a2}	0.0 &
	a3 0.0		0.0	\${a3} &
	basis 0.0	0.0	0.0 basis 0.5	0.5	0.0 &
	basis 0.0	0.33333 0.5 basis 0.5   0.83333 0.5 &
	basis 0.0	0.66667 0.5 basis 0.5   0.16667 0.5
		
region		mybox block 0 \${boxa} 0 \${boxb} 0 \${boxc} units box
create_box	2 mybox

#create atoms
create_atoms 	${idx["Ti"]} box basis 1 ${idx["Nb"]} basis 2 ${idx["Nb"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]}
 
mass		1 $mass1 
mass		2 $mass2

group		grp region mybox

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

#region frz1 sphere 0.0 0.0 0.0 0.1 units box
#region frz2 sphere 0.5 0.5 0.0 0.1 units box
#region frz union 2 frz1 frz2
#group frz region frz
#
#fix FREEZE freeze frz

thermo 1000
thermo_style custom step temp press pe ke etotal lx ly lz
timestep 0.001

variable a equal lx
variable coa equal lz/lx
variable Gamma equal 180-2*180*atan(ly/lx)/3.141592654

fix P all press/berendsen aniso $PRES $PRES 1000.0
run 10000

dump D all custom 1 $dir/tests/omgTi2Nb/abcvp.coords xs ys zs

print "$PRESS \${a} \${coa} \${Gamma}"
__

$lammps < $dir/tests/omgTi2Nb/abcvp.lin > $dir/tests/omgTi2Nb/abcvp.out
tail -1 $dir/tests/omgTi2Nb/abcvp.out >> $dir/tests/omgTi2Nb/abcvp.dat

done
# <><>><><><><---------------------------------------------------------

########################################################################################
# 	THERMODYNAMIC PROPERTIES
########################################################################################
if [ "$thermo_flag" == "TRUE" ]; then
echo -e "${grn}THERMODYNAMICS...${non}"
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BCC formation energy >>>>>>>>>>>>>>>>>>>
echo -e "BCC formation energy..."
PRESS=1.0 # 1 atm
# first calculate BCC Ti energy and enthalpy at given pressure
cat > $dir/tests/thermo/BCC_formation.lin << ++
# computes formation energy for solid solution BCC ti-nb compound
units metal
atom_style atomic
boundary p p p

variable SIZE equal $B2eqp

lattice bcc $B2eqp

region mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box 2 mybox
create_atoms ${idx["Ti"]} box

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

thermo 1
thermo_style custom step temp pe ke etotal enthalpy
run 0

variable peat equal pe/atoms 
variable heat equal enthalpy/atoms

thermo_style custom step temp etotal enthalpy v_peat v_heat

fix REL all box/relax iso $PRESS
min_style cg
minimize 1e-9 0.0 1000 100000

run 0
print "\${peat} \${heat}"
++

$lammps < $dir/tests/thermo/BCC_formation.lin > $dir/tests/thermo/BCC_formation.out
line=`tail -1 $dir/tests/thermo/BCC_formation.out`
ETi=`echo $line | awk '{print $1}'`
HTi=`echo $line | awk '{print $2}'`

cat > $dir/tests/thermo/BCC_formation.lin << ++
# computes formation energy for solid solution BCC ti-nb compound
units metal
atom_style atomic
boundary p p p

variable SIZE equal $B2eqp

lattice bcc $B2eqp

region mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box 2 mybox
create_atoms ${idx["Nb"]} box

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

thermo 1
thermo_style custom step temp pe ke etotal enthalpy
run 0

variable peat equal pe/atoms 
variable heat equal enthalpy/atoms

thermo_style custom step temp etotal enthalpy v_peat v_heat

fix REL all box/relax iso $PRESS
min_style cg
minimize 1e-9 0.0 1000 10000

run 0
print "\${peat} \${heat}"
++

$lammps < $dir/tests/thermo/BCC_formation.lin > $dir/tests/thermo/BCC_formation.out
line=`tail -1 $dir/tests/thermo/BCC_formation.out`
ENb=`echo $line | awk '{print $1}'`
HNb=`echo $line | awk '{print $2}'`

SIZE=1
NTIMES=10
CONC_PCTS=$(seq 5 5 95)

echo "# Concentration, formation energy" > $dir/tests/thermo/BCC_formation.dat
echo "0.0 0.00 0.00 0.00 0.00" >> $dir/tests/thermo/BCC_formation.dat
for CONC in $CONC_PCTS; do
PCT=`echo "$CONC/100" | bc -l`

rm -f tmp.dat; touch tmp.dat
for i in $(seq 1 $NTIMES); do

cat > $dir/tests/thermo/BCC_formation.lin << ++
# computes formation energy for solid solution BCC ti-nb compound
units metal
atom_style atomic
boundary p p p

variable SIZE equal $SIZE*$B2eqp
variable ETi  equal $ETi
variable ENb  equal $ENb
variable HTi  equal $HTi
variable HNb  equal $HNb

lattice bcc $B2eqp

region mybox block 0 \${SIZE} 0 \${SIZE} 0 \${SIZE} units box
create_box 2 mybox
create_atoms ${idx["Ti"]} box

set region mybox type/fraction 2 $PCT $RANDOM

group TI type 1
group NB type 2

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

thermo 1
thermo_style custom step temp pe ke etotal enthalpy press
run 0

variable Eform equal 1000*(pe-count(TI)*\${ETi}-count(NB)*\${ENb})/count(all)
variable Hform equal 1000*(pe-count(TI)*\${HTi}-count(NB)*\${HNb})/count(all)
# Gform uses solid solution entropy

thermo_style custom step temp pe ke etotal enthalpy press v_Eform v_Hform

fix REL all box/relax iso $PRESS
min_style cg
minimize 1e-9 0.0 1000 10000

run 0
print "\${Eform} \${Hform}"
++
$lammps < $dir/tests/thermo/BCC_formation.lin > $dir/tests/thermo/BCC_formation.out
tail -1 $dir/tests/thermo/BCC_formation.out >> tmp.dat

done

Eavg=`awk '{sum += $1; sumsq += $1*$1} END {if (NR>0) print sum/NR, sqrt((sumsq-sum*sum/NR)/NR)}' tmp.dat`
Havg=`awk '{sum += $2; sumsq += $2*$2} END {if (NR>0) print sum/NR, sqrt((sumsq-sum*sum/NR)/NR)}' tmp.dat`
rm tmp.dat

echo "$PCT $Eavg $Havg" >> $dir/tests/thermo/BCC_formation.dat 
done
echo "1.00 0.00 0.00 0.00 0.00" >> $dir/tests/thermo/BCC_formation.dat

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< HCP formation energy >>>>>>>>>>>>>>>>>>>
echo -e "HCP formation energy..."
PRESS=1.0 # 1 atm
# first calculate HCP Ti energy and enthalpy at given pressure
cat > $dir/tests/thermo/HCP_formation.lin << ++
# computes formation energy for solid solution HCP ti-nb compound
units metal
atom_style atomic
boundary p p p

variable SIZE equal $HCPtieqp
variable YIZE equal sqrt(3)*$HCPtieqp
variable ZIZE equal sqrt(8/3)*$HCPtieqp

lattice hcp $HCPtieqp

region mybox block 0 \${SIZE} 0 \${YIZE} 0 \${ZIZE} units box
create_box 2 mybox
create_atoms ${idx["Ti"]} box

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

thermo 1
thermo_style custom step temp pe ke etotal enthalpy
run 0

variable peat equal pe/atoms 
variable heat equal enthalpy/atoms

thermo_style custom step temp etotal enthalpy v_peat v_heat

fix REL all box/relax aniso $PRESS
#fix REL all box/relax iso $PRESS
min_style cg
minimize 1e-6 0.0 1000 100000

run 0
print "\${peat} \${heat}"
++

$lammps < $dir/tests/thermo/HCP_formation.lin > $dir/tests/thermo/HCP_formation_Ti.out
line=`tail -1 $dir/tests/thermo/HCP_formation_Ti.out`
ETi=`echo $line | awk '{print $1}'`
HTi=`echo $line | awk '{print $2}'`

cat > $dir/tests/thermo/HCP_formation.lin << ++
# computes formation energy for solid solution HCP ti-nb compound
units metal
atom_style atomic
boundary p p p

variable SIZE equal $HCPtieqp
variable YIZE equal sqrt(3)*$HCPtieqp
variable ZIZE equal sqrt(8/3)*$HCPtieqp

lattice hcp $HCPtieqp

region mybox block 0 \${SIZE} 0 \${YIZE} 0 \${ZIZE} units box
create_box 2 mybox
create_atoms ${idx["Nb"]} box

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

thermo 1
thermo_style custom step temp pe ke etotal enthalpy
run 0

variable peat equal pe/atoms 
variable heat equal enthalpy/atoms

thermo_style custom step temp etotal enthalpy v_peat v_heat

fix REL all box/relax aniso $PRESS
#fix REL all box/relax iso $PRESS
min_style cg
minimize 1e-6 0.0 1000 10000

run 0

print "\${peat} \${heat}"
++

$lammps < $dir/tests/thermo/HCP_formation.lin > $dir/tests/thermo/HCP_formation_Nb.out
line=`tail -1 $dir/tests/thermo/HCP_formation_Nb.out`
ENb=`echo $line | awk '{print $1}'`
HNb=`echo $line | awk '{print $2}'`

SIZE=7
NTIMES=10
CONC_PCTS=$(seq 5 5 95)

echo "# Concentration, formation energy" > $dir/tests/thermo/HCP_formation.dat
echo "0.0 0.00 0.00 0.00 0.00" >> $dir/tests/thermo/HCP_formation.dat
for CONC in $CONC_PCTS; do
PCT=`echo "$CONC/100" | bc -l`

rm -f tmp.dat; touch tmp.dat
for i in $(seq 1 $NTIMES); do

cat > $dir/tests/thermo/HCP_formation.lin << ++
# computes formation energy for solid solution HCP ti-nb compound
units metal
atom_style atomic
boundary p p p

variable SIZE equal $SIZE*$HCPtieqp
variable ETi  equal $ETi
variable ENb  equal $ENb
variable HTi  equal $HTi
variable HNb  equal $HNb
variable YIZE equal sqrt(3)*\${SIZE}
variable ZIZE equal sqrt(8/3)*\${SIZE}

lattice hcp $HCPtieqp

region mybox block 0 \${SIZE} 0 \${YIZE} 0 \${ZIZE} units box
create_box 2 mybox
create_atoms ${idx["Ti"]} box

set region mybox type/fraction 2 $PCT $RANDOM

group TI type 1
group NB type 2

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

thermo 1
thermo_style custom step temp pe ke etotal enthalpy press
run 0

variable Eform equal 1000*(pe-count(TI)*\${ETi}-count(NB)*\${ENb})/count(all)
variable Hform equal 1000*(pe-count(TI)*\${HTi}-count(NB)*\${HNb})/count(all)
# Gform uses solid solution entropy

thermo_style custom step temp pe ke etotal enthalpy press v_Eform v_Hform

fix REL all box/relax aniso $PRESS
#fix REL all box/relax iso $PRESS
min_style cg
minimize 1e-6 0.0 1000 10000

run 0
print "\${Eform} \${Hform}"
++
$lammps < $dir/tests/thermo/HCP_formation.lin > $dir/tests/thermo/HCP_formation.out
tail -1 $dir/tests/thermo/HCP_formation.out >> tmp.dat

done

Eavg=`awk '{sum += $1; sumsq += $1*$1} END {if (NR>0) print sum/NR, sqrt((sumsq-sum*sum/NR)/NR)}' tmp.dat`
Havg=`awk '{sum += $2; sumsq += $2*$2} END {if (NR>0) print sum/NR, sqrt((sumsq-sum*sum/NR)/NR)}' tmp.dat`
rm tmp.dat

echo "$PCT $Eavg $Havg" >> $dir/tests/thermo/HCP_formation.dat 
done
echo "1.00 0.00 0.00 0.00 0.00" >> $dir/tests/thermo/HCP_formation.dat
fi

########################################################################################
# 	TRANSFORMATION PATHWAYS
########################################################################################

if [ "$transition_flag" == "TRUE" ]; then

echo -e "${grn}COMPUTING PHASE TRANSITION PATHWAYS... ${non}"
#<<<<<<<<<<<<<<<<<<<<<<< alpha - alpha' - alpha'' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
echo "Martensite deformation..."
index=0
PDAMP=250.0
for CONC in 0.0 0.1 0.2 0.3 0.4 0.5; do
#for CONC in 0.0 0.125 0.25 0.375 0.40 0.425 0.45 0.475 0.48 0.485 0.49 0.495 0.50; do
index=$(($index+1))
touch $dir/tests/transitions/alpha_box.dat
cat > $dir/tests/transitions/alphas.lin << ++
units metal
boundary p p p
atom_style atomic

variable CONC   equal   $CONC
variable lat	equal	3.48
variable a1	equal	1.00
variable a2	equal	sqrt(3)
variable a3	equal	sqrt(8/3)
variable SIZE	equal	20
variable LX0	equal	\${SIZE}*\${lat}*\${a1}
variable LY0	equal	\${SIZE}*\${lat}*\${a2}
variable LZ0	equal	\${lat}*\${a3}
variable XY0	equal	0.0
variable XZ0	equal	0.0
variable YZ0	equal	0.0

variable CHIF   equal   0.5
variable NRUN	equal	100

lattice custom \${lat} &
	a1 \${a1}	0.0	0.0 &
	a2 0.0		\${a2}	0.0 &
	a3 0.0		0.0	\${a3} &
	basis 0.0	0.0	0.0 basis 0.5	0.5	0.0 &
	basis 0.0	0.33333 0.5 basis 0.5   0.83333 0.5

region box block 0 \${LX0} 0 \${LY0} 0 \${LZ0} units box
create_box 2 box

create_atoms ${idx["Nb"]} box

# set pair-style info
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

# define masses
mass 1 $mass1
mass 2 $mass2

timestep 0.001
thermo 100
thermo_style custom step temp pe ke etotal press lx ly lz xy xz yz zlo xlo xhi ylo yhi zlo zhi
run 0

variable tmp	equal pe/atoms
variable NbEAT  equal \${tmp}

set region box type/fraction 1 1.0 1

run 0 

variable tmp   equal pe/atoms
variable TiEAT equal \${tmp}

#region hbas1 sphere 0.00 0.33333 0.5 0.01 units lattice 
#region hbas2 sphere 0.50 0.83333 0.5 0.01 units lattice
#region hbas  union 2 hbas1 hbas2
#group hbas region hbas

#replicate \${SIZE} \${SIZE} 1
#replicate \${SIZE} \${SIZE} \${SIZE}

set region box type/fraction 2 \${CONC} $RANDOM

variable tmp equal lx
print \${tmp}
run 0
fix rel all press/berendsen aniso 0.0 0.0 $PDAMP
run 200
unfix rel
print \${tmp}

variable tmp equal pe/atoms
variable peat0 equal \${tmp}
group Ti type 1
group Nb type 2 
variable tmp equal pe-count(Ti)*\${TiEAT}-count(Nb)*\${NbEAT}
variable eform0 equal \${tmp} 

run 0

variable tmp equal zlo+lz*0.40
variable zmin equal \${tmp}
variable tmp equal zlo+lz*0.60
variable zmax equal \${tmp}
variable tmp  equal xhi
variable XHI   equal \${tmp}
variable tmp  equal xlo
variable XLO   equal \${tmp}
variable tmp  equal yhi
variable YHI  equal \${tmp}
variable tmp  equal ylo
variable YLO  equal \${tmp}

#region hbas block \${XLO} \${XHI} \${YLO} \${YHI} \${zmin} \${zmax} units box
region hbas plane 0.0 0.0 \${zmin} 0.0 0.0 \${zmin} side out units box
group hbas region hbas

variable Vx equal 0.0
variable Vy equal 100*ly*\${CHIF}/(\${SIZE}*\${NRUN})
variable Vz equal 0.0
variable Gamma equal 180-2*180*atan(ly/lx)/3.141592654
variable coa   equal lz/(lx/\${SIZE})
#variable dele  equal 1000*(pe/atoms-\${peat0})
#variable eform equal 1000*(pe-\${TiEAT}*count(Ti)-\${NbEAT}*count(Nb))/atoms
variable dele equal 1000*(pe-\${TiEAT}*count(Ti)-\${NbEAT}*count(Nb)-\${eform0})/atoms

reset_timestep 0

variable JRUN equal 100
variable DY    equal (0.10*ly/\${SIZE})/(\${NRUN})
variable chid  equal step*\${DY}/\${JRUN}
variable chi   equal v_chid/(lx/\${SIZE}) 

fix P all print \${JRUN} "\${chid} \${chi} \${dele} \${Gamma} \${coa}" file $dir/tests/transitions/chi_e_$CONC.dat screen no
#dump D all custom \${JRUN} chi_$CONC.atom id element xs ys zs
#dump_modify D element Ti Nb

label		jLOOP
 variable	j loop \${NRUN}
  
 displace_atoms hbas move 0.0 \${DY} 0.0 units box
 fix rel all press/berendsen aniso 0.0 0.0 $PDAMP 
 run \${JRUN}
 unfix rel

 next j
 if "\$j >= \${NRUN}" then "jump $dir/tests/transitions/alphas.lin jBREAK"
 jump $dir/tests/transitions/alphas.lin jLOOP

label		jBREAK
variable 	j delete

print "\${CONC} \${Gamma} \${coa}"
++

$lammps -in $dir/tests/transitions/alphas.lin > $dir/tests/transitions/alphas_log

#tail -1 log$CONC >> box.dat
coa=`cat $dir/tests/transitions/chi_e_$CONC.dat | sort -gk3 | awk 'NR==2{print $4}'`
gam=`cat $dir/tests/transitions/chi_e_$CONC.dat | sort -gk3 | awk 'NR==2{print $5}'`
echo "$CONC $gam $coa" >> $dir/tests/transitions/alpha_box.dat
done

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BETA-OMEGA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo "Beta to omega..."
if [ "$betaOmegaLammpsFlag" == "TRUE" ]; then
numpts=100

for PRESS in 0 25 50 75 100; do

cat >$dir/tests/transitions/beta_omega.lin << --
################################################
#  CALCULATES TRANSITION PATHWAY FOR BETA-L60
# TO OMEGA-LAZAR TI3NB
################################################

#define simulation region and bcc grid
units		metal
atom_style	atomic

#initialization variables
variable	NMESH equal $(($numpts+1))
variable	aOMGA equal ${omgPLAT[$PRESS]}
variable	aBETA equal ${G1PLAT[$PRESS]}

# a1 lattice constnats
variable 	aOMGAs equal \${aOMGA}
variable	aBETAs equal 2*sqrt(2)*\${aBETA}

# a2 lattice constants
variable	bOMGAs equal \${aOMGA}*sqrt(3)/2
variable	bBETAs equal sqrt(6)*\${aBETA} # hexagonalar 

# a3 lattice constants
variable	cOMGAs equal 0.3*\${aOMGA}
variable	cBETAs equal sqrt(3/8)*\${aBETAs}/2

# tilt factors
variable xyOMGAs equal \${aOMGAs}/2
variable xyBETAs equal \${aBETAs}/2

variable	sqrt3  equal sqrt(3)
variable	sqrt32  equal sqrt(3)/2
variable	sqrt6  equal sqrt(6)
variable	sqrt2  equal sqrt(2)
variable	sqrt38  equal sqrt(3/8)

lattice		custom \${aOMGA} &
		a1 1.0 0.0 0.0 &
		a2 0.5 \${sqrt32} 0.0 &
		a3 0.0 0.0 0.3 &
		basis 0.5 0.5 0.0 basis 0.0 0.5 0.0 basis 0.5 0.0 0.0 basis 0.8333333 0.3333333 0.5 &
		basis 0.3333333 0.8333333 0.5 basis 0.8333333 0.8333333 0.5 &
		basis 0.1666667 0.1666667 0.5 basis 0.6666667 0.1666667 0.5 &
		basis 0.1666667 0.6666667 0.5 basis 0.0 0.0 0.0 &
		basis 0.3333333 0.3333333 0.5 basis 0.6666667 0.6666667 0.5

# now create box and define atom basis
box tilt large
region		mybox prism 0 \${aOMGAs} 0 \${bOMGAs} 0 \${cOMGAs} \${xyOMGAs} 0.0 0.0 units box
create_box	2 mybox

create_atoms ${idx["Nb"]} box &
		basis 1 ${idx["Ti"]} basis 2 ${idx["Ti"]} basis 3 ${idx["Ti"]} basis 4 ${idx["Ti"]} basis 5 ${idx["Ti"]} basis 6 ${idx["Ti"]} &
		basis 7 ${idx["Ti"]} basis 8 ${idx["Ti"]} basis 9 ${idx["Ti"]} basis 10 ${idx["Nb"]} basis 11 ${idx["Nb"]} basis 12 ${idx["Nb"]}

group		whole region mybox

mass		1 $mass1 
mass		2 $mass2

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

variable dn1x equal 0.1666667*\${aOMGAs}+0.1666667*0.5*\${aOMGAs}
variable dn1y equal 0.1666667*(\${sqrt32}*\${aOMGAs})
variable dn1z equal 0.5*\${cOMGAs}
variable dn2x equal 0.1666667*\${aOMGAs}+0.666667*0.5*\${aOMGAs}
variable dn2y equal 0.6666667*(\${sqrt32}*\${aOMGAs})
variable dn2z equal 0.5*\${cOMGAs}
variable dn3x equal 0.6666667*\${aOMGAs}+0.1666667*0.5*\${aOMGAs}
variable dn3y equal 0.1666667*(\${sqrt32}*\${aOMGAs})
variable dn3z equal 0.5*\${cOMGAs}
variable dn4x equal 0.6666667*\${aOMGAs}+0.666667*0.5*\${aOMGAs}
variable dn4y equal 0.6666667*(\${sqrt32}*\${aOMGAs})
variable dn4z equal 0.5*\${cOMGAs}

region dn1 sphere \${dn1x} \${dn1y} \${dn1z} 0.1 units box
region dn2 sphere \${dn2x} \${dn2y} \${dn2z} 0.1 units box
region dn3 sphere \${dn3x} \${dn3y} \${dn3z} 0.1 units box
region dn4 sphere \${dn4x} \${dn4y} \${dn4z} 0.1 units box
group dn1 region dn1
group dn2 region dn2
group dn3 region dn3
group dn4 region dn4
variable up1x equal 0.3333333*\${aOMGAs}+0.3333333*0.5*\${aOMGAs}
variable up1y equal 0.3333333*(\${sqrt32}*\${aOMGAs})
variable up1z equal 0.5*\${cOMGAs}
variable up2x equal 0.3333333*\${aOMGAs}+0.8333333*0.5*\${aOMGAs}
variable up2y equal 0.8333333*(\${sqrt32}*\${aOMGAs})
variable up2z equal 0.5*\${cOMGAs}
variable up3x equal 0.8333333*\${aOMGAs}+0.3333333*0.5*\${aOMGAs}
variable up3y equal 0.3333333*(\${sqrt32}*\${aOMGAs})
variable up3z equal 0.5*\${cOMGAs}
variable up4x equal 0.8333333*\${aOMGAs}+0.8333333*0.5*\${aOMGAs}
variable up4y equal 0.8333333*(\${sqrt32}*\${aOMGAs})
variable up4z equal 0.5*\${cOMGAs}

region up1 sphere \${up1x} \${up1y} \${up1z} 0.1 units box
region up2 sphere \${up2x} \${up2y} \${up2z} 0.1 units box
region up3 sphere \${up3x} \${up3y} \${up3z} 0.1 units box
region up4 sphere \${up4x} \${up4y} \${up4z} 0.1 units box
group up1 region dn1
group up2 region dn2
group up3 region dn3
group up4 region dn4


region up union 4 up1 up2 up3 up4
region dn union 4 dn1 dn2 dn3 dn4

group up region up
group dn region dn

# define simulation variables
timestep 0.001
thermo 1
thermo_style custom step pe etotal enthalpy press vol fmax xlat ylat zlat 
dump D all custom 1 dump.atom type xs ys zs
run 0
# initialize output
variable 	eat	equal pe/atoms
variable	hat	equal enthalpy/atoms
variable	delj	equal 0
variable 	eat0    equal \${eat}
variable 	hat0    equal \${hat}
variable 	dele    equal 1000*(v_eat-\${eat0})
variable 	delh    equal 1000*(v_hat-\${hat0})

# output fix every two(!)
fix P all print 1 "\${delj} \${dele} \${delh}" screen no file $dir/tests/transitions/beta_omega_$PRESS.dat
#dump D all custom 1 dump_$PRESS.atom type xs ys zs
 
	# and reset variables for after j loop

 label		jLOOP
  variable	j loop \${NMESH}
  variable	delj equal (v_j-1)/(\${NMESH}-1)					    
  
  variable 	XHI	 equal ((1-v_delj)*\${aOMGAs}+v_delj*\${aBETAs})  # switching function for shear
  variable 	YHI	 equal ((1-v_delj)*\${bOMGAs}+v_delj*\${bBETAs}) 	    # switching function for shear
  variable 	ZHI	 equal ((1-v_delj)*\${cOMGAs}+v_delj*\${cBETAs})  	    # switching function for shear
  variable 	XYHI	 equal ((1-v_delj)*\${xyOMGAs}+v_delj*\${xyBETAs})  	    # switching function for shear

 
 # define shift variables 
 variable 	DZZ	equal 	-\${ZHI}/6/\${NMESH}				 	    
 variable 	UZZ	equal 	-\${DZZ}
 variable	RDZ	equal	\${ZHI}/6
 variable	RUZ	equal	-\${RDZ}
 
  change_box all x final 0 \${XHI} y final 0 \${YHI} z final 0 \${ZHI} xy final \${XYHI} remap units box
 
  displace_atoms up move 0.0 0.0 \${UZZ} units box
  displace_atoms dn move 0.0 0.0 \${DZZ} units box

  run 1
 
  next j
  jump $dir/tests/transitions/beta_omega.lin jLOOP

print "Success!"
--

$lammps < $dir/tests/transitions/beta_omega.lin > $dir/tests/transitions/beta_omega_$PRESS.out
done

else

numpts=100
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/transitions/beta_omega.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

numpts=100

for PRESS in 0 10 20 30; do

rm -f $dir/tests/transitions/beta_omega.conf; touch $dir/tests/transitions/beta_omega.conf
for iii in $(seq 0 $numpts); do
	deli=`echo "$iii/$numpts" | bc -l`
	python $betaOmegaAlloy ${G1PLAT["$PRESS"]} ${omgPLAT["$PRESS"]} $deli $deli conf >> $dir/tests/transitions/beta_omega.conf 
done

# compute with meamzilla
$meamz -p $dir/tests/meamz_params > $dir/tests/transitions/meamz_betaomega.out

# compute energy
echo "#del, E" > $dir/tests/transitions/beta_omega$PRESS.dat
e0=`awk 'NR==3{print $6}' data.energy`
awk -v e0=$e0 -v np=$numpts 'NR>1{
					   getline;
					   print $1/np, 1000*($6-e0);
				     }' data.energy >> $dir/tests/transitions/beta_omega$PRESS.dat
done
fi

## PURE TITANIUM BETA TO OMEGA ##
##  beta o omega pathway ##

for PRESS in 0 25 50 75 100; do


cat > $dir/tests/transitions/betaOmega.lin << --
###########################################
#       LAMMPS script for computing
#       melting temperature
###########################################

units metal
boundary p p p

lattice custom ${OMGtiPLAT["$PRESS"]} a1 1.0000 0.0 0.0 a2 0.00 1.732051 0.0 a3 0.0 0.0 0.612367 &
                    basis 0.0 0.0 0.0 basis 0.5 0.5 0.0 &
                    basis 0.0 0.3333 0.3333 basis 0.5 0.8333 0.3333 &
                    basis 0.0 0.6667 0.6667 basis 0.5 0.1667 0.6667

region whole block 0 1 0 1 0 1
create_box 2 whole
create_atoms ${idx["Ti"]} box
mass 1 $mass1
mass 2 $mass2
variable jmax equal 100

#define potential
pair_style	$pairstyle
pair_coeff	* * $potPath $elem1 $elem2

region mid11 sphere 0.0000 0.3333 0.3333 0.01 units lattice
region mid12 sphere 0.5000 0.8333 0.3333 0.01 units lattice
region mid21 sphere 0.0000 0.6667 0.6667 0.01 units lattice
region mid22 sphere 0.5000 0.1667 0.6667 0.01 units lattice

region mid1 union 2 mid11 mid12
region mid2 union 2 mid21 mid22

group mid1 region mid1
group mid2 region mid2

thermo 1
variable peat equal pe/atoms
variable delt equal v_j/\${jmax}

fix out all print 1 "\${delt} \${peat}" file $dir/tests/transitions/betaOmegaTi$PRESS.dat

variable xd equal (0.16667)/\${jmax}
variable j loop 0 \${jmax}

label loop

displace_atoms mid1 move 0.0 0.0 \${xd} units lattice
displace_atoms mid2 move 0.0 0.0 -\${xd} units lattice

run 1

next j
jump $dir/tests/transitions/betaOmega.lin loop
--

$lammps < $dir/tests/transitions/betaOmega.lin > $dir/tests/transitions/betaOmegaTi.out 

eW=`tail -1 "$dir/tests/transitions/betaOmegaTi$PRESS.dat" | awk '{print $2}'`
sed 1d $dir/tests/transitions/betaOmegaTi$PRESS.dat | awk -v ew=$eW '{print $1, $2, 1000*($2-ew)}' > tmp
mv -f tmp $dir/tests/transitions/betaOmegaTi$PRESS.dat
#sed -i '1 i\ ## react coord xhi, e(xhi), e(xhi)-e_omega' $dir/tests/transitions/betaOmega$PRESS.dat

done
## end beta to omega pathway ##

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BETA-alpha'' >>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo "beta to alpha''"
cat > $dir/tests/meamz_params <<@@
ngroups 1
optstyle powell

num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot end
tempfile temp
config $dir/tests/transitions/beta_app.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@

numpts=40
python $betaApp $G1eqp $APPeqp $numpts $numpts conf > $dir/tests/transitions/beta_app.conf

$meamz -p $dir/tests/meamz_params > $dir/tests/transitions/meamz_betaapp.out

# compute energy
rm -f $dir/tests/transitions/beta_app.dat
echo "#del, E" > $dir/tests/transitions/beta_app.dat
e0=`sort -gk6 data.energy | head -3 | tail -1 | awk '{print $6}'`
awk -v e0=$e0 -v np=$numpts 'NR>2{
					   i = int(($1)/(np+1))
					   j = ($1)%(np+1)
					   print i/np, j/np, 1000*($6-e0), $6;
				     }' data.energy >> $dir/tests/transitions/beta_app.dat 

# compute transition enthalpy
rm -f $dir/tests/transitions/beta_app_enth.dat
echo "#del, enth" > $dir/tests/transitions/beta_app_enth.dat
awk 'NR>2{press=$5; getline; press+=$5; getline; press+=$5; press/=3; print $1, press;
	 getline; getline; getline}' data.stress > tmp1
grep '##' $dir/tests/transitions/beta_app.conf | awk '{print $4}' > tmp2
echo "##" > ptmp
paste tmp1 tmp2 >> ptmp
#rm tmp1 tmp2
paste $dir/tests/transitions/beta_app.dat ptmp | awk 'NR>1{print $1, $2, 1000*($4+$6*$7)}' > tmp3
h0=`awk 'NR==1{print $3}' tmp3`
awk -v h0=$h0 '{print $1, $2, $3-h0}' tmp3 >> $dir/tests/transitions/beta_app_enth.dat

############ alpha to beta in pure Ti ###################

#generate config file
NMESH=40
python /n/jww-1/ehemann.2/testingScripts/burgers.py $BCCtieqp $HCPtieqp $HCPticoap $NMESH $NMESH conf > $dir/tests/transitions/alphabetaTi.conf
cat > $dir/tests/transitions/alphabetaTi_params <<@@
ngroups 1

optstyle powell
num_powell 0
init_scale 10.0
pop_size 1
cross_rate 0.0
mut_rate 0.0
fit_rate 0.0
rescale_rate 0.0
order_breed 1
gen_save 1

rescale 0
embed_extrap 0

startpot $mmzpot
endpot junk
tempfile temp
config $dir/tests/transitions/alphabetaTi.conf
lammpsfile lmp.pt

energy_weight 10.0
stress_weight 10.0

d_eps 0.0
max_steps 0

seed 1
@@
$meamz -p $dir/tests/transitions/alphabetaTi_params > $dir/tests/transitions/alphabetaTi.out
rm -f $dir/tests/transitions/alphaBetaTi.dat; echo "## lambda_1, lambda_2, e/atom" > $dir/tests/transitions/alphaBetaTi.dat

ehcpburg=`sort -gk6 data.energy | head -3 | tail -1 | awk '{print $6}'`
for i in $(seq 0 $NMESH); do
        for j in $(seq 0 $NMESH); do
                linnum=$((3 + $i*$(($NMESH + 1)) + $j))

                lami=`echo "$i/$NMESH" | bc -l`
                lamj=`echo "$j/$NMESH" | bc -l`

                en=`awk -v nr=$linnum -v ea=$ehcpburg 'NR==nr{print 1000*($6-ea)}' data.energy`

                echo "$lami $lamj $en" >> $dir/tests/transitions/alphaBetaTi.dat
        done
done

fi

echo "FINISHED TESTS!"
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#	P L O T T I N G     P O R T I O N
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
########################################################################################
#	SPLINE FUNCTIONS
########################################################################################
echo "plotting splines..."

if [ "$typ" == "MEAM" ]; then

	N_PAIR_AA=`awk 'NR==4{print $1}' $dir/lammps.pt`
	N_PAIR_AB=`awk -v nr=$(($N_PAIR_AA + 7)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_PAIR_BB=`awk -v nr=$(($N_PAIR_AA + $N_PAIR_AB + 10)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_PAIR=$(($N_PAIR_AA + $N_PAIR_AB + $N_PAIR_BB))

	N_RHO_A=`awk -v nr=$(($NTOT_PAIR + 13)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_RHO_B=`awk -v nr=$(($NTOT_PAIR + $N_RHO_A + 16)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_RHO=$(($NTOT_PAIR + $N_RHO_A + $N_RHO_B))

	N_EMBED_A=`awk -v nr=$(($NTOT_RHO + 19)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_EMBED_B=`awk -v nr=$(($NTOT_RHO + $N_EMBED_A + 22)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_EMB=$(($NTOT_RHO + $N_EMBED_A + $N_EMBED_B))

	N_F_AA=`awk -v nr=$(($NTOT_EMB + 25)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_F_AB=`awk -v nr=$(($NTOT_EMB + $N_F_AA + 28)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_F_BB=`awk -v nr=$(($NTOT_EMB + $N_F_AA + $N_F_AB + 31)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_F=$(($NTOT_EMB + $N_F_AA + $N_F_AB + $N_F_BB))
	
	N_G_A=`awk -v nr=$(($NTOT_F + 34)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_G_B=`awk -v nr=$(($NTOT_F + $N_G_A + 37)) 'NR==nr{print $1}' $dir/lammps.pt`

	LB_PAIR_AA=6;
	LE_PAIR_AA=$(($LB_PAIR_AA + $N_PAIR_AA - 1));
	sed -n "$(($LB_PAIR_AA-2))","$LE_PAIR_AA"p $dir/lammps.pt > phiaatmp
 	python $cspline phiaatmp $dir/tests/phiaaspline; rm -f phiaatmp
	echo "phi_aa..."

	LB_PAIR_AB=$(($LE_PAIR_AA + 4));
	LE_PAIR_AB=$(($LB_PAIR_AB + $N_PAIR_AB - 1));
	sed -n "$(($LB_PAIR_AB-2))","$LE_PAIR_AB"p $dir/lammps.pt > phiabtmp
 	python $cspline phiabtmp $dir/tests/phiabspline; rm -f phiabtmp
	echo "phi_ab..."

	LB_PAIR_BB=$(($LE_PAIR_AB + 4));
	LE_PAIR_BB=$(($LB_PAIR_BB + $N_PAIR_BB - 1));
	sed -n "$(($LB_PAIR_BB-2))","$LE_PAIR_BB"p $dir/lammps.pt > phibbtmp
 	python $cspline phibbtmp $dir/tests/phibbspline; rm -f phibbtmp
	echo "phi_bb..."

	LB_RHO_A=$(($LE_PAIR_BB + 4));
	LE_RHO_A=$(($LB_RHO_A + $N_RHO_A - 1));
	sed -n "$(($LB_RHO_A-2))","$LE_RHO_A"p $dir/lammps.pt > rhoatmp
	python $cspline rhoatmp $dir/tests/rhoaspline; rm -f rhoatmp
	echo "rho_a..."
	
	LB_RHO_B=$(($LE_RHO_A + 4));
	LE_RHO_B=$(($LB_RHO_B + $N_RHO_B - 1));
	sed -n "$(($LB_RHO_B-2))","$LE_RHO_B"p $dir/lammps.pt > rhobtmp
	python $cspline rhobtmp $dir/tests/rhobspline; rm -f rhobtmp
	echo "rho_b..."
	echo $LB_RHO_B $LE_RHO_B

	LB_EMBED_A=$(($LE_RHO_B + 4));
	LE_EMBED_A=$(($LB_EMBED_A + $N_EMBED_A - 1));
	sed -n "$(($LB_EMBED_A-2))","$LE_EMBED_A"p $dir/lammps.pt > Uatmp
	python $cspline Uatmp $dir/tests/Uaspline; rm -f Uatmp
	echo "U_a..."

	LB_EMBED_B=$(($LE_EMBED_A + 4));
	LE_EMBED_B=$(($LB_EMBED_B + $N_EMBED_B - 1));
	sed -n "$(($LB_EMBED_B-2))","$LE_EMBED_B"p $dir/lammps.pt > Ubtmp
	python $cspline Ubtmp $dir/tests/Ubspline; rm -f Ubtmp
	echo "U_b..."

	LB_F_AA=$(($LE_EMBED_B + 4));
	LE_F_AA=$(($LB_F_AA + $N_F_AA - 1));
	sed -n "$(($LB_F_AA-2))","$LE_F_AA"p $dir/lammps.pt > faatmp
	python $cspline faatmp $dir/tests/faaspline; rm -f faatmp
	echo "f_aa..."

	LB_F_AB=$(($LE_F_AA + 4));
	LE_F_AB=$(($LB_F_AB + $N_F_AB - 1));
	sed -n "$(($LB_F_AB-2))","$LE_F_AB"p $dir/lammps.pt > fabtmp
	python $cspline fabtmp $dir/tests/fabspline; rm -f fabtmp
	echo "f_ab..."

	LB_F_BB=$(($LE_F_AB + 4));
	LE_F_BB=$(($LB_F_BB + $N_F_BB - 1));
	sed -n "$(($LB_F_BB-2))","$LE_F_BB"p $dir/lammps.pt > fbbtmp
	python $cspline fbbtmp $dir/tests/fbbspline; rm -f fbbtmp
	echo "f_bb..."

	LB_G_A=$(($LE_F_BB + 4));
	LE_G_A=$(($LB_G_A + $N_G_A - 1));
	sed -n "$(($LB_G_A-2))","$LE_G_A"p $dir/lammps.pt > gatmp
	python $cspline gatmp $dir/tests/gaspline; rm -f gatmp
	echo "g_a..."

	LB_G_B=$(($LE_G_A + 4));
	LE_G_B=$(($LB_G_B + $N_G_B - 1));
	sed -n "$(($LB_G_B-2))","$LE_G_B"p $dir/lammps.pt > gbtmp
	python $cspline gbtmp $dir/tests/gbspline; rm -f gbtmp
	echo "g_a..."

	echo "set terminal postscript enhanced color font ',5'
	set encoding iso_8859_1
	set output '$dir/tests/SPLINES.ps'
	set size 8,6
	set multiplot layout 4,3
	
	##set title '$elem1-$elem1 Pair Potential'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol f}(r) [eV]'
	plot '$dir/lammps.pt' every ::$(($LB_PAIR_AA-2))::$LE_PAIR_AA u 1:2 w p pt 7 lc 1 notitle,\
	'$dir/tests/phiaaspline' u 1:2 w l lt 1 lc 1 notitle
	
	##set title '$elem1-$elem2 Pair Potential'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol f}(r) [eV]'
	plot '$dir/lammps.pt' every ::$(($LB_PAIR_AB-2))::$LE_PAIR_AB u 1:2 w p pt 7 lc rgb 'purple' notitle,\
	'$dir/tests/phiabspline' u 1:2 w l lt 1 lc rgb 'purple' notitle
	
	##set title '$elem2-$elem2 Pair Potential'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol f}(r) [eV]'
	plot '$dir/lammps.pt' every ::$(($LB_PAIR_BB-2))::$LE_PAIR_BB u 1:2 w p pt 7 lc rgb 'blue' notitle,\
	'$dir/tests/phibbspline' u 1:2 w l lt 1 lc rgb 'blue' notitle
	
	##set title '$elem1-$elem1 SW Weight Function'
	set xlabel 'r [{\305}]'
	set ylabel 'f(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_F_AA-2))::$LE_F_AA u 1:2 w p pt 7 lc 1 notitle,\
	'$dir/tests/faaspline' u 1:2 w l lt 1 lc 1 notitle
	
	##set title '$elem1-$elem2 SW Weight Function'
	set xlabel 'r [{\305}]'
	set ylabel 'f(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_F_AB-2))::$LE_F_AB u 1:2 w p pt 7 lc rgb 'purple' notitle,\
	'$dir/tests/fabspline' u 1:2 w l lt 1 lc rgb 'purple' notitle
	
	##set title '$elem2-$elem2 SW Weight Function'
	set xlabel 'r [{\305}]'
	set ylabel 'f(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_F_BB-2))::$LE_F_BB u 1:2 w p pt 7 lc rgb 'blue' notitle,\
	'$dir/tests/fbbspline' u 1:2 w l lt 1 lc rgb 'blue' notitle
	
	##set title '$elem1 Density Function'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol r}(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_RHO_A-2))::$LE_RHO_A u 1:2 w p pt 7 lc 1 notitle,\
	'$dir/tests/rhoaspline' u 1:2 w l lt 1 lc 1 notitle
	
	##set title '$elem1 Embedding Function'
	set xlabel '{/Symbol r}'
	set ylabel 'U({/Symbol r}) [eV]'
	##set xrange [0:1]
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_EMBED_A-2))::$LE_EMBED_A u 1:2 w p pt 7 lc 1 notitle,\
	'$dir/tests/Uaspline' u 1:2 w l lt 1 lc 1 notitle

	##set title '$elem1 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_A-2))::$LE_G_A u 1:2 w p pt 7 lc 1 notitle,\
	'$dir/tests/gaspline' u 1:2 w l lt 1 lc 1 notitle

	set xrange [*:*]
	set title '$elem2 Density Function'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol r}(r)'
	plot '$dir/lammps.pt' every ::$(($LB_RHO_B-2))::$LE_RHO_B u 1:2 w p pt 7 lc rgb 'blue' notitle,\
	'$dir/tests/rhobspline' u 1:2 w l lt 1 lc rgb 'blue' notitle
	
	##set title '$elem2 Embedding Function'
	set xlabel '{/Symbol r}'
	set ylabel 'U({/Symbol r}) [eV]'
	##set xrange [0:1]
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_EMBED_B-2))::$LE_EMBED_B u 1:2 w p pt 7 lc rgb 'blue' notitle,\
	'$dir/tests/Ubspline' u 1:2 w l lt 1 lc rgb 'blue' notitle

	##set title '$elem2 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_B-2))::$LE_G_B u 1:2 w p pt 7 lc rgb 'blue' notitle,\
	'$dir/tests/gbspline' u 1:2 w l lt 1 lc rgb 'blue' notitle
	
	unset multiplot" | gnuplot

#rm -f $dir/tests/*spline ## clean up

$ps2png $dir/tests/SPLINES.ps $dir/tests/SPLINES.png

elif [ "$typ" == "GMEAM" ]; then

	N_PAIR_AA=`awk 'NR==4{print $1}' $dir/lammps.pt`
	N_PAIR_AB=`awk -v nr=$(($N_PAIR_AA + 7)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_PAIR_BB=`awk -v nr=$(($N_PAIR_AA + $N_PAIR_AB + 10)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_PAIR=$(($N_PAIR_AA + $N_PAIR_AB + $N_PAIR_BB))

	N_RHO_A=`awk -v nr=$(($NTOT_PAIR + 13)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_RHO_B=`awk -v nr=$(($NTOT_PAIR + $N_RHO_A + 16)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_RHO=$(($NTOT_PAIR + $N_RHO_A + $N_RHO_B))

	N_EMBED_A=`awk -v nr=$(($NTOT_RHO + 19)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_EMBED_B=`awk -v nr=$(($NTOT_RHO + $N_EMBED_A + 22)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_EMB=$(($NTOT_RHO + $N_EMBED_A + $N_EMBED_B))

	N_F_AA=`awk -v nr=$(($NTOT_EMB + 25)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_F_AB=`awk -v nr=$(($NTOT_EMB + $N_F_AA + 28)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_F_BB=`awk -v nr=$(($NTOT_EMB + $N_F_AA + $N_F_AB + 31)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_F=$(($NTOT_EMB + $N_F_AA + $N_F_AB + $N_F_BB))
	
	N_G_AAA=`awk -v nr=$(($NTOT_F + 34)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_G_AAB=`awk -v nr=$(($NTOT_F + $N_G_AAA + 37)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_G_ABB=`awk -v nr=$(($NTOT_F + $N_G_AAA + $N_G_AAB + 40)) 'NR==nr{print $1}' $dir/lammps.pt`

	NTOT_GA=$(($NTOT_F + $N_G_AAA + $N_G_AAB + $N_G_ABB))

	N_G_BAA=`awk -v nr=$(($NTOT_GA + 43)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_G_BAB=`awk -v nr=$(($NTOT_GA + $N_G_BAA + 46)) 'NR==nr{print $1}' $dir/lammps.pt`
	N_G_BBB=`awk -v nr=$(($NTOT_GA + $N_G_BAA + $N_G_BAB + 49)) 'NR==nr{print $1}' $dir/lammps.pt`

	LB_PAIR_AA=6;
	LE_PAIR_AA=$(($LB_PAIR_AA + $N_PAIR_AA - 1));
	sed -n "$(($LB_PAIR_AA-2))","$LE_PAIR_AA"p $dir/lammps.pt > phiaatmp
 	python $cspline phiaatmp $dir/tests/phiaaspline; rm -f phiaatmp
	echo "phi_aa..."

	LB_PAIR_AB=$(($LE_PAIR_AA + 4));
	LE_PAIR_AB=$(($LB_PAIR_AB + $N_PAIR_AB - 1));
	sed -n "$(($LB_PAIR_AB-2))","$LE_PAIR_AB"p $dir/lammps.pt > phiabtmp
 	python $cspline phiabtmp $dir/tests/phiabspline; rm -f phiabtmp
	echo "phi_ab..."

	LB_PAIR_BB=$(($LE_PAIR_AB + 4));
	LE_PAIR_BB=$(($LB_PAIR_BB + $N_PAIR_BB - 1));
	sed -n "$(($LB_PAIR_BB-2))","$LE_PAIR_BB"p $dir/lammps.pt > phibbtmp
 	python $cspline phibbtmp $dir/tests/phibbspline; rm -f phibbtmp
	echo "phi_bb..."

	LB_RHO_A=$(($LE_PAIR_BB + 4));
	LE_RHO_A=$(($LB_RHO_A + $N_RHO_A - 1));
	sed -n "$(($LB_RHO_A-2))","$LE_RHO_A"p $dir/lammps.pt > rhoatmp
	python $cspline rhoatmp $dir/tests/rhoaspline; rm -f rhoatmp
	echo "rho_a..."
	
	LB_RHO_B=$(($LE_RHO_A + 4));
	LE_RHO_B=$(($LB_RHO_B + $N_RHO_B - 1));
	sed -n "$(($LB_RHO_B-2))","$LE_RHO_B"p $dir/lammps.pt > rhobtmp
	python $cspline rhobtmp $dir/tests/rhobspline; rm -f rhobtmp
	echo "rho_b..."
	echo $LB_RHO_B $LE_RHO_B

	LB_EMBED_A=$(($LE_RHO_B + 4));
	LE_EMBED_A=$(($LB_EMBED_A + $N_EMBED_A - 1));
	sed -n "$(($LB_EMBED_A-2))","$LE_EMBED_A"p $dir/lammps.pt > Uatmp
	python $cspline Uatmp $dir/tests/Uaspline; rm -f Uatmp
	echo "U_a..."

	LB_EMBED_B=$(($LE_EMBED_A + 4));
	LE_EMBED_B=$(($LB_EMBED_B + $N_EMBED_B - 1));
	sed -n "$(($LB_EMBED_B-2))","$LE_EMBED_B"p $dir/lammps.pt > Ubtmp
	python $cspline Ubtmp $dir/tests/Ubspline; rm -f Ubtmp
	echo "U_b..."

	LB_F_AA=$(($LE_EMBED_B + 4));
	LE_F_AA=$(($LB_F_AA + $N_F_AA - 1));
	sed -n "$(($LB_F_AA-2))","$LE_F_AA"p $dir/lammps.pt > faatmp
	python $cspline faatmp $dir/tests/faaspline; rm -f faatmp
	echo "f_aa..."

	LB_F_AB=$(($LE_F_AA + 4));
	LE_F_AB=$(($LB_F_AB + $N_F_AB - 1));
	sed -n "$(($LB_F_AB-2))","$LE_F_AB"p $dir/lammps.pt > fabtmp
	python $cspline fabtmp $dir/tests/fabspline; rm -f fabtmp
	echo "f_ab..."

	LB_F_BB=$(($LE_F_AB + 4));
	LE_F_BB=$(($LB_F_BB + $N_F_BB - 1));
	sed -n "$(($LB_F_BB-2))","$LE_F_BB"p $dir/lammps.pt > fbbtmp
	python $cspline fbbtmp $dir/tests/fbbspline; rm -f fbbtmp
	echo "f_bb..."

	LB_G_AAA=$(($LE_F_BB + 4));
	LE_G_AAA=$(($LB_G_AAA + $N_G_AAA - 1));
	sed -n "$(($LB_G_AAA-2))","$LE_G_AAA"p $dir/lammps.pt > gatmp
	python $cspline gatmp $dir/tests/gaaaspline; rm -f gatmp
	echo "g_aaa..."
	
	LB_G_AAB=$(($LE_G_AAA + 4));
	LE_G_AAB=$(($LB_G_AAB + $N_G_AAB - 1));
	sed -n "$(($LB_G_AAB-2))","$LE_G_AAB"p $dir/lammps.pt > gatmp
	python $cspline gatmp $dir/tests/gaabspline; rm -f gatmp
	echo "g_aab..."

	LB_G_ABB=$(($LE_G_AAB + 4));
	LE_G_ABB=$(($LB_G_ABB + $N_G_ABB - 1));
	sed -n "$(($LB_G_ABB-2))","$LE_G_ABB"p $dir/lammps.pt > gatmp
	python $cspline gatmp $dir/tests/gabbspline; rm -f gatmp
	echo "g_abb..."

	LB_G_BAA=$(($LE_G_ABB + 4));
	LE_G_BAA=$(($LB_G_BAA + $N_G_BAA - 1));
	sed -n "$(($LB_G_BAA-2))","$LE_G_BAA"p $dir/lammps.pt > gatmp
	python $cspline gatmp $dir/tests/gbaaspline; rm -f gatmp
	echo "g_baa..."

	LB_G_BAB=$(($LE_G_BAA + 4));
	LE_G_BAB=$(($LB_G_BAB + $N_G_BAB - 1));
	sed -n "$(($LB_G_BAB-2))","$LE_G_BAB"p $dir/lammps.pt > gatmp
	python $cspline gatmp $dir/tests/gbabspline; rm -f gatmp
	echo "g_bab..."

	LB_G_BBB=$(($LE_G_BAB + 4));
	LE_G_BBB=$(($LB_G_BBB + $N_G_BBB - 1));
	sed -n "$(($LB_G_BBB-2))","$LE_G_BBB"p $dir/lammps.pt > gbtmp
	python $cspline gbtmp $dir/tests/gbbbspline; rm -f gbtmp
	echo "g_bbb..."

	blue="#0000FF"
	red="#FF0000"
	purple="#800080"

	echo "set terminal postscript enhanced color font ',5'
	set encoding iso_8859_1
	set output '$dir/tests/SPLINES.ps'
	set size 8,6
	set multiplot layout 4,4
	
	##set title '$elem1-$elem1 Pair Potential'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol f}(r) [eV]'
	plot '$dir/lammps.pt' every ::$(($LB_PAIR_AA-2))::$LE_PAIR_AA u 1:2 w p pt 7 lc rgb '$red' notitle,\
	'$dir/tests/phiaaspline' u 1:2 w l lt 1 lc rgb '$red' notitle
	
	##set title '$elem1-$elem2 Pair Potential'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol f}(r) [eV]'
	plot '$dir/lammps.pt' every ::$(($LB_PAIR_AB-2))::$LE_PAIR_AB u 1:2 w p pt 7 lc rgb '$purple' notitle,\
	'$dir/tests/phiabspline' u 1:2 w l lt 1 lc rgb '$purple' notitle
	
	##set title '$elem2-$elem2 Pair Potential'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol f}(r) [eV]'
	plot '$dir/lammps.pt' every ::$(($LB_PAIR_BB-2))::$LE_PAIR_BB u 1:2 w p pt 7 lc rgb '$blue' notitle,\
	'$dir/tests/phibbspline' u 1:2 w l lt 1 lc rgb '$blue' notitle
	
	##set title '$elem1 Density Function'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol r}(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_RHO_A-2))::$LE_RHO_A u 1:2 w p pt 7 lc rgb '$red' notitle,\
	'$dir/tests/rhoaspline' u 1:2 w l lt 1 lc rgb '$red' notitle
	
	##set title '$elem1-$elem1 SW Weight Function'
	set xlabel 'r [{\305}]'
	set ylabel 'f(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_F_AA-2))::$LE_F_AA u 1:2 w p pt 7 lc rgb '$red' notitle,\
	'$dir/tests/faaspline' u 1:2 w l lt 1 lc rgb '$red' notitle
	
	##set title '$elem1-$elem2 SW Weight Function'
	set xlabel 'r [{\305}]'
	set ylabel 'f(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_F_AB-2))::$LE_F_AB u 1:2 w p pt 7 lc rgb '$purple' notitle,\
	'$dir/tests/fabspline' u 1:2 w l lt 1 lc rgb '$purple' notitle
	
	##set title '$elem2-$elem2 SW Weight Function'
	set xlabel 'r [{\305}]'
	set ylabel 'f(r)'
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_F_BB-2))::$LE_F_BB u 1:2 w p pt 7 lc rgb '$blue' notitle,\
	'$dir/tests/fbbspline' u 1:2 w l lt 1 lc rgb '$blue' notitle
	
	set xrange [*:*]
	##set title '$elem2 Density Function'
	set xlabel 'r [{\305}]'
	set ylabel '{/Symbol r}(r)'
	plot '$dir/lammps.pt' every ::$(($LB_RHO_B-2))::$LE_RHO_B u 1:2 w p pt 7 lc rgb '$blue' notitle,\
	'$dir/tests/rhobspline' u 1:2 w l lt 1 lc rgb '$blue' notitle
	
	##set title '$elem1-$elem1-$elem1 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_AAA-2))::$LE_G_AAA u 1:2 w p pt 7 lc rgb '$red' notitle,\
	'$dir/tests/gaaaspline' u 1:2 w l lt 1 lc rgb '$red' notitle

	##set title '$elem1-$elem1-$elem2 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_AAB-2))::$LE_G_AAB u 1:2 w p pt 7 lc rgb '$red' notitle,\
	'$dir/tests/gaabspline' u 1:2 w l lt 1 lc rgb '$purple' notitle

	##set title '$elem1-$elem2-$elem2 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_ABB-2))::$LE_G_ABB u 1:2 w p pt 7 lc rgb '$red' notitle,\
	'$dir/tests/gabbspline' u 1:2 w l lt 1 lc rgb '$blue' notitle

	##set title '$elem1 Embedding Function'
	set xlabel '{/Symbol r}'
	set ylabel 'U({/Symbol r}) [eV]'
	##set xrange [0:1]
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_EMBED_A-2))::$LE_EMBED_A u 1:2 w p pt 7 lc rgb '$red' notitle,\
	'$dir/tests/Uaspline' u 1:2 w l lt 1 lc rgb '$red' notitle

	##set title '$elem2-$elem1-$elem1 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_BAA-2))::$LE_G_BAA u 1:2 w p pt 7 lc rgb '$blue' notitle,\
	'$dir/tests/gbaaspline' u 1:2 w l lt 1 lc rgb '$red' notitle
	
	##set title '$elem2-$elem1-$elem2 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_BAB-2))::$LE_G_BAB u 1:2 w p pt 7 lc rgb '$blue' notitle,\
	'$dir/tests/gbabspline' u 1:2 w l lt 1 lc rgb '$purple' notitle
	
	##set title '$elem2-$elem2-$elem2 SW Angular Function'
	set xlabel 'cos {/Symbol q}'
	set ylabel 'g(cos {/Symbol q})'
	set xrange [-1:1]
	plot '$dir/lammps.pt' every ::$(($LB_G_BBB-2))::$LE_G_BBB u 1:2 w p pt 7 lc rgb '$blue' notitle,\
	'$dir/tests/gbbbspline' u 1:2 w l lt 1 lc rgb '$blue' notitle
	
	##set title '$elem2 Embedding Function'
	set xlabel '{/Symbol r}'
	set ylabel 'U({/Symbol r}) [eV]'
	##set xrange [0:1]
	set xrange [*:*]
	plot '$dir/lammps.pt' every ::$(($LB_EMBED_B-2))::$LE_EMBED_B u 1:2 w p pt 7 lc rgb '$blue' notitle,\
	'$dir/tests/Ubspline' u 1:2 w l lt 1 lc rgb '$blue' notitle

	unset multiplot" | gnuplot

$ps2png $dir/tests/SPLINES.ps $dir/tests/SPLINES.png
#rm -f $dir/tests/*spline ## clean up

fi
########################################################################################
#
########################################################################################
echo "plotting overall performance..."
########################################################################################
#	RADIAL PERFORMANCE PLOT
########################################################################################

numprops=25
dtheta=`echo "360/($numprops)" | bc -l`
theta=0

echo "## property, angle, error(%), potval, dftval" > $dir/tests/radperf.dat

# D03 lattice constant: 0
error=`python -c "print 100*(${D03PLAT[0]}-$D03lat)/$D03lat"`
lab0=`printf '%.2f' $error`
echo "a_hcp $theta $error ${D03PLAT[0]} $D03lat" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
D03eqp=`printf '%.3f' ${D03PLAT[0]}`
D03laterr=$lab0 

# G1 lattice constant: 1
error=`python -c "print 100*(${G1PLAT[0]}-$G1lat)/$G1lat"`
lab1=`printf '%.2f' $error`
echo "a_bcc $theta $error ${G1PLAT[0]} $G1lat" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
G1eqp=`printf '%.3f' ${G1PLAT[0]}` 
G1laterr=$lab1

# SQS lattice constant: 2
error=`python -c "print 100*(${SQS7525PLAT[0]}-$SQS7525lat)/$SQS7525lat"`
lab2=`printf '%.2f' $error`
echo "a_SQS $theta $error ${SQS7525PLAT[0]} $SQS7525lat" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
SQS7525eqp=`printf '%.3f' ${SQS7525PLAT[0]}` 
SQS7525laterr=$lab2

# APP lattice constant: 3
error=`python -c "print 100*(${APPPLAT[0]}-$APPlat)/$APPlat"`
lab3=`printf '%.2f' $error`
echo "a_APP $theta $error ${APPPLAT[0]} $APPlat" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
APPeqp=`printf '%.3f' ${APPPLAT[0]}` 
APPlaterr=$lab3

# B2 lattice constant: not on circle
error=`python -c "print 100*(${B2PLAT[0]}-$B2lat)/$B2lat"`
B2laterr=`printf '%.2f' $error`
B2eqp=`printf '%.3f' ${B2PLAT[0]}` 

# APP b/a ratio: 4
error=`python -c "print 100*($APPboap-$APPboa)/$APPboa"`
lab4=`printf '%.2f' $error`
echo "boa_APP $theta $error $APPboap $APPboa" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
APPboap=`printf '%.3f' $APPboap` 
APPboa=`printf '%.3f' $APPboa` 
APPboaerr=$lab4

# APP c/a ratio: 5
error=`python -c "print 100*($APPcoap-$APPcoa)/$APPcoa"`
lab5=`printf '%.2f' $error`
echo "a_APP $theta $error $APPcoap $APPcoa" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
APPcoap=`printf '%.3f' $APPcoap` 
APPcoa=`printf '%.3f' $APPcoa` 
APPcoaerr=$lab5

# AP lat : 6
error=`python -c "print 100*(${APPLAT[0]}-$APlat)/$APlat"`
lab6=`printf '%.2f' $error`
echo "a_AP $theta $error ${APPLAT[0]} $APlat" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
APlat=`printf '%.3f' $APlat`
APeqp=`printf '%.3f' ${APPLAT[0]}`
APlaterr=$lab6 

# omega lattice constant: 7
error=`python -c "print 100*(${omgPLAT[0]}-$omglat)/$omglat"`
lab7=`printf '%.2f' $error`
echo "a_omg $theta $error ${omgPLAT[0]} $omglat" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
omgeqp=`printf '%.3f' ${omgPLAT[0]}` 
omglat=`printf '%.3f' $omglat` 
omglaterr=$lab7

# omega c/a ratio: 8
error=`python -c "print 100*($omgcoap-$omgcoa)/$omgcoa"`
lab8=`printf '%.2f' $error`
echo "a_omg $theta $error $omgcoap $omgcoa" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
omgcoap=`printf '%.3f' $omgcoap` 
omgcoa=`printf '%.3f' $omgcoa` 
omgcoaerr=$lab8

# D03 ediff: 9
pdiff=`python -c "print 1000*($D03pote-$G1pote)"`
D03pd=`printf '%.3f' $pdiff`
ddiff=`python -c "print 1000*($D03ep-$G1ep)"`
D03dd=`printf '%.3f' $ddiff`
error=`python -c "print 100*($pdiff-$ddiff)/$ddiff"`
lab9=`printf '%.2f' $error`
echo "D03-G1 $theta $error $pdiff $ddiff" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
D03edifferr=$lab9

# B2 ediff: not on circle
pdiff=`python -c "print 1000*($B2pote-$G1pote)"`
B2pd=`printf '%.3f' $pdiff`
ddiff=`python -c "print 1000*($B2ep-$G1ep)"`
B2dd=`printf '%.3f' $ddiff`
error=`python -c "print 100*($pdiff-$ddiff)/$ddiff"`
B2edifferr=`printf '%.2f' $error`

# SQS ediff: 10
pdiff=`python -c "print 1000*($SQS7525pote-$G1pote)"`
SQS7525pd=`printf '%.3f' $pdiff`
ddiff=`python -c "print 1000*($SQS7525ep-$G1ep)"`
SQS7525dd=`printf '%.3f' $ddiff`
error=`python -c "print 100*($pdiff-$ddiff)/$ddiff"`
lab10=`printf '%.2f' $error`
echo "SQS-G1 $theta $error $pdiff $ddiff" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
SQS2575edifferr=$lab10

# APP ediff: 11
pdiff=`python -c "print 1000*($APPpote-$G1pote)"`
APPpd=`printf '%.3f' $pdiff`
ddiff=`python -c "print 1000*($APPep-$G1ep)"`
APPdd=`printf '%.3f' $ddiff`
error=`python -c "print 100*($pdiff-$ddiff)/$ddiff"`
lab11=`printf '%.2f' $error`
echo "app-G1 $theta $error $pdiff $ddiff" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
APPedifferr=$lab11

# AP ediff: 12
pdiff=`python -c "print 1000*($APpote-$G1pote)"`
APpd=`printf '%.3f' $pdiff`
ddiff=`python -c "print 1000*($APep-$G1ep)"`
APdd=`printf '%.3f' $ddiff`
error=`python -c "print 100*($pdiff-$ddiff)/$ddiff"`
lab12=`printf '%.2f' $error`
echo "ap-G1 $theta $error $pdiff $ddiff" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
APedifferr=$lab12

# omega ediff: 13
pdiff=`python -c "print 1000*($omgpote-$G1pote)"`
omgpd=`printf '%.3f' $pdiff`
ddiff=`python -c "print 1000*($omgep-$G1ep)"`
omgdd=`printf '%.3f' $ddiff`
error=`python -c "print 100*($pdiff-$ddiff)/$ddiff"`
lab13=`printf '%.2f' $error`
echo "omg-G1 $theta $error $pdiff $ddiff" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
omgedifferr=$lab13

# D03 Bulk modulus: 14
error=`python -c "print 100*($D03bulkp-$D03bulkd)/$D03bulkd"`
echo "B_D03 $theta $error $D03bulkp $D03bulkd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab14=`printf '%.2f' $error`
D03bulkp=`printf '%.0f' $D03bulkp` 
D03bulkd=`printf '%.0f' $D03bulkd` 
D03Berr=$lab14

# B2 Bulk modulus: not on circle
error=`python -c "print 100*($B2bulkp-$B2bulkd)/$B2bulkd"`
B2Berr=`printf '%.2f' $error`
B2bulkp=`printf '%.0f' $B2bulkp` 
B2bulkd=`printf '%.0f' $B2bulkd` 

# G1 Bulk modulus: 15
error=`python -c "print 100*($G1bulkp-$G1bulkd)/$G1bulkd"`
echo "B_G1 $theta $error $G1bulkp $G1bulkd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab15=`printf '%.2f' $error`
G1bulkp=`printf '%.0f' $G1bulkp` 
G1bulkd=`printf '%.0f' $G1bulkd` 
G1Berr=$lab15

# SQS Bulk modulus: 16
error=`python -c "print 100*($SQS7525bulkp-$SQS7525bulkd)/$SQS7525bulkd"`
echo "B_SQS $theta $error $SQS7525bulkp $SQS7525bulkd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab16=`printf '%.2f' $error`
SQS7525bulkp=`printf '%.0f' $SQS7525bulkp` 
SQS7525bulkd=`printf '%.0f' $SQS7525bulkd` 
SQS7525Berr=$lab16

# APP Bulk modulus: 17
error=`python -c "print 100*($APPbulkp-$APPbulkd)/$APPbulkd"`
echo "B_APP $theta $error $APPbulkp $APPbulkd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab17=`printf '%.2f' $error`
APPbulkp=`printf '%.0f' $APPbulkp` 
APPbulkd=`printf '%.0f' $APPbulkd` 
APPBerr=$lab17

# AP Bulk modulus: 18
error=`python -c "print 100*($APbulkp-$APbulkd)/$APbulkd"`
echo "B_AP $theta $error $APbulkp $APbulkd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab18=`printf '%.2f' $error`
APbulkp=`printf '%.0f' $APbulkp` 
APbulkd=`printf '%.0f' $APbulkd` 
APBerr=$lab18

# omega Bulk modulus: 19
error=`python -c "print 100*($omgbulkp-$omgbulkd)/$omgbulkd"`
echo "B_omg $theta $error $omgbulkp $omgbulkd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab19=`printf '%.2f' $error`
omgbulkp=`printf '%.0f' $omgbulkp`
omgbulkd=`printf '%.0f' $omgbulkd`
omgBerr=$lab19

# G1 C11: 20
error=`python -c "print 100*($G1C11-$C11G1d)/float($C11G1d)"`
echo "C11_G1 $theta $error $G1C11 $C11G1d" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab20=`printf '%.2f' $error`
G1C11err=$lab21

# G1 C12: 21
error=`python -c "print 100*($G1C12-$C12G1d)/float($C12G1d)"`
echo "C12_G1 $theta $error $G1C12 $C12G1d" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab21=`printf '%.2f' $error`
G1C12err=$lab22

# G1 C44: 22
error=`python -c "print 100*($G1C44-$C44G1d)/float($C44G1d)"`
echo "C44_G1 $theta $error $G1C44 $C44G1d" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab22=`printf '%.2f' $error`
G1C44err=$lab23

# HCP-Ti prism easy: 23
error=`python -c "print 100*($hcpPrisEasyp-$hcpPrisEasyd)/$hcpPrisEasyd"`
echo "prism_e $theta $error $hcpPrisEasyp $hcpP4isEasyd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab23=`printf '%.2f' $error`

# HCP-Ti prism hard: 24
error=`python -c "print 100*($hcpPrisHardp-$hcpPrisHardd)/$hcpPrisHardd"`
echo "prism_e $theta $error $hcpPrisHardp $hcpP4isHardd" >> $dir/tests/radperf.dat
theta=`echo "$theta + $dtheta" | bc -l`
lab24=`printf '%.2f' $error`


awk 'NR==2{print}' $dir/tests/radperf.dat >> $dir/tests/radperf.dat



labrad=230
echo "set terminal postscript enhanced color
set encoding iso_8859_1
set polar
set angles degrees
set size square
set xrange [-250:250]
set yrange [-250:250]

#unset tics
#unset grid
set output '$dir/tests/radperf_small.ps'
plot '$dir/tests/radperf.dat' u 2:(\$3+100) w lp lc rgb 'red' lw 3 lt 1 pt 7 notitle,\
#'' u 2:(\$3+100) w filledcurves above r=100 lc rgb '#A9A9A9' lw 3 lt 1 notitle,\
#'' u 2:(\$3+100) w filledcurves below r=100 lc rgb '#A9A9A9' lw 3 lt 1 notitle,\
100 w l lc 0 lw 3 notitle

unset border
set output '$dir/tests/radperf.ps'
set rrange [0:200]
set xrange [-250:250]# noextend
set yrange [-250:250]# noextend
#unset xtics
unset ytics
unset y2tics
set rtics format \"\"
set style fill solid 0.5
unset grid
set grid polar $dtheta
set xtics (\"-50\\%%\" 50, \"\" 100, \"50\\%%\" 150, \"\" 200) nomirror axis

set label \"$lab0  \% \n a_{D0_{3}}\"				 at ($labrad*cos(0)),($labrad*sin(0)) center front
set label \"$lab1  \% \n a_{G1}\" 				 at ($labrad*cos($dtheta)),($labrad*sin($dtheta)) center front
set label \"$lab2  \% \n a_{SQS}\" 				 at ($labrad*cos(2*$dtheta)),($labrad*sin(2*$dtheta)) center front
set label \"$lab3  \% \n a_{{/Symbol a}^{,,}}\" 		 at ($labrad*cos(3*$dtheta)),($labrad*sin(3*$dtheta)) center front
set label \"$lab4  \% \n (b/a)_{{/Symbol a}^{,,}}\" 		 at ($labrad*cos(4*$dtheta)),($labrad*sin(4*$dtheta)) center front
set label \"$lab5  \% \n (c/a)_{{/Symbol a}^{,,}}\" 		 at ($labrad*cos(5*$dtheta)),($labrad*sin(5*$dtheta)) center front
set label \"$lab6  \% \n a_{{/Symbol a}^{,,}}\" 		 at ($labrad*cos(6*$dtheta)),($labrad*sin(6*$dtheta)) center front
set label \"$lab7  \% \n a_{{/Symbol w}}\" 			 at ($labrad*cos(7*$dtheta)),($labrad*sin(7*$dtheta)) center front
set label \"$lab8  \% \n (c/a)_{{/Symbol w}}\" 			 at ($labrad*cos(8*$dtheta)),($labrad*sin(8*$dtheta)) center front
set label \"$lab9  \% \n {/Symbol D}E_{D03}\" 			 at ($labrad*cos(9*$dtheta)),($labrad*sin(9*$dtheta)) center front
set label \"$lab10 \% \n {/Symbol D}E_{SQS}\" 			 at ($labrad*cos(10*$dtheta)),($labrad*sin(10*$dtheta)) center front
set label \"$lab11 \% \n {/Symbol D}E_{{/Symbol a}^{,,}}\" 	 at ($labrad*cos(11*$dtheta)),($labrad*sin(11*$dtheta)) center front
set label \"$lab12 \% \n {/Symbol D}E_{{/Symbol a}^{,}}\" 	 at ($labrad*cos(12*$dtheta)),($labrad*sin(12*$dtheta)) center front
set label \"$lab13 \% \n {/Symbol D}E_{{/Symbol w}}\" 		 at ($labrad*cos(13*$dtheta)),($labrad*sin(13*$dtheta)) center front
set label \"$lab14 \% \n B_{D03}\" 				 at ($labrad*cos(14*$dtheta)),($labrad*sin(14*$dtheta)) center front
set label \"$lab15 \% \n B_{G1}\" 				 at ($labrad*cos(15*$dtheta)),($labrad*sin(15*$dtheta)) center front
set label \"$lab16 \% \n B_{SQS}\" 				 at ($labrad*cos(16*$dtheta)),($labrad*sin(16*$dtheta)) center front
set label \"$lab17 \% \n B_{{/Symbol a}^{,,}}\"  		 at ($labrad*cos(17*$dtheta)),($labrad*sin(17*$dtheta)) center front
set label \"$lab18 \% \n B_{{/Symbol a}^{,}}\"  		 at ($labrad*cos(18*$dtheta)),($labrad*sin(18*$dtheta)) center front
set label \"$lab19 \% \n B_{{/Symbol w}}\"			 at ($labrad*cos(19*$dtheta)),($labrad*sin(19*$dtheta)) center front
set label \"$lab20 \% \n C_{11}^{G1}\"				 at ($labrad*cos(20*$dtheta)),($labrad*sin(20*$dtheta)) center front
set label \"$lab21 \% \n C_{12}^{G1}\"				 at ($labrad*cos(21*$dtheta)),($labrad*sin(21*$dtheta)) center front
set label \"$lab22 \% \n C_{44}^{G1}\"				 at ($labrad*cos(22*$dtheta)),($labrad*sin(22*$dtheta)) center front
set label \"$lab23 \% \n {/Symbol g}_{pris}^{0.5-E}\"		 at ($labrad*cos(23*$dtheta)),($labrad*sin(23*$dtheta)) center front
set label \"$lab24 \% \n {/Symbol g}_{pris}^{0.5-H}\"		 at ($labrad*cos(24*$dtheta)),($labrad*sin(24*$dtheta)) center front

set key at 250,250

plot '$dir/tests/radperf.dat' u 2:(\$3+100) w lp lc rgb 'red' lw 3 lt 1 pt 7 notitle,\
# '' u 2:(\$3+100) w filledcurves above r=100 lc rgb '#A9A9A9' lw 3 lt 1 notitle,\
# '' u 2:(\$3+100) w filledcurves below r=100 lc rgb '#A9A9A9' lw 3 lt 1 notitle,\
 100 w l lc 0 lw 3 notitle,1/0 w l lc 0 lw 3 lt 2 title 'PAW-PBE',\
 1/0 w l lc rgb 'red' lw 3 lt 1 title '$typ'

unset polar" > $dir/tests/radperf_plot
echo "load '$dir/tests/radperf_plot'" | gnuplot
#cat $dir/tests/radperf.dat

$ps2png $dir/tests/radperf.ps $dir/tests/radperf.png
ps2png $dir/tests/radperf_small.ps $dir/tests/radperf_small.png 2>&1 /dev/null

########################################################################################
#
########################################################################################


########################################################################################
#    ERROR PLOTS
########################################################################################
echo "plotting fitting errors..."
# get total fitting error
if [ -f $dir/meamzilla.out ]; then
	fitError=`grep 'Final total sum' $dir/meamzilla.out | awk '{printf "%4f", $9}'`
elif [ -f $dir/meamzilla_powell.out ]; then
	fitError=`grep 'Final total sum' $dir/meamzilla_powell.out | awk '{printf "%4f", $9}'`
elif [ -f $dir/meamzilla_genalg.out ]; then
	fitError=`grep 'Final total sum' $dir/meamzilla_genalg.out | awk '{printf "%4f", $9}'`
else
	fitError="Not Found"
fi

#create force plot file
awk 'NR>1{ 
	fdft=0; feam=0; fdot=0;
	for (i=1; i<3; i++) {
		getline
		fdft += $6**2
		feam += $5**2
		fdot += $5*$6
	} fdft = sqrt(fdft); feam = sqrt(feam);
	if (fdft == 0) { 
		ferr = 1000
		fang = 0
	}
	else { 
		ferr = 100*sqrt((feam-fdft)**2)/sqrt(fdft**2)
		fdot = fdot/sqrt((feam*fdft)**2)
		fang = 180*atan2(sqrt(sqrt((1-fdot**2)**2), fdot)/atan2(0,-1); 	# inverse cosine!
	}
	{printf "%9f %9f %9f %9f\n", fdft, feam, ferr, fang}}' $dir/data.force > forcetmp

#create energy plot file
awk 'NR>1{
	print $7, $6}' $dir/data.energy > energytmp	

#create stress plot file
awk 'NR>1{
	pdft=0; peam=0;
	for (i=1; i<3; i++) {
		getline
		pdft -= $6
		peam -= $5
	} pdft = sqrt(pdft**2)/3; peam = sqrt(peam**2)/3;
	if (pdft == 0) {
		perr = 1000
	}
	else {
		perr = 100*sqrt((peam-pdft)**2)/sqrt(pdft**2)
	}
	sdft=0; seam=0;
	for (j=1; j<3; j++) {
		getline
		sdft += $6
		seam += $5
	} sdft = sqrt(sdft**2)/3; seam = sqrt(seam**2)/3;
	if (sdft == 0) {
		serr = 1000
	}
	else {
		serr = 100*sqrt((seam-sdft)**2)/sqrt(sdft**2)
	}
	{printf "%9f %9f %9f %9f %9f %9f\n", pdft, peam, perr, sdft, seam, serr}}' $dir/data.stress > stresstmp

echo "set terminal postscript enhanced color
set encoding iso_8859_1
set output '$dir/tests/ERRORS.ps'

set multiplot layout 2,3

unset key
f(x) = x

set lmargin 3.5
set rmargin 3.5
set tmargin 3.5
set bmargin 3.5

set title 'E_{MEAM} vs E_{DFT}'
set xtics -8.5,.5,-5.5
set ytics -8.5,.5,-5.5
set xlabel 'E_{DFT} (eV/atom)' offset 2.5
set ylabel 'E_{MEAM} (eV/atom)' offset 1.5
plot 'energytmp' u 1:2 w p pt 13 lc rgb 'dark-orange' ps 1,\
f(x) lc 0

set xtics auto
set ytics auto

set title 'Pressure Deviations'
set yrange [0:120]
set xlabel '|P_{DFT}| (eV/{\\305}^{3})' offset 2.5
set ylabel '|P_{MEAM}| error (%)' offset 2.5
plot 'stresstmp' u 1:3 w p pt 9 lc rgb 'dark-green' ps 1

set title 'Shear Stress Deviations'
set xlabel 'avg. |{/Symbol t}_{DFT}| (eV/{\\305}^{3})' offset 2.5
set ylabel 'avg. |{/Symbol t}_{MEAM}| error (%)' offset 2.5
plot 'stresstmp' u 4:6 w p pt 9 lc rgb 'dark-green' ps 1

set autoscale y
set xlabel '|F_{DFT}| (eV/{\\305})' offset 2.5
set title 'F_{MEAM} vs F_{DFT}'
set ylabel '|F_{MEAM}| (eV/{\\305})' offset 1.5
plot 'forcetmp' u 1:2 w p pt 7 lc rgb 'blue' ps 0.5,\
f(x) lc 0

set yrange [0:120]
set title 'Magnitude Deviations'
set ylabel '|F_{MEAM}| error (%)' offset 2.5
plot 'forcetmp' u 1:3 w p pt 7 lc rgb 'blue' ps 0.5

set title 'Angular Deviations'
set ylabel '{/Symbol q}_{err} (degrees)' offset 2.5
plot 'forcetmp' u 1:4 w p pt 7 lc rgb 'blue' ps 0.5" | gnuplot

$ps2png $dir/tests/ERRORS.ps $dir/tests/ERRORS.png

rm stresstmp energytmp forcetmp


########################################################################################
#	ENERGY-VOLUME AND PRESSURE-VOLUME CURVES
########################################################################################
echo "plotting E-V and P-V curves..."

EVminRANGE=14
#if [ -e $dftdat/$elem/hcp/evsv.dat ]; then
#	EVmaxRANGE=`tail -1 $dftdat/$elem/hcp/evsv.dat | awk '{print $1}'`
#else
#	EVmaxRANGE=21
#fi
EVmaxRANGE=21

EVminRANGEmeta=14
EVmaxRANGEmeta=21

echo "set terminal postscript enhanced color
set encoding iso_8859_1
set title 'Energy - Volume curves for Ti_{3}Nb {/Symbol b} Candidates'
set xlabel 'Volume [{\305}^{3}/atom]'
set ylabel 'Energy [meV/atom]'
set xrange [$EVminRANGE:$EVmaxRANGE]
set key top center

set output '$dir/tests/E_VS_V_BETA.ps'
plot '$dir/tests/D03/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'red' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/D03/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb 'red' pt 6 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 6 ps 1 lc rgb 'red' title 'D0_{3}',\
'$dir/tests/G1/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'dark-green' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/G1/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p pt 4 ps 1 lc rgb 'dark-green' notitle,\
1/0 w lp lt 1 lw 3 pt 4 ps 1 lc rgb 'dark-green' title 'G1',\
'$dir/tests/L60Ti3Nb/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb '#696969' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/L60Ti3Nb/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb '#696969' pt 8 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 8 ps 1 lc rgb '#696969' title 'L6_{0}',\
'$dir/tests/SQS7525/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'blue' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/SQS7525/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb 'blue' pt 10 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 10 ps 1 lc rgb 'blue' title 'SQS',\
'$dir/tests/SS/evsv_bcc_7525.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'black' lt 1 lw 2 title '{/Symbol b} S.S.',\
'$dir/tests/SS/evsv_hcp_7525.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'black' lt 0 lw 2 title '{/Symbol a} S.S.'

set xrange [$EVminRANGEmeta:$EVmaxRANGEmeta]
set title 'Energy - Volume curves for Ti_{3}Nb Metastable Phases'
set output '$dir/tests/E_VS_V_META.ps'
plot '$dir/tests/APP/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'dark-orange' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/APP/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb 'dark-orange' pt 12 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 12 ps 1 lc rgb 'dark-orange' title '{/Symbol a}^{,,}',\
'$dir/tests/A15Ti3Nb/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb '#32CD32' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/A15Ti3Nb/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb '#32CD32' pt 12 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 7 ps 1 lc rgb '#32CD32' title 'A15',\
'$dir/tests/omega/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'dark-turquoise' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/omega/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb 'dark-turquoise' pt 14 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 14 ps 1 lc rgb 'dark-turquoise' title '{/Symbol w}',\
'$dir/tests/AP/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb '#9932CC' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/AP/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb '#9932CC' pt 8 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 8 ps 1 lc rgb '#9932CC' title '{/Symbol a}^{,}',\
'$dir/tests/D019Ti3Nb/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb '#FF00FF' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/D019/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb '#FF00FF' pt 10 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 49 ps 1 lc rgb '#FF00FF' title '{/Symbol a}-D0_{19}',\
'$dir/tests/G1/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'dark-green' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/G1/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb 'dark-green' pt 4 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 4 ps 1 lc rgb 'dark-green' title '{/Symbol b}_{G1}',\
'$dir/tests/SS/evsv_bcc_7525.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'black' lt 1 lw 2 title '{/Symbol b} S.S.',\
'$dir/tests/SS/evsv_hcp_7525.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'black' lt 0 lw 2 title '{/Symbol a} S.S.'

set key top center
set title 'Energy - Volume curves for Pure $elem1'
set output '$dir/tests/E_VS_V_Ti.ps'
plot '$dir/tests/HCPti/evsv.dat' u 1:(1000*(\$2+(-1)*$HCPtipote)) w l lc rgb '#228B22' lt 1 lw 3 notitle,\
'$dftdat/$elem1/hcp/evsv.dat' u 1:(1000*(\$2+(-1)*$hcpTiep)) w p lc rgb '#228B22' pt 8 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 8 ps 1 lc rgb '#228B22' title 'hcp $elem1',\
'$dir/tests/OMGti/evsv.dat' u 1:(1000*(\$2+(-1)*$HCPtipote)) w l lc rgb '#B22222' lt 1 lw 3 notitle,\
'$dftdat/$elem1/omega/evsv.dat' u 1:(1000*(\$2+(-1)*$hcpTiep)) w p lc rgb '#B22222' pt 10 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 10 ps 1 lc rgb '#B22222' title '{/Symbol w}-$elem1',\
'$dir/tests/FCCti/evsv.dat' u 1:(1000*(\$2+(-1)*$HCPtipote)) w l lc rgb '#8A2BE2' lt 1 lw 3 notitle,\
'$dftdat/$elem1/fcc/evsv.dat' u 1:(1000*(\$2+(-1)*$hcpTiep)) w p lc rgb '#8A2BE2' pt 16 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 16 ps 1 lc rgb '#8A3BE2' title 'fcc $elem1',\
'$dir/tests/A15ti/evsv.dat' u 1:(1000*(\$2+(-1)*$HCPtipote)) w l lc rgb '#FF1493' lt 1 lw 3 notitle,\
'$dftdat/$elem1/A15/evsv.dat' u 1:(1000*(\$2+(-1)*$hcpTiep)) w p lc rgb '#FF1493' pt 14 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 14 ps 1 lc rgb '#FF1493' title 'A15 $elem1',\
'$dir/tests/BCCti/evsv.dat' u 1:(1000*(\$2+(-1)*$HCPtipote)) w l lc rgb '#1E90FF' lt 1 lw 3 notitle,\
'$dftdat/$elem1/bcc/evsv.dat' u 1:(1000*(\$2+(-1)*$hcpTiep)) w p lc rgb '#1E90FF' pt 12 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 12 ps 1 lc rgb '#1E90FF' title 'bcc $elem1'

set key top center
set xrange [15:22]
set title 'Energy - Volume curves for Pure $elem2'
set output '$dir/tests/E_VS_V_Nb.ps'
plot '$dir/tests/BCCnb/evsv.dat' u 1:(1000*(\$2+(-1)*$BCCnbpote)) w l lc rgb 'dark-violet' lt 1 lw 3 notitle,\
'$dftdat/$elem2/bcc/evsv.dat' u 1:(1000*(\$2+(-1)*$bccNbep)) w p lc rgb 'dark-violet' pt 4 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 4 ps 1 lc rgb 'dark-violet' title 'bcc $elem2',\
'$dir/tests/A15nb/evsv.dat' u 1:(1000*(\$2+(-1)*$BCCnbpote)) w l lc rgb '#4682B4' lt 1 lw 3 notitle,\
'$dftdat/$elem2/A15/evsv.dat' u 1:(1000*(\$2+(-1)*$bccNbep)) w p lc rgb '#4682B4' pt 6 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 6 ps 1 lc rgb '#4682B4' title 'A15 $elem2',\
'$dir/tests/FCCnb/evsv.dat' u 1:(1000*(\$2+(-1)*$BCCnbpote)) w l lc rgb '#FF1493' lt 1 lw 3 notitle,\
'$dftdat/$elem2/fcc/evsv.dat' u 1:(1000*(\$2+(-1)*$bccNbep)) w p lc rgb '#FF1493' pt 16 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 16 ps 1 lc rgb '#FF1493' title 'fcc $elem2',\
'$dir/tests/HCPnb/evsv.dat' u 1:(1000*(\$2+(-1)*$BCCnbpote)) w l lc rgb '#228B22' lt 1 lw 3 notitle,\
'$dftdat/$elem2/hcp/evsv.dat' u 1:(1000*(\$2+(-1)*$bccNbep)) w p lc rgb '#228B22' pt 8 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 8 ps 1 lc rgb '#228B22' title 'hcp $elem2',\
'$dir/tests/OMGnb/evsv.dat' u 1:(1000*(\$2+(-1)*$BCCnbpote)) w l lc rgb '#32CD32' lt 1 lw 3 notitle,\
'$dftdat/$elem2/omega/evsv.dat' u 1:(1000*(\$2+(-1)*$bccNbep)) w p lc rgb '#32CD32' pt 12 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 12 ps 1 lc rgb '#32CD32' title '{/Symbol w}_{Ti} $elem2'

set key top center
set xrange [$EVminRANGE:$EVmaxRANGE]
set title 'Energy - Volume curves for TiNb Phases'
set output '$dir/tests/E_VS_V_5050.ps'
plot '$dir/tests/B2/evsv.dat' u 1:(1000*(\$2+(-1)*$B2pote)) w l lc rgb 'orange-red' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/TiNb/B2/evsv.dat' u 1:(1000*(\$2+(-1)*$B2ep)) w p lc rgb 'orange-red' pt 4 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 8 ps 1 lc rgb 'orange-red' title 'B2-{/Symbol b}',\
'$dir/tests/A3/evsv.dat' u 1:(1000*(\$2+(-1)*$B2pote)) w l lc rgb 'dark-green' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/TiNb/A3/evsv.dat' u 1:(1000*(\$2+(-1)*$B2ep)) w p lc rgb 'dark-green' pt 10 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 10 ps 1 lc rgb 'dark-green' title 'A3-{/Symbol a}',\
'$dir/tests/bcc110/evsv.dat' u 1:(1000*(\$2+(-1)*$B2pote)) w l lc rgb 'blue' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/TiNb/bcc110/evsv.dat' u 1:(1000*(\$2+(-1)*$B2ep)) w p lc rgb 'blue' pt 10 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 10 ps 1 lc rgb 'blue' title '{110}-layered {/Symbol b}',\
'$dir/tests/L10/evsv.dat' u 1:(1000*(\$2+(-1)*$B2pote)) w l lc rgb 'purple' lt 1 lw 3 notitle,\
1/0 w lp lt 1 lw 3 pt 8 ps 1 lc rgb 'purple' title 'L1_{0}-FCT',\
'$dir/tests/SQS5050/evsv.dat' u 1:(1000*(\$2+(-1)*$B2pote)) w l lc rgb 'red' lt 1 lw 3 notitle,\
1/0 w p lt 1 lw 3 pt 14 ps 1 lc rgb 'red' title '{/Symbol b} SQS',\
'$dir/tests/SS/evsv_bcc_5050.dat' u 1:(1000*(\$2+(-1)*$B2pote)) w l lc rgb 'black' lt 1 lw 2 title '{/Symbol b} S.S.',\
'$dir/tests/SS/evsv_hcp_5050.dat' u 1:(1000*(\$2+(-1)*$B2pote)) w l lc rgb 'black' lt 0 lw 2 title '{/Symbol a} S.S.'

set key top center
set title 'Energy - Volume curves for TiNb_{3} Phases'
set output '$dir/tests/E_VS_V_2575.ps'
plot '$dir/tests/L60TiNb3/evsv.dat' u 1:(1000*(\$2+(-1)*$L60TiNb3pote)) w l lc rgb '#00008B' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/L60TiNb3/evsv.dat' u 1:(1000*(\$2+(-1)*$L60TiNb3ep)) w p lc rgb '#00008B' pt 8 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 8 ps 1 lc rgb '#00008B' title 'L6_{0}-{/Symbol b}',\
'$dir/tests/A15TiNb3/evsv.dat' u 1:(1000*(\$2+(-1)*$L60TiNb3pote)) w l lc rgb '#008080' lt 1 lw 3 notitle,\
'$dftdat/$elem1-$elem2/A15TiNb3/evsv.dat' u 1:(1000*(\$2+(-1)*$L60TiNb3ep)) w p lc rgb '#008080' pt 8 ps 1 notitle,\
1/0 w lp lt 1 lw 3 pt 7 ps 1 lc rgb '#008080' title 'A15',\
'$dir/tests/SQS2575/evsv.dat' u 1:(1000*(\$2+(-1)*$L60TiNb3pote)) w l lc rgb 'red' lt 1 lw 3 notitle,\
1/0 w p lt 1 lw 3 pt 14 ps 1 lc rgb 'red' title '{/Symbol b} SQS',\
'$dir/tests/SS/evsv_bcc_2575.dat' u 1:(1000*(\$2+(-1)*$L60TiNb3pote)) w l lc rgb 'black' lt 1 lw 2 title '{/Symbol b} S.S.',\
'$dir/tests/SS/evsv_hcp_2575.dat' u 1:(1000*(\$2+(-1)*$L60TiNb3pote)) w l lc rgb 'black' lt 0 lw 2 title '{/Symbol a} S.S.'

set key top center
set title 'Energy - Volume curves for Ti_{2}Nb Phases'
set output '$dir/tests/E_VS_V_6733.ps'
plot '$dir/tests/omgTi2Nb/evsv.dat' u 1:(1000*(\$2+(-1)*$omgTi2Nbpote)) w l lw 2 lc rgb 'red' notitle,\
'$dftdat/Ti-Nb/omgTi2Nb/evsv.dat' u 1:(1000*(\$2+(-1)*$omgTi2Nbep)) w p pt 4 ps 1 lc rgb 'red' notitle,\
'$dir/tests/SS/evsv_bcc_6733.dat' u 1:(1000*(\$2+(-1)*$omgTi2Nbpote)) w l lc rgb 'black' lt 1 lw 2 title '{/Symbol b} S.S.',\
'$dir/tests/SS/evsv_hcp_6733.dat' u 1:(1000*(\$2+(-1)*$omgTi2Nbpote)) w l lc rgb 'black' lt 0 lw 2 title '{/Symbol a} S.S.',\
1/0 w lp lt 1 lw 3 pt 4 lc rgb 'red' title '{/Symbol w}-Ti'

set key top center
set title 'Energy - Volume curves for TiNb_{2} Phases'
set output '$dir/tests/E_VS_V_3367.ps'
plot '$dir/tests/omgTiNb2/evsv.dat' u 1:(1000*(\$2+(-1)*$omgTiNb2pote)) w l lw 2 lc rgb '#0066FF' notitle,\
'$dftdat/Ti-Nb/omgTiNb2/evsv.dat' u 1:(1000*(\$2+(-1)*$omgTiNb2ep)) w p pt 6 ps 1 lc rgb '#0066FF' notitle,\
'$dir/tests/SS/evsv_bcc_3367.dat' u 1:(1000*(\$2+(-1)*$omgTiNb2pote)) w l lc rgb 'black' lt 1 lw 2 title '{/Symbol b} S.S.',\
'$dir/tests/SS/evsv_hcp_3367.dat' u 1:(1000*(\$2+(-1)*$omgTiNb2pote)) w l lc rgb 'black' lt 0 lw 2 title '{/Symbol a} S.S.',\
1/0 w lp lt 1 lw 3 pt 6 lc rgb '#0066FF' title '{/Symbol w}-Ti'" | gnuplot

#'$dir/tests/AP/evsv.dat' u 1:(1000*(\$2+(-1)*$G1pote)) w l lc rgb 'purple' lt 1 lw 3 title '{/Symbol a}^{,} - MEAM',\
#'$dftdat/$elem1-$elem2/AP/evsv.dat' u 1:(1000*(\$2+(-1)*$G1ep)) w p lc rgb 'purple' pt 8 ps 1 title '      DFT ',\

$ps2png $dir/tests/E_VS_V_BETA.ps $dir/tests/E_VS_V_BETA.png
$ps2png $dir/tests/E_VS_V_META.ps $dir/tests/E_VS_V_META.png
$ps2png $dir/tests/E_VS_V_Ti.ps $dir/tests/E_VS_V_Ti.png
$ps2png $dir/tests/E_VS_V_Nb.ps $dir/tests/E_VS_V_Nb.png
$ps2png $dir/tests/E_VS_V_5050.ps $dir/tests/E_VS_V_5050.png
$ps2png $dir/tests/E_VS_V_2575.ps $dir/tests/E_VS_V_2575.png
$ps2png $dir/tests/E_VS_V_6733.ps $dir/tests/E_VS_V_6733.png
$ps2png $dir/tests/E_VS_V_3367.ps $dir/tests/E_VS_V_3367.png

echo "set terminal postscript enhanced color font ',20'
set encoding iso_8859_1
set xlabel 'V/V_{0}'
set ylabel 'Pressure [GPa]'
set xrange [0.5:1.0]
unset key

set tmargin at screen 0.95
set bmargin at screen 0.05
set rmargin at screen 0.90
set lmargin at screen 0.10

set title 'Pressure-Volume curve for D0_{3} Ti_{3}Nb'
set output '$dir/tests/D03_PVSV.ps'
plot '$dir/tests/D03/pvsv.dat' u 1:2 w l lc rgb 'red' lt 1 lw 3,\
'$dftdat/$elem1-$elem2/D03/pvsv.dat' u 1:2 w p pt 6 ps 3 lc rgb 'red' 

set title 'Pressure-Volume curve for G1 Ti_{3}Nb'
set output '$dir/tests/G1_PVSV.ps'
plot '$dir/tests/G1/pvsv.dat' u 1:2 w l lc rgb 'dark-green' lt 1 lw 3,\
'$dftdat/$elem1-$elem2/G1/pvsv.dat' u 1:2 w p pt 4 ps 3 lc rgb 'dark-green' 

set title 'Pressure-Volume curve for SQS Ti_{3}Nb'
set output '$dir/tests/SQS_PVSV.ps'
plot '$dir/tests/SQS7525/pvsv.dat' u 1:2 w l lc rgb 'blue' lt 1 lw 3,\
'$dftdat/$elem1-$elem2/SQS7525/pvsv.dat' u 1:2 w p pt 10 ps 3 lc rgb 'blue'

set title 'Pressure-Volume curve for {/Symbol a}^{,,} Ti_{3}Nb'
set output '$dir/tests/APP_PVSV.ps'
plot '$dir/tests/APP/pvsv.dat' u 1:2 w l lc rgb 'dark-orange' lt 1 lw 3,\
'$dftdat/$elem1-$elem2/APP/pvsv.dat' u 1:2 w p pt 12 ps 3 lc rgb 'dark-orange'

set title 'Pressure-Volume curve for {/Symbol a}^{,} Ti_{3}Nb'
set output '$dir/tests/AP_PVSV.ps'
plot '$dir/tests/AP/pvsv.dat' u 1:2 w l lc rgb 'purple' lt 1 lw 3,\
'$dftdat/$elem1-$elem2/AP/pvsv.dat' u 1:2 w p pt 8 ps 3 lc rgb 'purple'

set title 'Pressure-Volume curve for {/Symbol w} Ti_{3}Nb'
set output '$dir/tests/OMEGA_PVSV.ps'
plot '$dir/tests/omega/pvsv.dat' u 1:2 w l lc rgb 'dark-turquoise' lt 1 lw 3,\
'$dftdat/$elem1-$elem2/omega/pvsv.dat' u 1:2 w p pt 14 ps 3 lc rgb 'dark-turquoise'

set title 'Pressure-Volume curve for B2 TiNb'
set output '$dir/tests/B2_PVSV.ps'
plot '$dir/tests/B2/pvsv.dat' u 1:2 w l lc rgb 'orange-red' lt 1 lw 3,\
'$dftdat/$elem1-$elem2/TiNb/B2/pvsv.dat' u 1:2 w p pt 4 ps 3 lc rgb 'orange-red'

" | gnuplot

$ps2png $dir/tests/D03_PVSV.ps $dir/tests/D03_PVSV.png
$ps2png $dir/tests/G1_PVSV.ps $dir/tests/G1_PVSV.png
$ps2png $dir/tests/SQS_PVSV.ps $dir/tests/SQS7525_PVSV.png

$ps2png $dir/tests/APP_PVSV.ps $dir/tests/APP_PVSV.png
$ps2png $dir/tests/AP_PVSV.ps $dir/tests/AP_PVSV.png
$ps2png $dir/tests/OMEGA_PVSV.ps $dir/tests/OMEGA_PVSV.png

$ps2png $dir/tests/B2_PVSV.ps $dir/tests/B2_PVSV.png

########################################################################################
#
########################################################################################

########################################################################################
#	Pressure dependence of lattice properties
########################################################################################

echo "set terminal postscript enhanced color
set encoding iso_8859_1
set output '$dir/tests/APP_ABCVSP.ps'

set title '{/Symbol a}^{,,} Ti_{3}Nb'
set ylabel 'Lattice constant, a [{\305}]'
set y2label 'Lattice ratios'
set xlabel 'Pressure [GPa]'
set y2tics

set key bottom left

set ytics nomirror
set y2tics 0.2

plot '$dir/tests/APP/abcvp.dat' u 1:2 axes x1y1 w l lw 3 lc rgb 'red' title 'a',\
'' u 1:3 axes x1y2 w l lw 3 lc rgb 'dark-green' title 'b/a',\
'' u 1:4 axes x1y2 w l lw 3 lc rgb 'blue' title 'c/a',\
'$dftdat/Ti-Nb/APP/abcvsp.dat' u 1:2 axes x1y1 w p pt 7 ps 2 lc rgb 'red' notitle,\
'' u 1:3 axes x1y2 w p pt 9 ps 2 lc rgb 'dark-green' notitle,\
'' u 1:4 axes x1y2 w p pt 11 ps 2 lc rgb 'blue' notitle

set output '$dir/tests/HCPti_ABCVSP.ps'
set title '{/Symbol a} Ti'
plot '$dir/tests/HCPti/abcvp.dat' u 1:2 axes x1y1 w l lw 3 lc rgb 'red' title 'a',\
'' u 1:3 axes x1y2 w l lw 3 lc rgb 'blue' title 'c/a',\
'$dftdat/Ti/hcp/abcvsp.dat' u 1:2 axes x1y1 w p pt 7 ps 2 lc rgb 'red' notitle,\
'' u 1:4 axes x1y2 w p pt 9 ps 2 lc rgb 'blue' notitle

set output '$dir/tests/OMGTi2Nb_ABCVSP.ps'
set title '{/Symbol w} Ti_{2}Nb'
plot '$dir/tests/omgTi2Nb/abcvp.dat' u 1:2 axes x1y1 w l lw 3 lc rgb 'red' title 'a',\
'' u 1:3 axes x1y2 w l lw 3 lc rgb 'blue' title 'c/a',\
'$dftdat/Ti-Nb/omgTi2Nb/abcvsp.dat' u 1:2 axes x1y1 w p pt 7 ps 2 lc rgb 'red' notitle,\
'' u 1:4 axes x1y2 w p pt 9 ps 2 lc rgb 'blue' notitle

set output '$dir/tests/OMGTiNb2_ABCVSP.ps'
set title '{/Symbol w} TiNb_{2}'
plot '$dir/tests/omgTiNb2/abcvp.dat' u 1:2 axes x1y1 w l lw 3 lc rgb 'red' title 'a',\
'' u 1:3 axes x1y2 w l lw 3 lc rgb 'blue' title 'c/a',\
'$dftdat/Ti-Nb/omgTiNb2/abcvsp.dat' u 1:2 axes x1y1 w p pt 7 ps 2 lc rgb 'red' notitle,\
'' u 1:4 axes x1y2 w p pt 9 ps 2 lc rgb 'blue' notitle" | gnuplot

$ps2png $dir/tests/APP_ABCVSP.ps $dir/tests/APP_ABCVSP.png
$ps2png $dir/tests/HCPti_ABCVSP.ps $dir/tests/HCPti_ABCVSP.png
$ps2png $dir/tests/OMGTiNb2_ABCVSP.ps $dir/tests/OMGTiNb2_ABCVSP.png
$ps2png $dir/tests/OMGTi2Nb_ABCVSP.ps $dir/tests/OMGTi2Nb_ABCVSP.png

########################################################################################
#	ELASTIC CONSTANT PLOTS
########################################################################################

#define colors
C11COLOR="rgb '#0000FF'"
C12COLOR="rgb '#DC153C'"
C13COLOR="rgb '#008B8B'"
C22COLOR="rgb '#1E90FF'"
C23COLOR="rgb '#FF8C00'"
C33COLOR="rgb '#006400'"
C44COLOR="rgb '#8B008B'"
C55COLOR="rgb '#696969'"
C66COLOR="rgb '#008080'"

## C vs P
#'$dftdat/$elem1-$elem2/D03/cvsp.dat' u 1:((\$2+\$5+\$7)/3) w lp lt 2 pt 6 lc 1 notitle,\
#'$dftdat/$elem1-$elem2/D03/cvsp.dat' u 1:((\$3+\$4+\$6)/3) w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/D03/cvsp.dat' u 1:((\$8+\$9+\$10)/3) w lp lt 2 pt 6 lc 5 notitle
#'$dftdat/$elem1-$elem2/G1/cvsp.dat' u 1:((\$2+\$5+\$7)/3) w lp lt 2 pt 6 lc 1 notitle,\
#'$dftdat/$elem1-$elem2/G1/cvsp.dat' u 1:((\$3+\$4+\$6)/3) w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/G1/cvsp.dat' u 1:((\$8+\$9+\$10)/3) w lp lt 2 pt 6 lc 5 notitle
#'$dftdat/$elem2/bcc/cvsp.dat' u 1:((\$2+\$5+\$7)/3) w lp lt 2 pt 6 lc 1 notitle,\
#'$dftdat/$elem2/bcc/cvsp.dat' u 1:((\$3+\$4+\$6)/3) w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem2/bcc/cvsp.dat' u 1:((\$8+\$9+\$10)/3) w lp lt 2 pt 6 lc 5 notitle
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:((\$2+\$5)/2) w lp lt 2 pt 6 lc 1 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:(\$3) w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:((\$4+\$6)/2) w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:(\$7) w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:((\$8+\$9+\$10)/3) w lp lt 2 pt 6 lc 5 notitle
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:2 w lp lt 2 pt 6 lc 1 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:3 w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:4 w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:5 w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:6 w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:7 w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:8 w lp lt 2 pt 6 lc 2 notitle,\
#'$dftdat/$elem1-$elem2/omega/cvsp.dat' u 1:9 w lp lt 2 pt 6 lc 5 notitle

#elastic constants versus pressure
echo "set terminal postscript enhanced color
set encoding iso_8859_1
set xlabel 'Pressure [GPa]'
set ylabel 'Elastic Constant [GPa]'
set key top left

set output '$dir/tests/D03_CVSP.ps'
set xrange [0:$PRESMAX]
set title 'D0_{3} C_{ij} vs P'
plot '$dir/tests/D03/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/D03/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/D03/C_VS_P.dat' u 1:4 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}'

set output '$dir/tests/G1_CVSP.ps'
set title 'G1 C_{ij} vs P'
plot '$dir/tests/G1/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/G1/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/G1/C_VS_P.dat' u 1:4 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}'

set output '$dir/tests/BCCnb_CVSP.ps'
set title 'bcc Nb C_{ij} vs P'
plot '$dir/tests/BCCnb/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/BCCnb/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/BCCnb/C_VS_P.dat' u 1:4 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}'
'$dftdat/$elem2/bcc/Nb_C11vsP.dat' u 1:2 w lp lw 2 lt 2 pt 6 lc $C11COLOR notitle,\
'$dftdat/$elem2/bcc/Nb_C12vsP.dat' u 1:2 w lp lw 2 lt 2 pt 6 lc $C12COLOR notitle,\
'$dftdat/$elem2/bcc/Nb_C44vsP.dat' u 1:2 w lp lw 2 lt 2 pt 6 lc $C44COLOR notitle

set output '$dir/tests/BCCti_CVSP.ps'
set title 'bcc Ti C_{ij} vs P'
plot '$dir/tests/BCCti/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/BCCti/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/BCCti/C_VS_P.dat' u 1:4 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}',\
'$dftdat/$elem1/bcc/cvsp.dat' u 1:((\$2+\$5+\$7)/3) w lp lt 2 pt 6 lc $C11COLOR notitle,\
'$dftdat/$elem1/bcc/cvsp.dat' u 1:((\$3+\$4+\$6)/3) w lp lt 2 pt 6 lc $C12COLOR notitle,\
'$dftdat/$elem1/bcc/cvsp.dat' u 1:((\$8+\$9+\$10)/3) w lp lt 2 pt 6 lc $C44COLOR notitle

set xrange [0:10]
set output '$dir/tests/HCPti_CVSP.ps'
set title 'hcp Ti C_{ij} vs P'
plot '$dir/tests/HCPti/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/HCPti/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/HCPti/C_VS_P.dat' u 1:4 w l lc $C13COLOR lw 2 lt 1 title 'C_{13}',\
'$dir/tests/HCPti/C_VS_P.dat' u 1:5 w l lc $C33COLOR lw 2 lt 1 title 'C_{33}',\
'$dir/tests/HCPti/C_VS_P.dat' u 1:6 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}',\
'$dftdat/$elem1/hcp/cvsp.dat' u 1:((\$2+\$5)/2) w lp lt 2 pt 6 lc $C11COLOR notitle,\
'$dftdat/$elem1/hcp/cvsp.dat' u 1:(\$3) w lp lt 2 pt 6 lc $C12COLOR notitle,\
'$dftdat/$elem1/hcp/cvsp.dat' u 1:((\$4+\$6)/2) w lp lt 2 pt 6 lc $C13COLOR notitle,\
'$dftdat/$elem1/hcp/cvsp.dat' u 1:(\$7) w lp lt 2 pt 6 lc $C33COLOR notitle,\
'$dftdat/$elem1/hcp/cvsp.dat' u 1:((\$8+\$9+\$10)/3) w lp lt 2 pt 6 lc $C44COLOR notitle

set xrange [0:$PRESMAX]
set output '$dir/tests/OMGti_CVSP.ps'
set title '{/Symbol w}-Ti C_{ij} vs P'
plot '$dir/tests/OMGti/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/OMGti/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/OMGti/C_VS_P.dat' u 1:4 w l lc $C13COLOR lw 2 lt 1 title 'C_{13}',\
'$dir/tests/OMGti/C_VS_P.dat' u 1:5 w l lc $C33COLOR lw 2 lt 1 title 'C_{33}',\
'$dir/tests/OMGti/C_VS_P.dat' u 1:6 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}',\
'$dftdat/$elem1/omega/cvsp.dat' u 1:((\$2+\$5)/2) w lp lt 2 pt 6 lc $C11COLOR notitle,\
'$dftdat/$elem1/omega/cvsp.dat' u 1:(\$3) w lp lt 2 pt 6 lc $C12COLOR notitle,\
'$dftdat/$elem1/omega/cvsp.dat' u 1:((\$4+\$6)/2) w lp lt 2 pt 6 lc $C13COLOR notitle,\
'$dftdat/$elem1/omega/cvsp.dat' u 1:(\$7) w lp lt 2 pt 6 lc $C33COLOR notitle,\
'$dftdat/$elem1/omega/cvsp.dat' u 1:((\$8+\$9+\$10)/3) w lp lt 2 pt 6 lc $C44COLOR notitle

set output '$dir/tests/OMG_CVSP.ps'
set title '{/Symbol w}-Ti_{3}Nb C_{ij} vs P'
plot '$dir/tests/omega/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/omega/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/omega/C_VS_P.dat' u 1:4 w l lc $C13COLOR lw 2 lt 1 title 'C_{13}',\
'$dir/tests/omega/C_VS_P.dat' u 1:5 w l lc $C33COLOR lw 2 lt 1 title 'C_{33}',\
'$dir/tests/omega/C_VS_P.dat' u 1:6 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}'

set output '$dir/tests/APP_CVSP.ps'
set title '{{/Symbol a}^{,,}}-Ti_{3}Nb C_{ij} vs P'
plot '$dir/tests/APP/C_VS_P.dat' u 1:2 w l lc $C11COLOR lw 2 lt 1 title 'C_{11}',\
'$dir/tests/APP/C_VS_P.dat' u 1:3 w l lc $C12COLOR lw 2 lt 1 title 'C_{12}',\
'$dir/tests/APP/C_VS_P.dat' u 1:4 w l lc $C13COLOR lw 2 lt 1 title 'C_{13}',\
'$dir/tests/APP/C_VS_P.dat' u 1:5 w l lc $C22COLOR lw 2 lt 1 title 'C_{22}',\
'$dir/tests/APP/C_VS_P.dat' u 1:6 w l lc $C23COLOR lw 2 lt 1 title 'C_{23}',\
'$dir/tests/APP/C_VS_P.dat' u 1:7 w l lc $C33COLOR lw 2 lt 1 title 'C_{33}',\
'$dir/tests/APP/C_VS_P.dat' u 1:8 w l lc $C44COLOR lw 2 lt 1 title 'C_{44}',\
'$dir/tests/APP/C_VS_P.dat' u 1:9 w l lc $C55COLOR lw 2 lt 1 title 'C_{55}',\
'$dir/tests/APP/C_VS_P.dat' u 1:10 w l lc $C66COLOR lw 2 lt 1 title 'C_{66}'
" | gnuplot

$ps2png $dir/tests/D03_CVSP.ps $dir/tests/D03_CVSP.png
$ps2png $dir/tests/G1_CVSP.ps $dir/tests/G1_CVSP.png
$ps2png $dir/tests/APP_CVSP.ps $dir/tests/APP_CVSP.png
$ps2png $dir/tests/OMG_CVSP.ps $dir/tests/OMG_CVSP.png
$ps2png $dir/tests/HCPti_CVSP.ps $dir/tests/HCPti_CVSP.png
$ps2png $dir/tests/OMGti_CVSP.ps $dir/tests/OMGti_CVSP.png
$ps2png $dir/tests/BCCti_CVSP.ps $dir/tests/BCCti_CVSP.png
$ps2png $dir/tests/BCCnb_CVSP.ps $dir/tests/BCCnb_CVSP.png

# determine errors for elastic constants etc
D03C11err=`python -c "print '%3.1f' % (100*($D03C11-$C11D03d)/$C11D03d)"`
D03C12err=`python -c "print '%3.1f' % (100*($D03C12-$C12D03d)/$C12D03d)"`
D03C44err=`python -c "print '%3.1f' % (100*($D03C44-$C44D03d)/$C44D03d)"`

G1C11err=`python -c "print '%3.1f' % (100*($G1C11-$C11G1d)/$C11G1d)"`
G1C12err=`python -c "print '%3.1f' % (100*($G1C12-$C12G1d)/$C12G1d)"`
G1C44err=`python -c "print '%3.1f' % (100*($G1C44-$C44G1d)/$C44G1d)"`

BCCnbC11err=`python -c "print '%3.1f' % (100*($BCCnbC11-$C11BCCnbd)/$C11BCCnbd)"`
BCCnbC12err=`python -c "print '%3.1f' % (100*($BCCnbC12-$C12BCCnbd)/$C12BCCnbd)"`
BCCnbC44err=`python -c "print '%3.1f' % (100*($BCCnbC44-$C44BCCnbd)/$C44BCCnbd)"`

BCCtiC11err=`python -c "print '%3.1f' % (100*($BCCtiC11-$C11BCCtid)/$C11BCCtid)"`
BCCtiC12err=`python -c "print '%3.1f' % (100*($BCCtiC12-$C12BCCtid)/$C12BCCtid)"`
BCCtiC44err=`python -c "print '%3.1f' % (100*($BCCtiC44-$C44BCCtid)/$C44BCCtid)"`

B2C11err=`python -c "print '%3.1f' % (100*($B2C11-$C11B2d)/$C11B2d)"`
B2C12err=`python -c "print '%3.1f' % (100*($B2C12-$C12B2d)/$C12B2d)"`
B2C44err=`python -c "print '%3.1f' % (100*($B2C44-$C44B2d)/$C44B2d)"`

HCPtiC11err=`python -c "print '%3.1f' % (100*($HCPtiC11-$C11HCPtid)/$C11HCPtid)"`
HCPtiC12err=`python -c "print '%3.1f' % (100*($HCPtiC12-$C12HCPtid)/$C12HCPtid)"`
HCPtiC13err=`python -c "print '%3.1f' % (100*($HCPtiC13-$C13HCPtid)/$C13HCPtid)"`
HCPtiC33err=`python -c "print '%3.1f' % (100*($HCPtiC33-$C33HCPtid)/$C33HCPtid)"`
HCPtiC44err=`python -c "print '%3.1f' % (100*($HCPtiC44-$C44HCPtid)/$C44HCPtid)"`

OMGtiC11err=`python -c "print '%3.1f' % (100*($OMGtiC11-$C11OMGtid)/$C11OMGtid)"`
OMGtiC12err=`python -c "print '%3.1f' % (100*($OMGtiC12-$C12OMGtid)/$C12OMGtid)"`
OMGtiC13err=`python -c "print '%3.1f' % (100*($OMGtiC13-$C13OMGtid)/$C13OMGtid)"`
OMGtiC33err=`python -c "print '%3.1f' % (100*($OMGtiC33-$C33OMGtid)/$C33OMGtid)"`
OMGtiC44err=`python -c "print '%3.1f' % (100*($OMGtiC44-$C44OMGtid)/$C44OMGtid)"`

omgC11err=`python -c "print '%3.1f' % (100*($omgC11-$C11omgd)/$C11omgd)"`
omgC12err=`python -c "print '%3.1f' % (100*($omgC12-$C12omgd)/$C12omgd)"`
omgC13err=`python -c "print '%3.1f' % (100*($omgC13-$C13omgd)/$C13omgd)"`
omgC33err=`python -c "print '%3.1f' % (100*($omgC33-$C33omgd)/$C33omgd)"`
omgC44err=`python -c "print '%3.1f' % (100*($omgC44-$C44omgd)/$C44omgd)"`
omgC66err=`python -c "print '%3.1f' % (100*($omgC66-$C66omgd)/$C66omgd)"`

APPC11err=`python -c "print '%3.1f' % (100*($APPC11-$C11APPd)/$C11APPd)"`
APPC12err=`python -c "print '%3.1f' % (100*($APPC12-$C12APPd)/$C12APPd)"`
APPC13err=`python -c "print '%3.1f' % (100*($APPC13-$C13APPd)/$C13APPd)"`
APPC22err=`python -c "print '%3.1f' % (100*($APPC22-$C22APPd)/$C22APPd)"`
APPC23err=`python -c "print '%3.1f' % (100*($APPC23-$C23APPd)/$C23APPd)"`
APPC33err=`python -c "print '%3.1f' % (100*($APPC33-$C33APPd)/$C33APPd)"`
APPC44err=`python -c "print '%3.1f' % (100*($APPC44-$C44APPd)/$C44APPd)"`
APPC55err=`python -c "print '%3.1f' % (100*($APPC55-$C55APPd)/$C55APPd)"`
APPC66err=`python -c "print '%3.1f' % (100*($APPC66-$C66APPd)/$C66APPd)"`

A15TiNb3bulkerr=`python -c "print '%3.1f' % (100*($A15TiNb3bulkp-$A15TiNb3bulkd)/$A15TiNb3bulkd)"`
A15TiNb3laterr=`python -c "print '%3.1f' % (100*($A15TiNb3eqp-$A15TiNb3lat)/$A15TiNb3lat)"`
A15TiNb3C11err=`python -c "print '%3.1f' % (100*($A15TiNb3C11-$C11A15TiNb3d)/$C11A15TiNb3d)"`
A15TiNb3C12err=`python -c "print '%3.1f' % (100*($A15TiNb3C12-$C12A15TiNb3d)/$C12A15TiNb3d)"`
A15TiNb3C44err=`python -c "print '%3.1f' % (100*($A15TiNb3C44-$C44A15TiNb3d)/$C44A15TiNb3d)"`

HCPtilaterr=`python -c "print '%3.1f' % (100*(${HCPtiPLAT[0]}-$hcpTilat)/$hcpTilat)"`
HCPticoaerr=`python -c "print '%3.1f' % (100*($HCPticoap-$hcpTicoa)/$hcpTicoa)"`
OMGtilaterr=`python -c "print '%3.1f' % (100*(${OMGtiPLAT[0]}-$omgTilat)/$omgTilat)"`
OMGticoaerr=`python -c "print '%3.1f' % (100*($OMGticoap-$omgTicoa)/$omgTicoa)"`
BCCtilaterr=`python -c "print '%3.1f' % (100*(${BCCtiPLAT[0]}-$bccTilat)/$bccTilat)"`
BCCnblaterr=`python -c "print '%3.1f' % (100*(${BCCnbPLAT[0]}-$bccTilat)/$bccTilat)"`

HCPtieqp=`printf '%.3f' ${HCPtiPLAT[0]}`
OMGtieqp=`printf '%.3f' ${OMGtiPLAT[0]}`
BCCtieqp=`printf '%.3f' ${BCCtiPLAT[0]}`
BCCnbeqp=`printf '%.3f' ${BCCnbPLAT[0]}`

HCPticoap=`printf '%.3f' $HCPticoap`
OMGticoap=`printf '%.3f' $OMGticoap`

## Anisotropy for G1
#echo "Plotting G1 Anisotropy..."
#echo "Scubic[C11_, C12_, C44_] :=    Inverse[{{C11, C12, C12, 0, 0, 0},
#    							{C12, C11, C12, 0, 0, 0},
#    							{C12, C12, C11, 0, 0, 0},
#    							{0, 0, 0, C44, 0, 0},
#    							{0, 0, 0, 0, C44, 0},
#    							{0, 0, 0, 0, 0, C44}}];
#Ycubic[C11_, C12_, 
#  C44_, th_, ph_] := (Scubic[C11, C12, C44][[1,1]] - 
#   2*(Scubic[C11, C12, C44][[1,1]] - 
#      Scubic[C11, C12, C44][[1,2]] - 
#      Scubic[C11, C12, C44][[4,4]]/
#       2)*(alph^2 * beta^2 + beta^2 * gamm^2 + alph^2 * gamm^2))^-1;
#Gcubic[C11_, C12_, 
#  C44_, th_, ph_] := (Scubic[C11, C12, C44][[4,4]] - 
#   4*(Scubic[C11, C12, C44][[1,1]] - 
#      Scubic[C11, C12, C44][[1,2]] - 
#      Scubic[C11, C12, C44][[4,4]]/
#       2)*(alph^2 * beta^2 + beta^2 * gamm^2 + alph^2i * gamm^2))^-1;
#
#h = Sin[th]*Cos[ph];
#k = Sin[th]*Sin[ph];
#l = Cos[th];
#alph = h/Sqrt[h^2 + k^2 + l^2];
#beta = k/Sqrt[h^2 + k^2 + l^2];
#gamm = l/Sqrt[h^2 + k^2 + l^2];
#
#C11g1  = $C11G1d; C12g1  = $C12G1d; C44g1  = $C44G1d;
#C11g1p = $G1C11; C12g1p = $G1C12; C44g1p = $G1C44;
#
#col1lab = 
# Style[\"DFT\", Black, Bold, 20, FontFamily -> \"Times\"]; col2lab = 
# Style[\"MEAM\", Black, Bold, 20, FontFamily -> \"Times\"];
#row1lab = 
# Rotate[Style[\"Young's Modulus E [GPa]\", Black, Bold, 20, 
#   FontFamily -> \"Times\"], Pi/2]; row2lab = 
# Rotate[Style[\"Shear Modulus G [GPa]\", Black, Bold, 20, 
#   FontFamily -> \"Times\"], Pi/2]; colorData = \"Rainbow\";
#
#axes[x_, y_, z_, f_, a_] := 
#  Show[Graphics3D[
#    Join[{Thick, {Arrowheads[a]}, 
#      Arrow[{{0, 0, 0}, #}] & /@ {{1*x, 0, 0}, {0, -1*y, 0}, {0, 0, 
#         1*z}}}, {Text[
#       Style[\"[010]\", Bold, FontSize -> Scaled[f]], {1.11*x, 0, 0}], 
#      Text[Style[\"[100]\", Bold, FontSize -> Scaled[f]], {0, -1.1*y, 
#        0}], Text[
#       Style[\"[001]\", Bold, FontSize -> Scaled[f]], {0, 0, 
#        1.1*z}]}]], Boxed -> False, 
#   Method -> {\"ShrinkWrap\" -> False}];
#MYcubic[c11_, c12_, c44_] := 
#  Max[{Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/2},    						
#    Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/2},    							
#    Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/4},    							
#    Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/4},    							
#    Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> 0},   							
#    Ycubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0/4},    							
#    Ycubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0}}];
#
#mYcubic[c11_, c12_, c44_] := 
# Min[{Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/2},   							
#   Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/2},   							
#   Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/4},   							
#   Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/4},   							
#   Ycubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> 0},  							
#   Ycubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0/4},   							
#   Ycubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0}}]; i
#
#MGcubic[c11_, c12_, c44_] := 
# Max[{Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/2},   							
#   Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/2},   							
#   Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/4},   							
#   Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/4},   							
#   Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> 0},  							
#   Gcubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0/4},   							
#   Gcubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0}}];
#
#mGcubic[c11_, c12_, c44_] := 
#  Min[{Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/2},    							
#    Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/2},    							
#    Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/2, ph -> Pi/4},    							
#    Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> Pi/4},    							
#    Gcubic[c11, c12, c44, th, ph] /. {th -> Pi/4, ph -> 0}, 							
#    Gcubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0/4},    							
#    Gcubic[c11, c12, c44, th, ph] /. {th -> 0, ph -> 0}}];
#
#Ycplot[c11_, c12_, c44_] := 
# Legended[Show[
#   axes[0.7*MYcubic[c11, c12, c44], 0.7*MYcubic[c11, c12, c44], 
#    0.7*MYcubic[c11, c12, c44], 0.05, 0.025], 
#   SphericalPlot3D[
#    Ycubic[c11, c12, c44, th, ph], {th, 0, Pi}, {ph, 0, 2*Pi}, 
#    ColorFunction -> (ColorData[colorData][#6] &), 
#    ImageSize -> {300, 300}, PlotPoints -> 50, Boxed -> False, 
#    Mesh -> 25, Axes -> None, Lighting -> \"Neutral\"], 
#   Lighting -> \"Neutral\", ImageSize -> {300, 300}, 
#   ViewAngle -> 27*Degree], 
#  BarLegend[{colorData, {mYcubic[c11, c12, c44], MYcubic[c11, c12, c44]}}]]; 
#
#Gcplot[c11_, c12_, c44_] := 
# Legended[Show[
#   axes[0.7*MGcubic[c11, c12, c44], 0.7*MGcubic[c11, c12, c44], 
#    0.7*MGcubic[c11, c12, c44], 0.05, 0.025], 
#   SphericalPlot3D[
#    Gcubic[c11, c12, c44, th, ph], {th, 
#     0, Pi}, {ph, 0, 2*Pi}, 
#    ColorFunction -> (ColorData[colorData][#6] &), 
#    ImageSize -> {300, 300}, PlotPoints -> 50, Boxed -> False, 
#    Mesh -> 25, Axes -> None, Lighting -> \"Neutral\"], 
#   Lighting -> \"Neutral\", ImageSize -> {300, 300}, 
#   ViewAngle -> 27*Degree],
#  BarLegend[{colorData, {mGcubic[c11, c12, c44], MGcubic[c11, c12, c44]}}]];
#
#Export[\"$dir/tests/G1_anisotropy.png\",
# TableForm[{{Ycplot[C11g1, C12g1, C44g1],
#    Ycplot[C11g1p, C12g1p, C44g1p]}, {Gcplot[C11g1, C12g1, C44g1],
#    Gcplot[C11g1p, C12g1p, C44g1p]}}, 
#  TableHeadings -> {{row1lab, row2lab}, {col1lab, col2lab}}, 
#  TableAlignments -> Center], ImageResolution -> 200] " | math -noprompt

########################################################################################
#	MIXING ENTHALPY PLOTS
########################################################################################

echo "set terminal postscript enhanced color
set encoding iso_8859_1
set output '$dir/tests/THERMO.ps'

set title 'Mixing Enthalpies at P = 0 atm'
set xlabel 'x_{Nb} [at. %]'
set ylabel '{/Symbol D}H_{mix} [eV/atom]'

set key top center

plot '$dir/tests/thermo/BCC_formation.dat' u 1:2:3 with yerrorbars lt 1 pt 1 ps 0 lc rgb 'blue' notitle,\
'$dir/tests/thermo/BCC_formation.dat' u 1:2 smooth csplines lt 1 lw 2 lc rgb 'blue' title 'BCC',\
'$dir/tests/thermo/HCP_formation.dat' u 1:2:3 with yerrorbars lt 1 pt 2 ps 0 lc rgb 'red' notitle,\
'$dir/tests/thermo/HCP_formation.dat' u 1:2 smooth csplines lt 1 lw 2 lc rgb 'red' title 'HCP' " | gnuplot

$ps2png $dir/tests/THERMO.ps $dir/tests/THERMO.png

########################################################################################
#	PHASE TRANSITION PLOTS
########################################################################################

COLOR1="rgb '#FFFFFF'"
COLOR2="rgb '#0000FF'"
COLOR3="rgb '#FF8000'"
COLOR4="rgb '#FF0000'"
COLOR5="rgb '#00CC00'"
COLOR6="rgb '#990099'"


# alpha-prime-dprime
echo "set terminal postscript enhanced color
set encoding iso_8859_1
set output '$dir/tests/ALPHAS.ps'

set multiplot layout 1,2

#set title 'Solid Solution {/Symbol a} - {/Symbol a}^{,} - {/Symbol a}^{,,} Transition'
set key top left title 'Nb Fraction' 
set ylabel 'Change in formation energy [meV/atom]'
set xlabel 'Shift along [1~1{.8-}00] [{\305}]'
set xtics 0,0.1,0.5
#set mxtics 0.05,0.1,0.45
plot '$dir/tests/transitions/chi_e_0.0.dat' u 1:3 w l lt 1 lw 2 lc $COLOR1 title '0.0',\
'$dir/tests/transitions/chi_e_0.1.dat' u 1:3 w l lt 1 lw 2 lc $COLOR2 title '0.1',\
'$dir/tests/transitions/chi_e_0.2.dat' u 1:3 w l lt 1 lw 2 lc $COLOR3 title '0.2',\
'$dir/tests/transitions/chi_e_0.3.dat' u 1:3 w l lt 1 lw 2 lc $COLOR4 title '0.3',\
'$dir/tests/transitions/chi_e_0.4.dat' u 1:3 w l lt 1 lw 2 lc $COLOR7 title '0.4',\
'$dir/tests/transitions/chi_e_0.5.dat' u 1:3 w l lt 1 lw 2 lc $COLOR6 title '0.5'

set key top right # title ''
set ylabel '{/Symbol g} [{\260}]'
set y2label 'c/a ratio'
set y2tics 0.05
set ytics 1 nomirror
set xlabel 'Nb Fraction'
#set title 'Cell Parameters'

plot '$dir/tests/transitions/alpha_box.dat' u 1:3 w l lt 1 lw 2 lc rgb 'red' title '{/Symbol g}',\
'' u 1:2 w l axes x1y2 lt 1 lw 2 lc rgb 'blue' title 'c/a'" | gnuplot

$ps2png $dir/tests/ALPHAS.ps $dir/tests/ALPHAS.png

#plot beta-omega
echo "set terminal postscript enhanced color
set encoding iso_8859_1
set output '$dir/tests/BETA_TO_OMEGA.ps'

set xlabel 'Reaction Coordinate {/Symbol x}'
set ylabel 'Energy Barrier {/Symbol D}E [meV/atom]'
set title '{/Symbol b}_{G1} to {/Symbol w} Transition via <111> Planar Collapse'

set label 1 '{/Symbol w}' at 0.0,-2 
set label 1 '{/Symbol b}_{G1}' at 1.0,-2 

set key top left

plot '$dir/tests/transitions/beta_omega_0.dat' u 1:2 w l lc rgb 'blue' lw 3 title '$typ',\
'$dftdat/$elem1-$elem2/transitions/beta_omega.dat' u 1:3 w p lc rgb 'blue' pt 8 title 'PAW-PBE',\
'$dftdat/$elem1-$elem2/transitions/beta_omega.dat' u 1:3 smooth csplines lc rgb 'blue' lt 2 lw 1 notitle
#'$dir/tests/transitions/beta_omega_25.dat' u 1:2 w l lc rgb 'orange' lw 3 title '25 GPa',\
#'$dir/tests/transitions/beta_omega_50.dat' u 1:2 w l lc rgb 'purple' lw 3 title '50 GPa',\
#'$dir/tests/transitions/beta_omega_75.dat' u 1:2 w l lc rgb 'green' lw 3 title '75 GPa',\
#'$dir/tests/transitions/beta_omega_100.dat' u 1:2 w l lc rgb 'red' lw 3 title '100 GPa'" | gnuplot

$ps2png $dir/tests/BETA_TO_OMEGA.ps $dir/tests/BETA_TO_OMEGA.png

#plot beta-omega in pure Ti
echo "set terminal postscript enhanced color
set encoding iso_8859_1
set output '$dir/tests/BETA_TO_OMEGA_Ti.ps'

set xlabel 'Reaction Coordinate {/Symbol x}'
set ylabel 'Energy Barrier {/Symbol D}E [meV/atom]'
set title '{/Symbol b} to {/Symbol w} Transition via <111> Planar Collapse in pure Ti'

set label 1 '{/Symbol w}' at 0.0,-2 
set label 1 '{/Symbol b}' at 1.0,-2 

plot \
'$dir/tests/transitions/betaOmegaTi0.dat' u 1:3 w l lc rgb 'red' lw 3 title '0 GPa',\
'$dftdat/$elem1/transitions/BETA-OMEGA.dat' u 1:2 w p lc rgb 'red' pt 8 title 'PAW-PBE',\
'$dftdat/$elem1/transitions/BETA-OMEGA.dat' u 1:2 smooth csplines lc rgb 'red' lt 2 lw 1 notitle,\
'$dir/tests/transitions/betaOmegaTi25.dat' u 1:3 w l lc rgb 'blue' lw 3 title '25 GPa',\
'$dftdat/$elem1/transitions/BETA-OMEGA-25.dat' u 1:3 w p lc rgb 'blue' pt 8 title 'PAW-PBE',\
'$dftdat/$elem1/transitions/BETA-OMEGA-25.dat' u 1:3 smooth csplines lc rgb 'blue' lt 2 lw 1 notitle,\
'$dir/tests/transitions/betaOmegaTi50.dat' u 1:3 w l lc rgb 'green' lw 3 title '50 GPa',\
'$dftdat/$elem1/transitions/BETA-OMEGA-50.dat' u 1:3 w p lc rgb 'green' pt 8 title 'PAW-PBE',\
'$dftdat/$elem1/transitions/BETA-OMEGA-50.dat' u 1:3 smooth csplines lc rgb 'green' lt 2 lw 1 notitle,\
'$dir/tests/transitions/betaOmegaTi75.dat' u 1:3 w l lc rgb 'orange' lw 3 title '75 GPa',\
'$dftdat/$elem1/transitions/BETA-OMEGA-75.dat' u 1:3 w p lc rgb 'orange' pt 8 title 'PAW-PBE',\
'$dftdat/$elem1/transitions/BETA-OMEGA-75.dat' u 1:3 smooth csplines lc rgb 'orange' lt 2 lw 1 notitle,\
'$dir/tests/transitions/betaOmegaTi100.dat' u 1:3 w l lc rgb 'purple' lw 3 title '100 GPa',\
'$dftdat/$elem1/transitions/BETA-OMEGA-100.dat' u 1:3 w p lc rgb 'purple' pt 8 title 'PAW-PBE',\
'$dftdat/$elem1/transitions/BETA-OMEGA-100.dat' u 1:3 smooth csplines lc rgb 'purple' lt 2 lw 1 notitle" | gnuplot

$ps2png $dir/tests/BETA_TO_OMEGA_Ti.ps $dir/tests/BETA_TO_OMEGA_Ti.png

mmzay="/n/jww-1/ehemann.2/pots/$typ/$elem1-$elem2/set.$setnum/pot.$potnum/tests/transitions/beta_app.dat"
dftay="/n/jww-1/ehemann.2/testingScripts/DFT_DATA/$elem1-$elem2/transitions/burgers.dat"

mmzti="/n/jww-1/ehemann.2/pots/$typ/$elem1-$elem2/set.$setnum/pot.$potnum/tests/transitions/alphaBetaTi.dat"
dftti="/n/jww-1/ehemann.2/testingScripts/DFT_DATA/Ti/transitions/burgers.dat"

meamin=`cat $mmzay | awk 'NR>1{print $3}' | sort -g | head -1` 
dftmin=`cat $dftay | awk 'NR>1{print 1000*$4}' | sort -g | head -1` 
meamax=`cat $mmzay | awk 'NR>1{print $3}' | sort -g | tail -1` 
dftmax=`cat $dftay | awk 'NR>1{print 1000*$4}' | sort -g | tail -1` 
dftminay=$dftmin
echo $dftmin $dftmax
echo $meamin $meamax
minay=`echo -e "$meamin \n $dftmin" | sort -g | head -1`
maxay=`echo -e "$meamax \n $dftmax" | sort -g | tail -1`

meamin=`cat $mmzti | awk 'NR>1{print $3}' | sort -g | head -1` 
dftmin=`cat $dftti | awk 'NR>1{print 1000*$4}' | sort -g | head -1` 
meamax=`cat $mmzti | awk 'NR>1{print $3}' | sort -g | tail -1` 
dftmax=`cat $dftti | awk 'NR>1{print 1000*$4}' | sort -g | tail -1` 
dftminti=$dftmin
echo $dftmin $dftmax
echo $meamin $meamax
minti=`echo -e "$meamin \n $dftmin" | sort -g | head -1`
maxti=`echo -e "$meamax \n $dftmax" | sort -g | tail -1`

echo "
set terminal postscript enhanced color size 7,7
set encoding iso_8859_1
set output '$dir/tests/ALPHA_BETA_BURGERS.ps'

set pm3d
set view map
set tmargin 0.
set bmargin 0.
set lmargin 0.
set rmargin 0.

set multiplot layout 2,2

set xrange [0.000:1.000]
set yrange [0.000:1.000]
#set size 0.50,0.44
#set size 0.60,0.60
#set size 0.75,0.66
#set size 1,1

set style line 1 lt 1 lc rgb 'black'
set style line 2 lt 1 lc rgb 'white'


######## meam alloy
f(u,v) = R*u**5 + S*(u**4)*v + T*(u**3)*(v**2) + U*(u**2)*(v**3) + V*u*(v**4) + W*v**5 + M*u**4 + N*u*v**3 + O*(u**2)*(v**2) + P*(u**3)*v + Q*v**4 + A*u**3 + B*(u**2)*v + C*u*(v**2) + D*v**3 + E*u*v + F*u**2 + G*v**2 + H*u + K*v + L
A=1.1250
B=-1.125
C=10000
D=-10000
E=100
F=-100
G=1000
H=-1000
I=100
J=-8
set isosamples 100
#splot f(x,y)
FIT_LIMIT=1e-8
fit [-0.1:1.1] [-0.1:1.1] f(x,y) '$mmzay' using 1:2:3:(1) via A,B,C,D,E,F,G,H,K,L,M,N,O,P,Q,R,S,T,U,V,W
set tmargin 0.
set bmargin 0.
set lmargin 0.
set rmargin 0.
set cbrange [$minay:$maxay]
set cblabel 'E-E_{{/Symbol a}^{,,}} [meV/atom]'
set colorbox vertical user origin 1.08,0.58115 size 0.075,0.4681
#set origin 0.185,0.5
set origin 0.0,0.5
set title '$typ' offset 0,-0.5
set ytics 0.1,0.2,0.9 out scale 0.5 nomirror
unset xlabel
unset ylabel
unset xtics
set label 111 '{/Symbol b}_{G1}' at 0.075,0.075 center textcolor rgb 'white'
set label 112 '{/Symbol a}^{,,}' at 0.925,0.925 center textcolor rgb 'white'
set label 33 'Ti_{3}Nb' at -0.25,0.5 center rotate by 90
set label 2 'Shuffle {/Symbol h}' at -0.183,0.0 center rotate by 90
#set size square
set size 0.70,0.70
set xrange [0:1]
set yrange [0:1]
#set pm3d
set contour
unset clabel
set cntrparam levels 25
set view map
splot f(x,y) w pm3d notitle

######## dft alloy
f(u,v) = R*u**5 + S*(u**4)*v + T*(u**3)*(v**2) + U*(u**2)*(v**3) + V*u*(v**4) + W*v**5 + M*u**4 + N*u*v**3 + O*(u**2)*(v**2) + P*(u**3)*v + Q*v**4 + A*u**3 + B*(u**2)*v + C*u*(v**2) + D*v**3 + E*u*v + F*u**2 + G*v**2 + H*u + K*v + L
A=1.1250
B=-1.125
C=10000
D=-10000
E=100
F=-100
G=1000
H=-1000
I=100
J=-8
set isosamples 100
#splot f(x,y)
FIT_LIMIT=1e-8
fit [-0.1:1.1] [-0.1:1.1] f(x,y) '$dftay' using 1:2:(1000*\$4):(1) via A,B,C,D,E,F,G,H,K,L,M,N,O,P,Q,R,S,T,U,V,W

set tmargin 0.
set bmargin 0.
set lmargin 0.
set rmargin 0.
unset label 2
unset label 33
unset label 111
unset label 112
unset ytics
set origin 0.5,0.5
set title 'PAW-PBE' offset 0,-0.5
set label 121 '{/Symbol b}_{G1}' at 0.075,0.075 center textcolor rgb 'white'
set label 122 '{/Symbol a}^{,,}' at 0.925,0.925 center textcolor rgb 'white'
show label

#set size square
set size 0.70,0.70
unset xlabel
unset ylabel
set cbrange [$minay:$maxay]
set xrange [0:1]
set yrange [0:1]
#set pm3d
set contour
unset clabel
set cntrparam levels 25
set view map
splot f(x,y) w pm3d notitle

######## meam pure ti
f(u,v) = R*u**5 + S*(u**4)*v + T*(u**3)*(v**2) + U*(u**2)*(v**3) + V*u*(v**4) + W*v**5 + M*u**4 + N*u*v**3 + O*(u**2)*(v**2) + P*(u**3)*v + Q*v**4 + A*u**3 + B*(u**2)*v + C*u*(v**2) + D*v**3 + E*u*v + F*u**2 + G*v**2 + H*u + K*v + L
A=1.1250
B=-1.125
C=10000
D=-10000
E=100
F=-100
G=1000
H=-1000
I=100
J=-8
set isosamples 100
#splot f(x,y)
FIT_LIMIT=1e-8
fit [-0.1:1.1] [-0.1:1.1] f(x,y) '$mmzti' using 1:2:3:(1) via A,B,C,D,E,F,G,H,K,L,M,N,O,P,Q,R,S,T,U,V,W

set tmargin 0.
set bmargin 0.
set lmargin 0.
set rmargin 0.
unset colorbox
unset label 121
unset label 122
set cbrange [$minti:$maxti]
set cblabel 'E-E_{{/Symbol a}} [meV/atom]'
#set colorbox vertical user origin 1.05,0.05 size 0.05,0.50
set colorbox vertical user origin 1.08,0.08500 size 0.075,0.4885
set origin 0.0,0.0
unset title
unset xlabel
unset ylabel
set xtics 0.1,0.2,0.9 out scale 0.5 nomirror
set ytics 0.1,0.2,0.9 out scale 0.5 nomirror
set label 211 '{/Symbol b}' at 0.075,0.075 center textcolor rgb 'white'
set label 212 '{/Symbol a}' at 0.925,0.925 center textcolor rgb 'white'
set label 1 'Shear {/Symbol x}' at 1.0,-0.155 center
set label 44 'Ti' at -0.25,0.5 center rotate by 90

#set size square
set size 0.70,0.70
set xrange [0:1]
set yrange [0:1]
#set pm3d
set contour
unset clabel
set cntrparam levels 25
set view map
splot f(x,y) w pm3d notitle

######## dft pure ti
f(u,v) = R*u**5 + S*(u**4)*v + T*(u**3)*(v**2) + U*(u**2)*(v**3) + V*u*(v**4) + W*v**5 + M*u**4 + N*u*v**3 + O*(u**2)*(v**2) + P*(u**3)*v + Q*v**4 + A*u**3 + B*(u**2)*v + C*u*(v**2) + D*v**3 + E*u*v + F*u**2 + G*v**2 + H*u + K*v + L
A=1.1250
B=-1.125
C=10000
D=-10000
E=100
F=-100
G=1000
H=-1000
I=100
J=-8
set isosamples 100

#splot f(x,y)
FIT_LIMIT=1e-8
fit [-0.1:1.1] [-0.1:1.1] f(x,y) '$dftti' using 1:2:(1000*\$4):(1) via A,B,C,D,E,F,G,H,K,L,M,N,O,P,Q,R,S,T,U,V,W

set tmargin 0.
set bmargin 0.
set lmargin 0.
set rmargin 0.
set cbrange [$minti:$maxti]
unset label 1
unset label 44
unset label 211
unset label 212
unset ytics
set origin 0.5,0.00
set label 221 '{/Symbol b}' at 0.075,0.075 center textcolor rgb 'white'
set label 222 '{/Symbol a}' at 0.925,0.925 center textcolor rgb 'white'
show label
unset xlabel
unset ylabel
set size 0.70,0.70
set xrange [0:1]
set yrange [0:1]
#set pm3d
set contour
unset clabel
set cntrparam levels 25
set view map
splot f(x,y) w pm3d notitle
" > $dir/tests/fitsurf.gp

echo "load '$dir/tests/fitsurf.gp'" | gnuplot

$fixbb $dir/tests/ALPHA_BETA_BURGERS.ps
$ps2png $dir/tests/ALPHA_BETA_BURGERS.ps $dir/tests/ALPHA_BETA_BURGERS.png

# plot BCT deformations
echo "set terminal postscript enhanced color
set encoding iso_8859_1
set output '$dir/tests/BCT_DEFORMATION.ps'

set title 'BCT deformation for {/Symbol b} structures'
set xlabel 'c/a ratio'
set ylabel '{/Symbol D}E [meV/atom]'

set xrange [0.75:1.65]
set yrange [-50:1800]
set xtics 0,1,0
set xtics add (\"0.8\" 0.8)
set xtics add (\"BCC\" 1.0)
set xtics add (\"1.2\" 1.2)
set xtics add (\"FCC\" 1.41421356)
set xtics add (\"1.6\" 1.6)
set ytics 0,200,1800
set key top center

set arrow 1 from 0.0,1800 to 0.0,-50 nohead lt 1 lw 2 lc 0
set arrow 2 from 1.41421356,1800 to 1.41421356,-50 nohead lt 1 lw 2 lc 0

show arrow 1
show arrow 2

plot '$dir/tests/BCCti/bct.dat' u 1:(1000*\$2) w l lw 3 lt 1 lc rgb '#1E90FF' title 'Ti',\
'$dir/tests/D03/bct.dat' u 1:(1000*\$2) w l lw 3 lt 2 lc rgb 'red' title 'D0_{3} Ti_{3}Nb',\
'$dir/tests/G1/bct.dat' u 1:(1000*\$2) w l lw 3 lt 3 lc rgb 'dark-green' title 'G1 Ti_{3}Nb',\
'$dir/tests/B2/bct.dat' u 1:(1000*\$2) w l lw 3 lt 1 lc rgb 'dark-orange' title 'B2 TiNb',\
'$dir/tests/BCCnb/bct.dat' u 1:(1000*\$2) w l lw 3 lt 1 lc rgb 'dark-violet' title 'Nb',\
'$dftdat/$elem1/bcc/tetra.dat' u 1:3 w lp lw 1 lt 2 lc rgb '#1E90FF' notitle,\
'$dftdat/$elem1-$elem2/D03/tetra.dat' u 1:3 w lp lw 1 lt 2 lc rgb 'red' notitle,\
'$dftdat/$elem1-$elem2/TiNb/B2/tetra.dat' u 1:3 w lp lw 1 lt 2 lc rgb 'dark-orange' notitle,\
'$dftdat/$elem2/bcc/tetra.dat' u 1:3 w lp lw 1 lt 2 lc rgb 'dark-violet' notitle
" | gnuplot

$ps2png $dir/tests/BCT_DEFORMATION.ps $dir/tests/BCT_DEFORMATION.png

########################################################################################
#  stsacking fault plots
########################################################################################

echo "set terminal postscript enhanced color
set encoding iso_8859_1
set xrange [0:1]
set key top right

set output '$dir/tests/SF_BCC_Nb.ps'
set title 'BCC Nb Stacking Faults'
set xlabel 'Fractional Displacement Along <111>'
set ylabel '{/Symbol g} [meV/{\305}^{2}]'
set key bottom center
plot '$dir/tests/BCCnb/0SF_110_111.dat' u 1:2 w l lc rgb 'dark-orange' lt 1 lw 3 title '\\{110\\} - unrelaxed',\
'$dir/tests/BCCnb/0SF_112_111.dat' u 1:2 w l lc rgb 'dark-violet' lt 1 lw 3 title '\\{112\\} - unrelaxed',\
'$dir/tests/BCCnb/0SF_110_111_relaxed.dat' u 1:2 w l lc rgb 'dark-orange' lt 0 lw 3 title '\\{110\\} - relaxed',\
'$dir/tests/BCCnb/0SF_112_111_relaxed.dat' u 1:2 w l lc rgb 'dark-violet' lt 0 lw 3 title '\\{112\\} - relaxed',\
'$dftdat/Nb/bcc/gamma110_111.dat' u 1:2 w lp pt 6 lt 2 lc rgb 'dark-orange' notitle,\
'$dftdat/Nb/bcc/gamma112_111.dat' u 1:2 w lp pt 6 lt 2 lc rgb 'dark-violet' notitle

set output '$dir/tests/SF_HCP_011b0_112b0.ps
set title 'HCP Ti \\{01~1{.8-}0\\}<11~2{.8-}0> Stacking Fault (prismatic)'
set xlabel 'Fractional Displacement Along <11~2{.8-}0>'
set ylabel '{/Symbol g} [meV/{\305}^{2}]'
set key bottom center
plot '$dir/tests/HCPti/0SF_011b0e_112b0.dat' u 1:2 w l lc rgb 'red' lt 1 lw 3 title 'Easy - unrelaxed',\
'$dir/tests/HCPti/0SF_011b0h_112b0.dat' u 1:2 w l lc rgb 'blue' lt 1 lw 3 title 'Hard - unrelaxed',\
'$dir/tests/HCPti/0SF_011b0e_112b0_relaxed.dat' u 1:2 w l lc rgb 'red' lt 0 lw 3 title 'Easy - relaxed',\
'$dir/tests/HCPti/0SF_011b0h_112b0_relaxed.dat' u 1:2 w l lc rgb 'blue' lt 0 lw 3 title 'Hard - relaxed',\
'$dftdat/$elem1/hcp/gamma01-10e11-20.dat' u 1:2 w lp pt 6 lt 2 lc rgb 'red' notitle,\
'$dftdat/$elem1/hcp/gamma01-10h11-20.dat' u 1:2 w lp pt 6 lt 2 lc rgb 'blue' notitle" | gnuplot

$ps2png $dir/tests/SF_HCP_011b0_112b0.ps $dir/tests/SF_HCP_011b0_112b0.png
$ps2png $dir/tests/SF_BCC_Nb.ps $dir/tests/SF_BCC_Nb.png




echo "FINISHED PLOTTING!"
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#	HTML 	GENERATION		PORTION
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
echo "copying files and generating .html..."

# copy lammps potential
cp $dir/lammps.pt ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/$elem1$elem2$typ$setnum$potnum.pt

# performance / potential / error plots
cp $dir/tests/radperf.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/radperf_small.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/ERRORS.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/SPLINES.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# E-V curves
cp $dir/tests/E_VS_V_BETA.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/E_VS_V_META.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/E_VS_V_Ti.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/E_VS_V_Nb.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/E_VS_V_5050.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/E_VS_V_2575.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/E_VS_V_3367.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/E_VS_V_6733.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# copy P-V curves
cp $dir/tests/D03_PVSV.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/G1_PVSV.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/SQS7525_PVSV.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/APP_PVSV.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/AP_PVSV.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/OMEGA_PVSV.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# pressure dependence
cp $dir/tests/APP_ABCVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/HCPti_ABCVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/OMGTi2Nb_ABCVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/OMGTiNb2_ABCVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# elastic constants
cp $dir/tests/D03_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/G1_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/APP_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/OMG_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/BCCnb_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/BCCti_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/HCPti_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/OMGti_CVSP.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# anisotropy
cp $dir/tests/G1_anisotropy.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# phase transitions
cp $dir/tests/BETA_TO_OMEGA.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/BETA_TO_OMEGA_Ti.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/BETA_TO_DPRIM.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/ALPHA_BETA_BURGERS.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/BCT_DEFORMATION.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/ALPHAS.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/


# stacking faults
cp $dir/tests/SF_HCP_011b0_112b0.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/SF_BCC_Nb.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# thermodynamics
cp $dir/tests/THERMO.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# phonons
cp $dir/tests/hcpTi_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/omgTi_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/bccTi_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/bccNb_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/A3_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/B2_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/G1_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/D03_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/APP_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/omega_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/
cp $dir/tests/A15TiNb3_phonons.png ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/

# now generate html file for results page
cat > ~/public_html/POTS/$elem1-$elem2/$typ/"set.$setnum"/"pot.$potnum"/results.html << !!
<!DOCTYPE html>
<html>

<head>
<link rel="stylesheet" type="text/css" href="../../../../../pot_style.css">
</head>

<body>

<div id="container">

<div id="buffer"></div>

<div id="header">
	<img src="../../../../../DFT/Ti/titanium.png" style="position:relative;width:auto;height:215px" align="left">
	<img src="../../../../../DFT/Ti3Nb/niobium.png" style="position:relative;width:auto;height:215px" align="right">
	<br><br>
	<h1><a href ="./$elem1$elem2$typ$setnum$potnum.pt"> $elem1-$elem2 $typ Potential $setnum.$potnum </a></h1>
</div>
<div id="nav">
	<ul>
		<li><a href="http://www.physics.ohio-state.edu/~ehemann.2"> Home </a></li>
		<li><a href="#"> Density Functional Theory </a>
			<ul>
				<li><a href="http://www.physics.ohio-state.edu/~ehemann.2/DFT/Ti/tidft.html"> Ti </a></li>
				<li><a href="http://www.physics.ohio-state.edu/~ehemann.2/DFT/Ti3Nb/ti3nbdft.html"> Ti<sub>3</sub>Nb </a></li>
			</ul>
		</li>
		<li><a href="#"> Interatomic Potentials </a>
			<ul>
				<li><a href="http://www.physics.ohio-state.edu/~ehemann.2/POTS/Ti/list.html"> Ti </a></li>
				<li><a href="http://www.physics.ohio-state.edu/~ehemann.2/POTS/Ti-Nb/list.html"> Ti-Nb </a></li>
			</ul>
		</li>
		<li><a href="#"> Molecular Dynamics </a>
			<ul>
				<li><a href="http://www.physics.ohio-state.edu/~ehemann.2/MD/Ti"> Ti </a></li>
			</ul>
		</li>
		</ul>
</div>
<div id="section">

	<h4> Overall Performance </h4>
	<figure align="center" style="width:75%">
		<img src="./radperf.png" style="width:100%">
	</figure>

	<h4> $typ Functions </h4>
	<figure align="center" style="width:85%">
		<img src="./SPLINES.png" style="width:100%">
	</figure>

	<h4> Fitting Errors: Sum total fitting error = $fitError </h4>
	<img src="./ERRORS.png" style="width:85%">

	<h4> Structural Properties for pure elements </h4>
	<table align="center" style="width:50%" border="1">
	<tr>
		<th>	</th><th>$typ</th><th>PAW-PBE</th><th>error (%)</th>
	</tr>
	<tr>
		<th> a<sub>bcc-Ti</sub> [&Aring;]</td><td> $BCCtieqp </td><td> $bccTilat </td><td> $BCCtilaterr </td>
	</tr>
	<tr>
		<th> a<sub>bcc-Nb</sub> [&Aring;]</td><td> $BCCnbeqp </td><td> $bccNblat </td><td> $BCCnblaterr </td>
	</tr>
	<tr>
		<th> a<sub>hcp-Ti</sub> [&Aring;]</td><td> $HCPtieqp </td><td> $hcpTilat </td><td> $HCPtilaterr </td>
	</tr>
	<tr>
		<th> (c/a)<sub>hcp-Ti</sub> [&Aring;]</td><td> $HCPticoap </td><td> $hcpTicoa </td><td> $HCPticoaerr </td>
	</tr>
	<tr>
		<th> a<sub>&omega;-Ti</sub> [&Aring;]</td><td> $OMGtieqp </td><td> $omgTilat </td><td> $OMGtilaterr </td>
	</tr>
	<tr>
		<th> (c/a)<sub>&omega;-Ti</sub> [&Aring;]</td><td> $OMGticoap </td><td> $omgTicoa </td><td> $OMGticoaerr </td>
	</tr>
	</table>
	</table>

	<h4> Structural Properties for &beta; Supercells </h4>
	<table align="center" style="width:50%" border="1">
	<tr>
		<th>	</th><th>$typ</th><th>PAW-PBE</th><th>error (%)</th>
	</tr>
	<tr>
		<th> a<sub>D0<sub>3</sub></sub> [&Aring;]</td><td> $D03eqp </td><td> $D03lat </td><td> $D03laterr </td>
	</tr>	
	<tr>
		<th> a<sub>G1</sub> [&Aring;]</td><td> $G1eqp </td><td> $G1lat </td><td> $G1laterr </td>
	</tr>	
	<tr>
		<th> a<sub>SQS</sub> [&Aring;]</td><td> $SQS7525eqp </td><td> $SQS7525lat </td><td> $SQS7525laterr </td>
	</tr>	
	<tr>
		<th> E<sub>D0<sub>3</sub></sub>-E<sub>G1</sub> [meV/atom]</td><td> $D03pd </td><td> $D03dd </td><td> $D03edifferr </td>
	</tr>	
	<tr>
		<th> E<sub>SQS</sub>-E<sub>G1</sub> [meV/atom]</td><td> $SQS7525pd </td><td> $SQS7525dd </td><td> $SQS2575edifferr </td>
	</tr>	
	<tr>
		<th> B<sub>D0<sub>3</sub></sub> [GPa]</td><td> $D03bulkp </td><td> $D03bulkd </td><td> $D03Berr </td>
	</tr>
	<tr>
		<th> B<sub>G1</sub> [GPa]</td><td> $G1bulkp </td><td> $G1bulkd </td><td> $G1Berr </td>
	</tr>
	<tr>
		<th> B<sub>SQS</sub> [GPa]</td><td> $SQS7525bulkp </td><td> $SQS7525bulkd </td><td> $SQS7525Berr </td>
	</tr>
	</table>	

	<h4> Structural Properties for Metastable Phases </h4>
	<table align="center" style="width:50%" border="1">
	<tr>
		<th>	</th><th>$typ</th><th>PAW-PBE</th><th>error (%)</th>
	</tr>
	<tr>
		<th> a<sub>&alpha;<sup>,,</sup></sub> [&Aring;]</td><td> $APPeqp </td><td> $APPlat </td><td> $APPlaterr </td>
	</tr>	
	<tr>
		<th> a<sub>&alpha;<sup>,</sup></sub> [&Aring;]</td><td> $APeqp </td><td> $APlat </td><td> $APlaterr </td>
	</tr>	
	<tr>
		<th> a<sub>&omega;</sub> [&Aring;]</td><td> $omgeqp </td><td> $omglat </td><td> $omglaterr </td>
	</tr>	
	<tr>
		<th> E<sub>&alpha;<sup>,,</sup></sub>-E<sub>G1</sub> [meV/atom]</td><td> $APPpd </td><td> $APPdd </td><td> $APPedifferr </td>
	</tr>	
	<tr>
		<th> E<sub>&alpha;<sup>,</sup></sub>-E<sub>G1</sub> [meV/atom]</td><td> $APpd </td><td> $APdd </td><td> $APedifferr </td>
	</tr>	
	<tr>
		<th> E<sub>&omega;</sub>-E<sub>G1</sub> [meV/atom]</td><td> $omgpd </td><td> $omgdd </td><td> $omgedifferr </td>
	</tr>	
	<tr>
		<th> B<sub>&alpha;<sup>,,</sup></sub> [GPa]</td><td> $APPbulkp </td><td> $APPbulkd </td><td> $APPBerr </td>
	</tr>
	<tr>
		<th> B<sub>&alpha;<sup>,</sup></sub> [GPa]</td><td> $APbulkp </td><td> $APbulkd </td><td> $APBerr </td>
	</tr>
	<tr>
		<th> B<sub>&omega;</sub> [GPa]</td><td> $omgbulkp </td><td> $omgbulkd </td><td> $omgBerr </td>
	</tr>
	</table>	

	<h4> Structural Properties for 50/50 Phases </h4>
	<table align="center" style="width:50%" border="1">
	<tr>
		<th>	</th><th>$typ</th><th>PAW-PBE</th><th>error (%)</th>
	</tr>
	<tr>
		<th> a<sub>B2</sub> [&Aring;]</td><td> $B2eqp </td><td> $B2lat </td><td> $B2laterr </td>
	</tr>	
	<tr>
		<th> E<sub>B2</sub>-E<sub>G1</sub> [meV/atom]</td><td> $B2pd </td><td> $B2dd </td><td> $B2edifferr </td>
	</tr>	
	<tr>
		<th> B<sub>B2</sub> [GPa]</td><td> $B2bulkp </td><td> $B2bulkd </td><td> $B2Berr </td>
	</tr>
	</table>

	<h4> A15 TiNb<sub>3</sub> </h4>
	
	<table align="center" style="width:50%" border="1">
	<tr>
		<th>	</th><th>$typ</th><th>PAW-PBE</th><th>error (%)</th>
	</tr> <tr>
		<th> a [&Aring;]</th><td> $A15TiNb3eqp </td><td> $A15TiNb3lat </td> <td> $A15TiNb3laterr </td>
	</tr> <tr>
		<th> B [GPa] </th><td> $A15TiNb3bulkp </td><td> $A15TiNb3bulkd </td><td> $A15TiNb3bulkerr </td>
	</tr> <tr>
		<th> C<sub>11</sub> [GPa] </th><td> $A15TiNb3C11p </td><td> $C11A15TiNb3d </td><td> $A15TiNb3C11err </td>
	</tr> <tr>
		<th> C<sub>12</sub> [GPa] </th><td> $A15TiNb3C12p </td><td> $C12A15TiNb3d </td><td> $A15TiNb3C12err </td>
	</tr> <tr>
		<th> C<sub>44</sub> [GPa] </th><td> $A15TiNb3C44p </td><td> $C44A15TiNb3d </td><td> $A15TiNb3C44err </td>
	</tr> </table>


	<h4> Energy-volume curves </h4>
	<table style="width:100%">
	</tr><tr>
		<td><img src="./E_VS_V_Ti.png" style="width:100%"></td>
		<td><img src="./E_VS_V_Nb.png" style="width:100%"></td>
	</tr><tr>	
		<td><img src="./E_VS_V_BETA.png" style="width:100%"></td>
		<td><img src="./E_VS_V_META.png" style="width:100%"></td>
	</tr><tr>
		<td><img src="./E_VS_V_5050.png" style="width:100%"></td>
		<td><img src="./E_VS_V_2575.png" style="width:100%"></td>
	</tr><tr>
		<td><img src="./E_VS_V_6733.png" style="width:100%"></td>
		<td><img src="./E_VS_V_3367.png" style="width:100%"></td>
	</tr></table>

	<h4> Pressure-volume curves </h4>

	<table style="width:100%">
	<tr>
		<td> <img src="./D03_PVSV.png" style="width:100%"> </td><td> <img src="./G1_PVSV.png" style="width:100%"> </td><td> <img src="./SQS7525_PVSV.png" style="width:100%"> </td>
	</tr>
	<tr>
		<td> <img src="./APP_PVSV.png" style="width:100%"> </td><td> <img src="./AP_PVSV.png" style="width:100%"> </td><td> <img src="./OMEGA_PVSV.png" style="width:100%"> </td>
	</tr>
	</table>

	<h4> High-Pressure Lattice Constants </h4>
	<table style="width:100%">
	<tr>
		<td> <img src="./APP_ABCVSP.png" style="width:100%" </td>
		<td> <img src="./HCPti_ABCVSP.png" style="width:100%" </td>
	</tr><tr>
		<td> <img src="./OMGTi2Nb_ABCVSP.png" style="width:100%" </td>
		<td> <img src="./OMGTiNb2_ABCVSP.png" style="width:100%" </td>
	</tr>
	</table>

	<h4> Elastic Properties </h4>
	<h5> Alloy Phases </h5>
	<table border="1" width="50%" align="center">
		<tr> <th></th><th>[GPa]</th><th> $typ </th><th> PAW-PBE </th><th> Error (%) </th> </tr>
		<tr></tr>
		<tr> <th rowspan="3"> D0<sub>3</sub> </td><td> C<sub>11</sub> </td><td> $D03C11 </td><td> $C11D03d </td><td> $D03C11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $D03C12 </td><td> $C12D03d </td><td> $D03C12err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $D03C44 </td><td> $C44D03d </td><td> $D03C44err </td></tr>
		<tr></tr>
		<tr> <th rowspan="3"> G1 </td><td> C<sub>11</sub> </td><td> $G1C11 </td><td> $C11G1d </td><td> $G1C11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $G1C12 </td><td> $C12G1d </td><td> $G1C12err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $G1C44 </td><td> $C44G1d </td><td> $G1C44err </td></tr>
		<tr></tr>
		<tr> <th rowspan="9"> &alpha;'' </td><td> C<sub>11</sub> </td><td> $APPC11 </td><td> $C11APPd </td><td> $APPC11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $APPC12 </td><td> $C12APPd </td><td> $APPC12err </td></tr>
		<tr> <td> C<sub>13</sub> </td><td> $APPC13 </td><td> $C13APPd </td><td> $APPC13err </td></tr>
		<tr> <td> C<sub>22</sub> </td><td> $APPC22 </td><td> $C22APPd </td><td> $APPC22err </td></tr>
		<tr> <td> C<sub>23</sub> </td><td> $APPC23 </td><td> $C23APPd </td><td> $APPC23err </td></tr>
		<tr> <td> C<sub>33</sub> </td><td> $APPC33 </td><td> $C33APPd </td><td> $APPC33err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $APPC44 </td><td> $C44APPd </td><td> $APPC44err </td></tr>
		<tr> <td> C<sub>55</sub> </td><td> $APPC55 </td><td> $C55APPd </td><td> $APPC55err </td></tr>
		<tr> <td> C<sub>66</sub> </td><td> $APPC66 </td><td> $C66APPd </td><td> $APPC66err </td></tr>
		<tr></tr>
		<tr> <th rowspan="6"> &omega; </td><td> C<sub>11</sub> </td><td> $omgC11 </td><td> $C11omgd </td><td> $omgC11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $omgC12 </td><td> $C12omgd </td><td> $omgC12err </td></tr>
		<tr> <td> C<sub>13</sub> </td><td> $omgC13 </td><td> $C13omgd </td><td> $omgC13err </td></tr>
		<tr> <td> C<sub>33</sub> </td><td> $omgC33 </td><td> $C33omgd </td><td> $omgC33err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $omgC44 </td><td> $C44omgd </td><td> $omgC44err </td></tr>
		<tr> <td> C<sub>66</sub> </td><td> $omgC66 </td><td> $C66omgd </td><td> $omgC66err </td></tr>
		<tr></tr>
	</table>
	
	<h5> 50/50 Phases </h5>
	<table border="1" width="50%" align="center">
		<tr> <th></th><th>[GPa]</th><th> $typ </th><th> PAW-PBE </th><th> Error (%) </th> </tr>
		<tr></tr>
		<tr> <th rowspan="3"> B2 </td><td> C<sub>11</sub> </td><td> $B2C11 </td><td> $C11B2d </td><td> $B2C11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $B2C12 </td><td> $C12B2d </td><td> $B2C12err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $B2C44 </td><td> $C44B2d </td><td> $B2C44err </td></tr>
		<tr></tr>
	</table>

	<h5> Pure Phases </h5>
	<table border="1" width="50%" align="center">
		<tr> <th></th><th>[GPa]</th><th> $typ </th><th> PAW-PBE </th><th> Error (%) </th> </tr>
		<tr></tr>
		<tr> <th rowspan="3"> BCC Nb </td><td> C<sub>11</sub> </td><td> $BCCnbC11 </td><td> $C11BCCnbd </td><td> $BCCnbC11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $BCCnbC12 </td><td> $C12BCCnbd </td><td> $BCCnbC12err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $BCCnbC44 </td><td> $C44BCCnbd </td><td> $BCCnbC44err </td></tr>
		<tr></tr>
		<tr> <th rowspan="3"> BCC Ti </td><td> C<sub>11</sub> </td><td> $BCCtiC11 </td><td> $C11BCCtid </td><td> $BCCtiC11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $BCCtiC12 </td><td> $C12BCCtid </td><td> $BCCtiC12err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $BCCtiC44 </td><td> $C44BCCtid </td><td> $BCCtiC44err </td></tr>
		<tr></tr>
		<tr> <th rowspan="5"> HCP Ti </td><td> C<sub>11</sub> </td><td> $HCPtiC11 </td><td> $C11HCPtid </td><td> $HCPtiC11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $HCPtiC12 </td><td> $C12HCPtid </td><td> $HCPtiC12err </td></tr>
		<tr> <td> C<sub>13</sub> </td><td> $HCPtiC13 </td><td> $C13HCPtid </td><td> $HCPtiC13err </td></tr>
		<tr> <td> C<sub>33</sub> </td><td> $HCPtiC33 </td><td> $C33HCPtid </td><td> $HCPtiC33err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $HCPtiC44 </td><td> $C44HCPtid </td><td> $HCPtiC44err </td></tr>
		<tr></tr>
		<tr> <th rowspan="5"> &omega; Ti </td><td> C<sub>11</sub> </td><td> $HCPtiC11 </td><td> $C11HCPtid </td><td> $HCPtiC11err </td></tr>
		<tr> <td> C<sub>12</sub> </td><td> $HCPtiC12 </td><td> $C12HCPtid </td><td> $HCPtiC12err </td></tr>
		<tr> <td> C<sub>13</sub> </td><td> $HCPtiC13 </td><td> $C13HCPtid </td><td> $HCPtiC13err </td></tr>
		<tr> <td> C<sub>33</sub> </td><td> $HCPtiC33 </td><td> $C33HCPtid </td><td> $HCPtiC33err </td></tr>
		<tr> <td> C<sub>44</sub> </td><td> $HCPtiC44 </td><td> $C44HCPtid </td><td> $HCPtiC44err </td></tr>
		<tr></tr>
	</table>

	<h5> C<sub>ij</sub> vs P </h5>
	<table>
		<tr>
			<td> <img src="./D03_CVSP.png" style="width:100%"> </td>
			<td> <img src="./G1_CVSP.png" style="width:100%"> </td>
		</tr><tr>
			<td> <img src="./APP_CVSP.png" style="width:100%"> </td>
			<td> <img src="./OMG_CVSP.png" style="width:100%"> </td>
		</tr><tr>
			<td> <img src="./BCCnb_CVSP.png" style="width:100%"> </td>
			<td> <img src="./BCCti_CVSP.png" style="width:100%"> </td>
		</tr><tr>
			<td> <img src="./HCPti_CVSP.png" style="width:100%"> </td>
			<td> <img src="./OMGti_CVSP.png" style="width:100%"> </td>
		</tr>
	</table>

 <!--	<img src="./D03_CVSP.png" style="width:75%">

	<h4> Elastic Anisotropy </h4>

	<h5> G1 </h5>
	<img src="./G1_anisotropy.png" style="width:75%"> -->

	<h4> Phase Transitions </h4>

	<h5> &beta; to &omega; </h5>
	<table> <tr>
	<td><img src="./BETA_TO_OMEGA.png" style="width:100%"></td>
	<td><img src="./BETA_TO_OMEGA_Ti.png" style="width:100%"></td>
	</tr> </table>

	<h5> &alpha; Distortion </h5>
	<img src="./ALPHAS.png" style="width:75%">
	
	<h5> &alpha; to &beta; </h5>
	<img src="./ALPHA_BETA_BURGERS.png" style="width:75%">

	<h5> &beta;' (BCT) Deformation </h5>
	<img src="./BCT_DEFORMATION.png" style="width:75%">

	<h4> Phonon Dispersion </h4>
	<table>
	<tr> 
		<td><img src="./bccNb_phonons.png" style="width:100%"></td>
		<td><img src="./G1_phonons.png" style="width:100%"></td>
	</tr>
	<tr> 
		<td><img src="./bccTi_phonons.png" style="width:100%"></td>
		<td><img src="./D03_phonons.png" style="width:100%"></td>
	</tr>
	<tr> 
		<td><img src="./omgTi_phonons.png" style="width:100%"></td>
		<td><img src="./omega_phonons.png" style="width:100%"></td>
	</tr>
	<tr> 
		<td><img src="./hcpTi_phonons.png" style="width:100%"></td>
		<td><img src="./APP_phonons.png" style="width:100%"></td>
	</tr>
	<tr> 
		<td><img src="./B2_phonons.png" style="width:100%"></td>
		<td><img src="./A3_phonons.png" style="width:100%"></td>
	</tr>
	<tr>
		<td><img src="./A15TiNb3_phonons.png" style="width:100%"></td>
	</tr>
	</table>
	
	<h4> Mixing enthalpies </h4>

	<img src="./THERMO.png" style="width:75%">

	<h4> Stacking Faults </h4>

	<table style="width:100%"> <tr>
		<td><img src="./SF_BCC_Nb.png" style="width:100%"></td>
		<td><img src="./SF_HCP_011b0_112b0.png" style="width:100%"></td>
	</tr></table>

	<h4> Point defects </h4>

	<table border="1" width="75%" align="center">
		<tr><th colspan="5"> B2 TiNb intrinsic and complex defect formation enthalpies in eV </th></tr>
		<tr><th> Defect </th><th> Symbol </th><th> MEAM </th><th> PAW-PBE </th><th> Error (%) </th></tr>
		<tr><th> Nb in Ti antisite </th><th> Nb<sub>Ti</sub> </th>
			<td> `printf %.3f ${B2_pds['Ti-antisite']}` </td><td> `printf %.3f ${B2_pdd['Ti-antisite']}`</td><td> ${B2_pderr['Ti-antisite']} </td></tr>
		<tr><th> Ti in Nb antisite </th><th> Ti<sub>Nb</sub> </th>
			<td> `printf %.3f ${B2_pds['Nb-antisite']}` </td><td> `printf %.3f ${B2_pdd['Nb-antisite']}`</td><td> ${B2_pderr['Nb-antisite']} </td></tr>
		<tr><th> Vacancy in Ti </th><th> Va<sub>Ti</sub> </th>
			<td> `printf %.3f ${B2_pds['Ti-vacancy']}` </td><td> `printf %.3f ${B2_pdd['Ti-vacancy']}`</td><td> ${B2_pderr['Ti-vacancy']} </td></tr>
		<tr><th> Vacancy in Nb </th><th> Va<sub>Nb</sub> </th>
			<td> `printf %.3f ${B2_pds['Nb-vacancy']}` </td><td> `printf %.3f ${B2_pdd['Nb-vacancy']}`</td><td> ${B2_pderr['Nb-vacancy']} </td></tr>
		<tr><th> Divacancy </th><th> Va<sub>Ti</sub> + Va<sub>Nb</sub> </th>
			<td> `printf %.3f ${B2_pds['divacancy']}` </td><td> `printf %.3f ${B2_pdd['divacancy']}`</td><td> ${B2_pderr['divacancy']} </td></tr>
		<tr><th> Antisite exchange </th><th> Ti<sub>Nb</sub> + Nb<sub>Ti</sub> </th>
			<td> `printf %.3f ${B2_pds['exchange']}` </td><td> `printf %.3f ${B2_pdd['exchange']}`</td><td> ${B2_pderr['exchange']} </td></tr>
		<tr><th> Triple Ti </th><th> 2Va<sub>Ti</sub> + Ti<sub>Nb</sub> </th>
			<td> `printf %.3f ${B2_pds['triple-Ti']}` </td><td> `printf %.3f ${B2_pdd['triple-Ti']}`</td><td> ${B2_pderr['triple-Ti']} </td></tr>
		<tr><th> Triple Nb </th><th> 2Va<sub>Nb</sub> + Nb<sub>Ti</sub> </th>
			<td> `printf %.3f ${B2_pds['triple-Nb']}` </td><td> `printf %.3f ${B2_pdd['triple-Nb']}`</td><td> ${B2_pderr['triple-Nb']} </td></tr>
		<tr><th> Interbranch Ti excitation </th><th> 2Va<sub>Nb</sub> - Ti<sub>Nb</sub> </th>
			<td> `printf %.3f ${B2_pds['inter-Ti']}` </td><td> `printf %.3f ${B2_pdd['inter-Ti']}`</td><td> ${B2_pderr['inter-Ti']} </td></tr>
		<tr><th> Interbranch Nb excitation </th><th> Nb<sub>Ti</sub> - 2Va<sub>Ti</sub> </th>
			<td> `printf %.3f ${B2_pds['inter-Nb']}` </td><td> `printf %.3f ${B2_pdd['inter-Nb']}`</td><td> ${B2_pderr['inter-Nb']} </td></tr>
	</table>

</div>

<div id="footer">
	Last modified: <script> document.write(document.lastModified); </script>
	<img src="../../../../../home/physdept.gif" style="width:200px;height:auto" align="center">
</div>
</div>

</body>
</html>
!!

done # end pot loop
done # end set loop

cd ~/public_html/POTS/$elem1-$elem2/
	./rehash.sc
cd $wd
echo "ALL DONE!"
