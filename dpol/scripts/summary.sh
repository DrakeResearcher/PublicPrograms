# convert a .POW file into the corresponding .POL file and store in dpoldl.dat
case $1 in
  -c) set `cat current.dat`
      echo input to summary from current.dat = $*
      if [ ! -f notdone ]; then echo " "$1 $2 $3 $4 POL >current.dat;fi
      ;;
   *) if [ "$1" = "" ]; then
         echo usage: summary Z n S L PO[WL] or
         echo "       "summary -c to use current values in current.dat
         exit
      fi
      ;;
esac

set `echo $* | tr '[a-z]' '[A-Z]'`
#if [ ! -f $1$2$3$4???.$5 ]; then
#  echo "-----------------------------------------------"
#  echo "   From summary: File $1$2$3$4???.$5 not found "
# echo "-----------------------------------------------"
# exit 1
#fi
if [ ! -f notdone ]; then 
   echo " "$1 $2 $3 $4 $5 summary >current.dat
fi
echo Z = $1 $2 $3$4 $5 BASIS SETS >summary.dat
#  echo .T. to include mass scaling of the Rydberg.
echo .F. >>summary.dat
#ls -l $1$2$3$4???.$5 |cut -c55- >>summary.dat
for i in $1$2$3$4[0-9][0-9][0-9].$5; do 
   echo $i >>summary.dat
done
for i in $1$2$3$4[0-9][0-9][0-9][0-9].$5; do 
   echo $i >>summary.dat
done

# Removed as we no longer use this naming scheme
#==========================================================#
# for i in $1$2$3$4[A-Z][0-9][0-9][0-9][0-9].$5; do  
#    echo ${#i} " " $i 
#    echo $i >>summary.dat
# done
#==========================================================#


echo EXIT >>summary.dat
#date |sed "s/\ CST//"|cut -c5-11,21-24 >DATE.DAT
date +"%b %d %Y" >DATE.DAT
if [ -f $1$2$3$4$5.MAT ]; then rm $1$2$3$4$5.?AT; fi
summarv.exe
