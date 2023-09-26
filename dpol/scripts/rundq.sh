prog=dnewdq.exe
wd=`pwd | tail -c 3`
echo $wd
cp ../$prog $prog$wd
nohup $prog$wd &
date >>nohup.out
echo Program $prog$wd started.
echo job number = $! > jobno
cat jobno
exit

