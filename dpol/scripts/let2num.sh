#  script to change file names with letter basis et izes to numerical basis sizes
#  for all files in the directory of the form *.PO?
for f in *.PO?; do
head=`echo $f |cut -c 1-4`
tail=`echo $f |cut -c 6-12`
x=`echo $f |cut -c 5`
p=""
n=$x
# get character code using printf %d\\n \'$x (drop \\n for "new line")
gx=`printf %d \'$x`
if [ $gx -gt 64 ]; then
  if [ $gx -lt 91 ]; then
  if [[ "$x" < "K" ]]; then
   n=`echo $x |tr '[A-J}' '[0-9]'`
   p="1"
  elif [[ "$x" > "J" ]] && [[ "$x" < "U" ]]; then
   n=`echo $x |tr '[K-T]' '[0-9]'`
   p="2"
  elif [[ "$x" > "T" ]]; then
   n=`echo $x |tr '[U-Z]' '[0-5]'`
   p="3"
  fi
else
  x=`echo $x | tr '[a-z]' '[A-Z]'`
  if [[ "$x" < "E" ]] && [[ $gx -gt 96 ]]; then
   n=`echo $x |tr '[A-D]' '[6-9]'`
   p="3"
   echo hello
  elif [[ "$x" > "D" ]] && [[ "$x" < "O" ]]; then
   n=`echo $x |tr '[E-N]' '[0-9]'`
   p="4"
  elif [[ "$x" > "N" ]] && [[ "$x" < "Y" ]]; then
   n=`echo $x |tr '[O-X]' '[0-9]'`
   p="5"
  elif [[ "$x" > "X" ]]; then
   n=`echo $x |tr '[Y-Z]' '[0-1]'`
   p="6"
  fi
fi
if [ ! "$p" = "" ]; then
echo mv $f $head$p$n$tail
mv $f $head$p$n$tail
fi
fi
done
