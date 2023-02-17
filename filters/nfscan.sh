make
for i in $(seq -12 6 36)
do
  testfilter nfscan lp2 $i > S$i
  testfilter lnfscan lp4 $i > L$i
done
