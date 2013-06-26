ARCHIVE=`awk '/^__ARCHIVE_BELOW__/ {print NR + 1; exit 0; }' $0`
TMPDIR=`mktemp -d`

echo ""
echo ""
echo ""
echo "Extracting"

tail -n+$ARCHIVE $0 | base64 -d | tar xzv -C $TMPDIR

CDIR=`pwd`

$PTEST_BINARY $CLOPTS > $TMPDIR/stdout.test 

#out="$(diff -ur $TMPDIR/stdout.test $TMPDIR/stdout.gold)"
ktest=`tail -n1 $TMPDIR/stdout.test | cut -d' ' -f6`
kgold=`tail -n1 $TMPDIR/stdout.gold | cut -d' ' -f6`

echo ""
echo ""
echo ""
echo "Expected"
cat $TMPDIR/stdout.gold

echo ""
echo ""
echo ""
echo "Got"
cat $TMPDIR/stdout.test

rm -rf $TMPDIR

if [[ $ktest == $kgold ]]
then
exit 0
else
exit 1
fi

__ARCHIVE_BELOW__
