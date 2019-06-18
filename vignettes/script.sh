for f in *.pdf
do
    echo $f;
    convert -compress Zip $f  $f;
done