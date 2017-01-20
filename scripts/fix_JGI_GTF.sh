for file in `cat fix.txt`; do 
 pref=$(echo $file | perl -p -e 's/(\S+)\.(\S+)\.v\d+(\.\d+)?\.gff3/$2/')
 perl  ~/src/genome-scripts/data_format/gtf2gff3_3level.pl --JGI -p $pref $file > ../$file
done
