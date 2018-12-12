Inp=$1
OutFil=`basename $Inp|sed s/.avinp//g`
echo $OutFil
table_annovar.pl $Inp /share/data/resources/hg19/ANNOVAR_humandb/ -buildver hg19 -out $OutFil -remove -protocol gnomad_exome,dbnsfp35a -operation f,f -nastring . -csvout -polish -xref $OutFil.anno.csv
