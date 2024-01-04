while read file;do 
 	media=$(basename $(dirname $file));
 	model=$(basename $file|sed 's/.xml_summary.txt//g');
 	paste $file|grep -v Reactions|grep -v Metabolites|sed 's/ /\t/g'|sed "s/^/${model}\t${media}\t/g";
done< <(find models/ -name "*summary.txt") > model_summary.tsv