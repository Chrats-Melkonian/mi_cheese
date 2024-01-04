while read model;do 
	media=$(basename $(dirname $model));
	species=$(basename ${model%.*});
	less $model|grep GAP_FILL -B 12|grep "reaction id"|sed 's/ fast.*$//g'|sed 's/^.*<//g'|sed 's/reaction id="//g'|sed 's/ name="//g'|sed 's/ reversible="//g'|sed 's/"$//g'|tr '"' '\t'|sed "s/^/$media\t$species\t/g";
done< <(find models/ -name "*.xml")