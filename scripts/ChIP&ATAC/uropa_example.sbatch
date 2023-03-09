#! /bin/bash
#SBATCH --job-name=uropa
#SBATCH --time=24:00:00
#SBATCH --mem=24g

DIR1=INFOLDER
DIR2="Uropa"

module load uropa/4.0.2

RE="$DIR1/macsNarrow/(.*)/.*narrowPeak"

for BED in $DIR1/macsNarrow/*/*narrowPeak;
do

if [[ $BED =~ $RE ]]; then
  sample=${BASH_REMATCH[1]}
fi

out="$DIR1/$DIR2/${sample}_macsNarrow_uropa4_protTSS"

echo '{"queries":[ ' > $DIR1/$DIR2/$sample.protTSS.json
echo '{ "feature":"gene","distance":[3000,1000],"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start"},' >> $DIR1/$DIR2/$sample.protTSS.json
echo '{ "feature":"gene","distance":10000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start"},' >> $DIR1/$DIR2/$sample.protTSS.json
echo '{ "feature":"gene","distance":100000,"filter.attribute":"gene_type","attribute.value":"protein_coding","feature.anchor":"start"}],' >> $DIR1/$DIR2/$sample.protTSS.json
echo '"show_attributes":["gene_id", "gene_name","gene_type"],' >> $DIR1/$DIR2/$sample.protTSS.json
echo '"priority":"Yes",' >> $DIR1/$DIR2/$sample.protTSS.json
echo '"gtf":"genes.gtf",' >> $DIR1/$DIR2/$sample.protTSS.json
echo "'bed': '$BED' }" | sed 's/'\''/'\"'/g' >> $DIR1/$DIR2/$sample.protTSS.json

uropa -i $DIR1/$DIR2/$sample.protTSS.json -p $out -s

done
