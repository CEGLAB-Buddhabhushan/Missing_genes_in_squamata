curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat


for i in `cat Squamata_Reference_scaffold_chromososme.sorted.lst`
do
accession=`echo $i|cut -f1 -d'-'`
Species_name=`echo $i|cut -f2 -d'-'`

./datasets download genome accession $accession --include genome --filename "$i".zip
rm -rf $i
unzip -q -o "$i".zip -d $i
rm "$i".zip
cd $i/ncbi_dataset/data/$accession/

/media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_genome_blast/BLASTn.sh $i

cd /media/lokdeep/sdf/BUDDHA/BMCBIO/IL34_genome_blast/
rm -r $i
echo $Species_name
done


