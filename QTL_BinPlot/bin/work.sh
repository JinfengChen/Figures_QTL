echo "bin distance matrix"
python MapDist.py --input ../input/MPR.geno.bin.uniq > ../input/MPR.geno.bin.uniq.dist

echo "bin map to genome position"
python MapGenome.py --input ../input/MPR.geno.bin.uniq

echo "draw bin map"
perl QTL_BinMap.py --input ../input/MPR.geno.bin.uniq.dist --bin ../input/MPR.geno.bin.uniq.new


echo "274"
python MapDist.py --input ../input/MPR.geno.bin.uniq > ../input/MPR.geno.bin.uniq.dist
python MapGenome.py --input ../input/MPR.geno.bin.uniq
perl QTL_BinMap.py --input ../input/MPR.geno.bin.uniq.dist --bin ../input/MPR.geno.bin.uniq.new


