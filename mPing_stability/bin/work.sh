perl confident2gff.pl --mping ../input/10092013.mpings.txt
python AlleleFrq.py --input ../input/10092013.mpings.gff > mping.ril.frquency

grep "Reference-only" ../input/HEG4.mping.all_inserts.gff > ../input/HEG4.mping.reference.gff


python Excision.py --input ../input/10092013.mpings.gff > excision.txt
python NonRef_RILmPing.py --input ../input/10092013.mpings.gff | sort -k1,1 -k3,3n > commonNonrefmPing.list

