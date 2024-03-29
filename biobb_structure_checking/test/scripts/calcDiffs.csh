#!/bin/csh
echo "Calculating diffs"
foreach f (*log *json *pdb *pqr *pdbqt *cmip *fasta)
diff -b -I Andrej -I Processor -I MODELLER -I 'Date and time of compilation' -I 'Job starting time' -I '#DEBUG' $f ref/$f > $f.diff
if !(-z $f.diff) then
	echo $f.diff
	more $f.diff
endif
end

