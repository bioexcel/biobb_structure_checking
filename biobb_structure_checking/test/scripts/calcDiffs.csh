#!/bin/csh
echo "Calculating diffs"
foreach f (*log *json *pdb)
diff -Zb -I Andrej -I Processor -I MODELLER -I 'Date and time of compilation' -I 'Job starting time' $f ref/$f > $f.diff
if !(-z $f.diff) then
	echo $f.diff
	more $f.diff
endif
end

