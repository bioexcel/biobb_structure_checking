#!/bin/csh
echo "Calculating diffs"
foreach f (*pdb *json *log)
diff $f ref/$f > $f.diff
if !(-z $f.diff) then
	echo $f.diff
	more $f.diff
endif
end

