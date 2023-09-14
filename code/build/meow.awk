#!/bin/awk -f
# Concatenate based on Delimiters
# https://learnxinyminutes.com/docs/awk/
# http://www.grymoire.com/Unix/Awk.html
# https://stackoverflow.com/questions/17988756/how-to-select-lines-between-two-marker-patterns-which-may-occur-multiple-times-w

BEGIN {
	counter = -1
	dflag = 0
}
/\#>>>/ {
	dflag = 1
	counter = counter + 2
	next
}

/\#<<</ {
	dflag = 0
	next
}

dflag {
    if (counter < 10)
    	print $0 > "0"counter"_join.R"
    else
        print $0 > counter"_join.R"
}
