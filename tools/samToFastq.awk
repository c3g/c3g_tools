
awk ' BEGIN { FS="\t" } { 
	print "@"$1; 
	print $10; 
	print "+";  
	print $11;
}'

