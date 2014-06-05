awk -v dir=$2 ' BEGIN {
	chrName=0
}
{
	if ($1 == ">") {
		if (chrName != 0) {
			close(dir"/"chrName".fa")
		}
		chrName=$2
		print chrName
	} else {
		x=split($1,nam,">")
		if (x > 1) {
			if (chrName != 0) {
				close(dir"/"chrName".fa")
			}
			chrName=nam[2]
			print chrName
		}
	}
	print $0 > dir"/"chrName".fa"
} 
END {
	close(dir"/"chrName".fa")
}' $1
	