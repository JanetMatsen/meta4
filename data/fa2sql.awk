BEGIN {
	s = -1;
	if (table == "")
		table = "sequences";
	if (create_table == 1)
		printf("DROP TABLE IF EXISTS %s; CREATE TABLE %s ( locus VARCHAR(32) PRIMARY KEY, description VARCHAR(256), sequence LONGTEXT );", table, table);
}

/^>/ {
	++s;
	sub(">", "", $1);
	locus[s] = $1;
	description[s] = $2;
	for (i = 3; i <= NF; ++i)
		description[s] = description[s] " " $i; 
	sequence[s] = "";
}
{
	if (substr($1, 1, 1) != ">")
		sequence[s] = sequence[s] $1;
}
END {
	for (i = 0; i <= s; ++i) { 
		printf("INSERT INTO %s (locus, description, sequence) VALUES (\"%s\", \"%s\", \"%s\");\n", table, locus[i], description[i], sequence[i]);
	}
} 
