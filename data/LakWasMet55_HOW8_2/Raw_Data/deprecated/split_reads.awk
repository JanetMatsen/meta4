BEGIN {
        line = 0;
        r1File = "r1.fastq";
        r2File = "r2.fastq";
        currentFile = r1File;
}
{
        printf("%s\n", $0) >> currentFile;
        ++line;
        if (line % 8 >= 4) {
                currentFile = r2File;
        } else {
                currentFile = r1File;
        }
		if (line % 40000 == 0)
			printf("processed %d reads\n", line / 4);
}
END {
	printf("found %d total reads\n", line / 4);
}
