DROP TABLE IF EXISTS summary_LakWasMet55_HOW8_2_bwa;
CREATE TABLE summary_LakWasMet55_HOW8_2_bwa (
        locus_tag VARCHAR(32) PRIMARY KEY,
	reads_mapped FLOAT,
	rpkm FLOAT
);

LOAD DATA INFILE "/gscratch/lidstrom/meta4/workspace/LakWasMet55_HOW8_2/bwa/LakWasMet55_HOW8_2.summary.dat" INTO TABLE summary_LakWasMet55_HOW8_2_bwa FIELDS TERMINATED BY '\t' (locus_tag, reads_mapped);

# count reads mapped to coding sequences
SELECT @sum_to_CDS:=SUM(reads_mapped) FROM summary_LakWasMet55_HOW8_2_bwa AS s INNER JOIN genes_isolate_genomes AS g ON s.locus_tag=g.locus_tag WHERE g.type='CDS';

# compute RPKM
UPDATE summary_LakWasMet55_HOW8_2_bwa AS s INNER JOIN genes_isolate_genomes AS g ON s.locus_tag=g.locus_tag SET rpkm=
        reads_mapped /
# length of gene in kilobases
                (((g.end_coord + 1) - g.start_coord) / 1000) /
# millions of reads mapped to coding sequences
                (@sum_to_CDS / 1000000)
;
