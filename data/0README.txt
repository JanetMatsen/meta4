# original files as received from JGI
isolate_genomes.original.fasta
isolate_genomes.original.gbk

# fixed files with duplicate copies of scaffold1 removed and renamed
isolate_genomes.fasta
isolate_genomes.gbk

# this runs fixup_fast and fixup_genbank and saves output
preprocess_input

# the GFF3 file for the scaffolds
isolate_genomes.gbk.gff
# this generates the gff from the gbk using bioperl
setup_genes_table

# data sets
LakWasMet55_HOW8_2

# index genome for bbmap.  Makes a dir called ref.
bbmap_index.sh

# order of operations:
./preprocess_input
./setup_genes_table
./bbmap_index.sh
