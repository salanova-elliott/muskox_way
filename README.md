# Scripts for metagenomic greenland
Scripts used in the processing of metagenomic data from

**The Muskox Way: Humans, plants and animals in northeast Greenland through the Holocene**
Inger Greve Alsos, Marieke Beaulieu, Nickolas L. Balascio, Lucas Elliott, Raymond S. Bradley, Sergi Pla, 5, Santiago Giralt, Luisa Deppe, Youri Lammers, Jamie Alumbaugh, Peter D. Heintzman, Love Dalén, Sandra Nogué, Frederik Fuuja Larsen, Christian Koch Madsen, Bianca Perren, William J. D'Andrea, Jostein Bakke

*sam_filter_mismatch.py* : Uses MD:Z string output from bowtie2 to calculate percent identity of each mapping and removes matches < 0.97.

*lca_bowtie_continuous.py* : Outputs least common ancestor of each read from all mappings. Requires a headless, sorted sam file and rocksDB of accessionID : ncbi taxon ID.

*kmer_count.py* : Calculates the percentage of unique k-mers 
