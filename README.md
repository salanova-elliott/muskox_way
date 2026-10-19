# Scripts for metagenomic greenland
Scripts used in the processing of metagenomic data from

**The Muskox Way: Humans, plants and animals in northeast Greenland through the Holocene**

*sam_filter_mismatch.py* : Uses MD:Z string output from bowtie2 to calculate percent identity of each mapping and removes matches < 0.97.

*lca_bowtie_continuous.py* : Outputs least common ancestor of each read from all mappings. Requires a headless, sorted sam file and rocksDB of accessionID : ncbi taxon ID.

*kmer_count.py* : Calculates the percentage of unique k-mers 
