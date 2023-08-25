#!/bin/bash
#run two scripts on 100K datasets to assess time difference

#v1
time python ../bin/demux_index_optim.py --fastq_path '/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/041323_RBenjamin_gigaassay/clustering_test/SRR20707787_0.1.fastq.gz' \
--barcode_path '/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/041323_RBenjamin_gigaassay/clustering_test/starcode_SRR20707797_0.1_cluster.txt.gz' \
--sample_name 'v1_test' --output_dir ./v1_test

#v2
time python ../bin/demux_index_optim_v2.py --fastq_path '/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/041323_RBenjamin_gigaassay/clustering_test/SRR20707787_0.1.fastq.gz' \
--barcode_path '/Users/martingordon/Documents/projects/041223_Ronald_B_Babu/041323_RBenjamin_gigaassay/clustering_test/starcode_SRR20707797_0.1_cluster.txt.gz' \
--sample_name 'v2_test' --output_dir ./v2_test
