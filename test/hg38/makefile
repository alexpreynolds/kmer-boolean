all:

parallelize:
	./parallelize.kmer-boolean.sh ${PWD}/../../kmer-boolean ${PWD}/results /net/seq/data/genomes/human/GRCh38/noalts/chromosomes ${PWD}

collapse_times:
	./collapse-times.kmer-boolean.sh ${PWD}/results/logs ${PWD}/results/logs-collapsed ${PWD}

graph:
	./graph.Rscript
