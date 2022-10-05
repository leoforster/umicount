# UMIcount

A tool for counting umis in SmartSeq3 libraries

# Concept

This tool works with two cli scripts: `process_smartseq` and `umicount`.
First `process_smartseq` is run on the fastq files to extract umis and attach them to readnames.
The resulting fastq files are run through a normal alignment pipeline.
After this, resulting alignments are quantified using `umicount`.

