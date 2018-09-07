# Running TRAPeS on other genomes

TRAPeS can be used for any genome and any organism, as long as the genomic annotations for the V/J/C segments are available. In order to run TRAPeS on genomes besides mm10, hg19 or hg38, the user must create those annotations and provide the correct value for the '-genome' parameter. <br /> <br />

Instructions to create new annotations for TRAPeS (see existing genome folders for examples): <br />

Under the Data folder, create a new folder. The name of the folder is the value you will provide the ‘-genome’ parameter. The folder should include only these 4 files: <br />
1.	*.gene.id.mapping.TCR.txt (file name must end with ‘gene.id.mapping.BCR.txt’): A file with no header and two columns: first column is the ID of the V/J/C segment as it appears in the *.TCR.bed file (see below) and the second column is the name of the segment (e.g. TRBV20-1, TRAJ54 etc.). This file can be created, for example, from the ensemble download page, selecting only Transcript stable ID and Gene name, and then filtering the file to include only V/J/C segments. No need to include the beta chain ‘D’ segment annotations.     
2.	*.TCR.bed: Genomic annotations in bed format of the V/J/C segments. Note that TRAPeS will only attempt to reconstruct TCRs from the V and J segments in this file.
3.	*.TCR.fa: A fasta file with the sequences of all V/J/C segments (can be generated from the TCR.bed file with the bedtools getfasta commend).
4.	*.conserved.AA.txt: A file with two columns and no header. The first column includes the name of the V/J segment (same name as in the *.gene.id.mapping.TCR.txt file). The second column includes the amino acid sequence of the protein up to the conserved amino acid from which the CDR3 is defined. I.e, for V segments it includes the amino acid sequence of the protein up until the conserved Cysteine at position 104 (including the Cysteine). For the J segments it has the amino acid sequence from the conserved W/F amino acid until the end of the J segment. The file needs to include the V and J sequences from all chains (alpha and beta).  
One way to generate these files is by copying the protein sequences from the “Protein display” page of the “IMGT repertoire” [(link)](http://www.imgt.org/IMGTrepertoire/Proteins/), and to manually arrange the sequences from all chains to fit the format above.
