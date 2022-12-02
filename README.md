# peptide-MHCassess
Application and assessment of peptide-MHC binding affinity prediction

We start with protein FASTA file downloaded from the National Center of Biotechnology Information (NCBI) Reference Sequence Database (https://www.ncbi.nlm.nih.gov/refseq/)

Using the above FASTA sequences, we used netchop v3.0 "C-term" model with a cleavage threshold of 0.1 to filter peptides that were not predicted to undergo canonical MHC class I antigen processing via proteasomal cleavage.
To assess binding affinity, we kmerized the protein FASTA into peptides of length 8-12.

From FASTA file:

``` ./preprocess_fasta_for_netchop.sh YOURFASTA.fasta OUTPUTFASTA.fasta``` 

This results in a deliminted fasta file for netchop.

Run netchop with default settings.

After netchop finishes, we postprocess the netchop output "netchop_out.txt" with: 

``` ./postprocess_netchopoutput.sh ``` 

After running postprocess_netchopoutput.sh we ran a custom python script (faster_kmer.py) to generate a list of peptides.
``` python3 faster_kmer.py > output.pep ``` 
``` sed '/|/d' output.pep > output_filt.pep ``` 


Each output file is a list of peptides that can be used for peptide-MHC binding predictions for HLAthena, MHCflurry, and MHCnuggets.
Random.pep is 1 million peptides used for analysis which did not undergo any filtering.


For peptide features analysis we sequence feature tools using:
``` python3 sequence_featurization_tools.py PEPFILE.pep OUTPUT.csv ```

Note that for producing figures, PCA and UMAP approximations are highly similar in first 2 dimensions so can be replaced using PCA_peps.R file if computational time for UMAP is too high.
