# Analyzing TgBDP4 through BLAST and Phylogenomic Analyses to Access Evolutionary Relationships Amongst Apicomplexians

# Materials
  * PC: Lenova Yoga 730
  * Coding Terminal: PuTTy
  * Preformance Cluster : Ron
  * Programs Utalized:
      * MAFFT
      * iqTree
      * FigTree v1.4.4
  * Databases Utalized:
  	* EUPATHdb
	* Uniprotkb
# Results from BLAST Analysis
BLAST is a local alignment anlysis which was preformed through EUPATHdb using the predicted protein sequence for TgBDP4. This BLAST analysis was preformed using a BLOSUM62 Matrix to find similar protein sequences within the Apicomplexians. This BLAST analysis yielded the top 50 apicomplexans that appear to have a similar protein sequence within their genome. This list was further narrowed down through the elimination of multiple samples from the same species.
These sequences are:
     
# Setting up for the Phylogenomic Analysis:
First I created a new directory to begin my analysis in:
	
	mkdir TgBDP4
		
Next, I moved into the newly created directed to begin my analysis:
	
	cd TgBDP4 

# Phylogenomics Analysis
A phylogenomic analysis is being utilized in this study to compare the relationships of apicomplexans to Toxoplasma gondii's TgBDP4. The amplicomplexans being compared to the T.gondii's TgBDP were chosen based on a previous blast analysis.
	
	# Phylogenomics database
		#TGME49_ChrVIIb
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		
			
	Cat   > TGBDP4_Phylogenetics_Tree.fa
		
		
A MAFFT analysis is then completed on all the data downloading
	
	mafft --auto TGBDP4_Phylogenetics_Tree.fa > Output_TGBDP4_Phylogenetics_Tree.fa
Once the MAFFT analysis is complete, the phylogenomic tree can be built utilizing iqtree.

	iqtree -s Output_TGBDP4_Phylogenetics_Tree.fa -m JC -bb 1000 -pre output
	Cat output.contree
The output from the cat command will then be copy/pasted into Figtree to yield the phylogenetic trees.
