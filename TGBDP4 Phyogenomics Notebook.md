# Analyzing TgBDP4 through BLAST and Phylogenomic Analyses to Access Evolutionary Relationships Amongst Amplicomplexians

# Materials
  * PC: Lenova Yoga 730
  * Coding Terminal: PuTTy
  * Preformance Cluster : Ron
  * Programs:
      * MAFFT
      * iqTree
      * FigTree v1.4.4
     
# Setting up:
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
		#TGVEG
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#KL998088
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#KX792174.1
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#ML1130729
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#KN000275
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#KNO44666
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#KN042473
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#KN042408
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}''
		
		
A MAFFT analysis is then completed on all the data downloading
	
	mafft --auto TGBDP4_Phylogenetics_Tree.fa > Output_TGBDP4_Phylogenetics_Tree.fa
Once the MAFFT analysis is complete, the phylogenomic tree can be built utilizing iqtree.

	iqtree -s Output_TGBDP4_Phylogenetics_Tree.fa -m JC -bb 1000 -pre output
	Cat output.contree
The output from the cat command will then be copy/pasted into Figtree to yield the phylogenetic trees.
