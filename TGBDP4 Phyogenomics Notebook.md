# Analyzing TgBDP4 through BLAST and Phylogenomic Analyses to Access Evolutionary Relationships Amongst Apicomplexians
This is a workbook detailing the methods and steps utilized in this bioinformatic project to access if the TgBDP4 sequence is conserved across the Apicomplexan phylum.
# Background

# Materials 
  * PC: Lenova Yoga 730
  * Coding Terminal: PuTTy
  * Preformance Cluster : Ron
  * Programs Utalized:
      * MAFFT
      * iqTree
      * FigTree v1.4.4 ( https://github.com/rambaut/figtree/releases)
  * Databases Utalized:
  	* VEuPathDB (https://veupathdb.org/veupathdb/app)
	* UniProtKB (https://www.uniprot.org/)
# Results from BLAST Analysis
BLAST is a local alignment anlysis which was preformed through EUPATHdb using the predicted protein sequence for TgBDP4. This BLAST analysis was preformed using a BLOSUM62 Matrix to find similar protein sequences within the Apicomplexians. This BLAST analysis yielded the top 50 apicomplexans that appear to have a similar protein sequence within their genome. This list was further narrowed down through the elimination of multiple samples from the same species.
These sequences are:
* Toxoplasma gondii ME49 : TGME49_306460
* Hammondia hammondi strain H.H.34: HHA_306460
* Neospora caninum Liverpool: 	NCLIV_044660
* Besnoitia besnoiti strain Bb-Ger1: BESB_018420
* Cystoisospora suis strain Wien I : CSUI_003223
* Sarcocystis neurona SO SN1: SRCN_4078
* Eimeria mitis Houghton: EMH_0040220
* Eimeria necatrix Houghton: ENH_00043160
* Eimeria praecox Houghton: EPH_0031840
* Cyclospora cayetanensis strain CHN_HEN01: cyc_01595
* Eimeria acervulina Houghton: 	EAH_00010620
* Cryptosporidium ubiquitum isolate 39726: cubi_02786
* Cryptosporidium parvum IOWA-ATCC: CPATCC_0006870
* Cryptosporidium tyzzeri isolate UGA55: CTYZ_00000707
* Cryptosporidium hominis isolate TU502_2012: ChTU502y2012_407g0560
* Cryptosporidium meleagridis strain UKMEL1: 	CmeUKMEL1_10795
* Cryptosporidium andersoni isolate 30847: cand_029210
* Plasmodium falciparum 3D7: PF3D7_1475600
* Plasmodium praefalciparum strain G01: PPRFG01_1476100
     
# Setting up for the Phylogenomic Analysis:
First I created a new directory to begin my analysis in:
	
	mkdir TgBDP4
		
Next, I moved into the newly created directed to begin my analysis:
	
	cd TgBDP4 

# Phylogenomics Analysis of TGBDP4 predicted protein sequence within Aplicomplexians 
A phylogenomic analysis is being utilized in this study to compare the relationships of apicomplexans to Toxoplasma gondii's TgBDP4. The amplicomplexans being compared to the T.gondii's TgBDP were chosen based on a previous blast analysis.

	# Phylogenomics database
		#Toxoplasma gondii ME49
			 awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TGME49_306460.fa > ToxoME49.fa
		# Hammondia hammondi strain H.H.34
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Neospora caninum Liverpool
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Besnoitia besnoiti strain Bb-Ger1
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cystoisospora suis strain Wien I
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Sarcocystis neurona SO SN1
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Eimeria mitis Houghton
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Eimeria necatrix Houghton
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Eimeria praecox Houghton
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cyclospora cayetanensis strain CHN_HEN01
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Eimeria acervulina Houghton
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cryptosporidium ubiquitum isolate 39726
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cryptosporidium parvum IOWA-ATCC
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cryptosporidium tyzzeri isolate UGA55
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cryptosporidium hominis isolate TU502_2012
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cryptosporidium meleagridis strain UKMEL1
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Cryptosporidium andersoni isolate 30847
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		# Plasmodium falciparum 3D7
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
		#  Plasmodium praefalciparum strain G01
			curl -LO 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'
Once all of the sequences where downloaded and formated, the files are then combined `cat` into a single file:
			
	Cat   > TGBDP4_Phylogenetics_Tree.fa
		
A MAFFT analysis is then completed on the combined data file:
	
	mafft --auto TGBDP4_Phylogenetics_Tree.fa > Output_TGBDP4_Phylogenetics_Tree.fa
Once the MAFFT analysis is complete, the phylogenomic tree was built utilizing iqtree.

	iqtree -s Output_TGBDP4_Phylogenetics_Tree.fa -m JC -bb 1000 -pre output
`JC` is the model utilized to create the phylogenomic tree, which is based on equal substitution rate and equal base frequencies. 

`-bb` is the bootstraping subset utilized to create this phylogenomic tree(1000). 

This iqtree analysis outputs a Tiled called : `output.contree`

	Cat output.contree
	
The output from the cat command will then be copy/pasted into Figtree to yield the phylogenetic trees.

# Results of the Phylogenomics Analysis of TGBDP4 predicted protein sequence within Aplicomplexians 
