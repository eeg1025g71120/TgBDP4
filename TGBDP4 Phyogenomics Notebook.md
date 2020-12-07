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
			 awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TGME49_306460.fa > Toxoplasma_gondii.fa
			 awk '/^>/{print ">Toxoplasma_gondii_ME49" ++i; next}{print}' Toxoplasma_gondii.fa > header_Toxoplasma_gondii.fa
		# Hammondia hammondi strain H.H.34
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' HHA_306460.fa > Hammondia_hammondi.fa
			awk '/^>/{print ">Hammondia_hammondi" ++i; next}{print}' Hammondia_hammondi.fa > header_Hammondia_hammondi.fa
		# Neospora caninum Liverpool
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NCLIV_044660.fa > Neospora_caninum.fa
			awk '/^>/{print ">Neospora_caninum" ++i; next}{print}' Neospora_caninum.fa > header_Neospora_caninum.fa
		# Besnoitia besnoiti strain Bb-Ger1
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BESB_018420.fa > Besnoitia_besnoiti.fa
			awk '/^>/{print ">Besnoitia_besnoiti" ++i; next}{print}' Besnoitia_besnoiti.fa > header_Besnoitia_besnoiti.fa
		# Cystoisospora suis strain Wien I
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CSUI_003223.fa > Cystoisospora_suis.fa
			awk '/^>/{print ">Cystoisospora_suis" ++i; next}{print}' Cystoisospora_suis.fa > header_Cystoisospora_suis.fa
		# Sarcocystis neurona SO SN1
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SRCN_4078 > Sarcocystis_neurona.fa
			awk '/^>/{print ">Sarcocystis_neurona" ++i; next}{print}' Sarcocystis_neurona.fa > header_Sarcocystis_neurona.fa
		# Eimeria mitis Houghton
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EMH_0040220 > Eimeria_mitis_Houghton.fa
			awk '/^>/{print ">Eimeria_mitis_Houghton" ++i; next}{print}' Eimeria_mitis_Houghton.fa > header_Eimeria_mitis_Houghton.fa
		# Eimeria necatrix Houghton
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENH_00043160.fa > Eimeria_necatrix_Houghton.fa
			awk '/^>/{print ">Eimeria_necatrix_Houghton" ++i; next}{print}' Eimeria_necatrix_Houghton.fa > header_Eimeria_necatrix_Houghton.fa
		# Eimeria praecox Houghton
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENH_0031840.fa > Eimeria_praecox_Houghton.fa
			awk '/^>/{print ">Eimeria_praecox_Houghton" ++i; next}{print}' Eimeria_praecox_Houghton.fa > header_Eimeria_praecox_Houghton.fa
		# Cyclospora cayetanensis strain CHN_HEN01
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cyc_01595.fa > Cyclospora_cayetanensis.fa
			awk '/^>/{print ">Cyclospora_cayetanensis" ++i; next}{print}' Cyclospora_cayetanensis.fa > header_Cyclospora_cayetanensis.fa
		# Eimeria acervulina Houghton
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EAH_00010620.fa > Eimeria_acervulina_Houghton.fa
			awk '/^>/{print ">Eimeria_acervulina_Houghton" ++i; next}{print}' Eimeria_acervulina_Houghton.fa > header_Eimeria_acervulina_Houghton.fa
		# Cryptosporidium ubiquitum isolate 39726
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cubi_02786.fa > Cryptosporidium_ubiquitum.fa
			awk '/^>/{print ">Cryptosporidium_ubiquitum" ++i; next}{print}' Cryptosporidium_ubiquitum.fa > header_Cryptosporidium_ubiquitum.fa
		# Cryptosporidium parvum IOWA-ATCC
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CPATCC_0006870.fa > Cryptosporidium_parvum.fa
			awk '/^>/{print ">Cryptosporidium_parvum" ++i; next}{print}' Cryptosporidium_parvum.fa > header_Cryptosporidium_parvum.fa
		# Cryptosporidium tyzzeri isolate UGA55
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CTYZ_00000707.fa > Cryptosporidium_tyzzeri.fa
			awk '/^>/{print ">Cryptosporidium_tyzzeri" ++i; next}{print}' Cryptosporidium_tyzzeri.fa > header_Cryptosporidium_tyzzeri.fa
		# Cryptosporidium hominis isolate TU502_2012
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' GY17_00001902.fa > Cryptosporidium_hominis.fa
			awk '/^>/{print ">Cryptosporidium_hominis" ++i; next}{print}' Cryptosporidium_hominis.fa > header_Cryptosporidium_hominis.fa
		# Cryptosporidium meleagridis strain UKMEL1
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CmeUKMEL1_10795.fa > Cryptosporidium_meleagridis.fa
			awk '/^>/{print ">Cryptosporidium_meleagridis" ++i; next}{print}' Cryptosporidium_meleagridis.fa > header_Cryptosporidium_meleagridis.fa
		# Cryptosporidium andersoni isolate 30847
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cand_029210.fa > Cryptosporidium_andersoni.fa
			awk '/^>/{print ">Cryptosporidium_andersoni" ++i; next}{print}' Cryptosporidium_andersoni.fa > header_Cryptosporidium_andersoni.fa
		# Plasmodium falciparum 3D7
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PF3D7_1475600.fa > Plasmodium_falciparum.fa
			awk '/^>/{print ">Plasmodium_falciparum" ++i; next}{print}' Plasmodium_falciparum.fa > header_Plasmodium_falciparum.fa
		#  Plasmodium praefalciparum strain G01
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PPRFG01_1476100.fa > Plasmodium_praefalciparum.fa
			awk '/^>/{print ">Plasmodium_praefalciparum" ++i; next}{print}' Plasmodium_praefalciparum.fa > header_Plasmodium_praefalciparum.fa
Once all of the sequences where downloaded and formated, the files are then combined `cat` into a single file:
			
	cat header_Toxoplasma_gondii.fa header_Hammondia_hammondi.fa header_Neospora_caninum.fa header_Besnoitia_besnoiti.fa header_Cystoisospora_suis.fa header_Sarcocystis_neurona.fa header_Eimeria_mitis_Houghton.fa header_Eimeria_necatrix_Houghton.fa header_Eimeria_praecox_Houghton.fa header_Cyclospora_cayetanensis.fa header_Eimeria_acervulina_Houghton.fa header_Cryptosporidium_ubiquitum.fa header_Cryptosporidium_parvum.fa header_Cryptosporidium_tyzzeri.fa header_Cryptosporidium_hominis.fa header_Cryptosporidium_meleagridis.fa header_Cryptosporidium_andersoni.fa header_Plasmodium_falciparum.fa header_Plasmodium_praefalciparum.fa > TGBDP4_Phylogenetics_Tree.fa
		
A MAFFT analysis is then completed on the combined data file:
	
	mafft --auto TGBDP4_Phylogenetics_Tree.fa > Output_TGBDP4_Phylogenetics_Tree.fa
Once the MAFFT analysis is complete, the phylogenomic tree was built utilizing iqtree.

	iqtree -s Output_TGBDP4_Phylogenetics_Tree.fa -m LG -bb 1000 -pre output
`LG` is the model utilized to create the phylogenomic tree, which is based on a general matrix model. 

`-bb` is the bootstraping subset utilized to create this phylogenomic tree(1000). 

This iqtree analysis outputs a Tiled called : `output.contree`

	cat output.contree
	
The output from the cat command will then be copy/pasted into Figtree to yield the phylogenetic trees.
The Output:

		(Toxoplasma_gondii_ME491:0.047447,Hammondia_hammondi1:0.061230,(Neospora_caninum1:0.292945,(Besnoitia_besnoiti1:0.490714,((Cystoisospora_suis1:0.600032,Sarcocystis_neurona1:0.576512)43:0.077945,((((Eimeria_mitis_Houghton1:0.202561,Eimeria_acervulina_Houghton1:0.121668)99:0.098321,Cyclospora_cayetanensis1:0.460964)81:0.074075,(Eimeria_necatrix_Houghton1:0.000003,Eimeria_praecox_Houghton1:0.000003)100:0.179929)100:0.680839,(((Cryptosporidium_ubiquitum1:0.161968,(((Cryptosporidium_parvum1:0.012607,Cryptosporidium_tyzzeri1:0.018416)96:0.006900,Cryptosporidium_hominis1:0.020187)100:0.023457,Cryptosporidium_meleagridis1:0.047238)100:0.178021)100:0.286836,Cryptosporidium_andersoni1:0.630007)100:0.599680,(Plasmodium_falciparum1:0.000003,Plasmodium_praefalciparum1:0.002628)100:1.298869)100:0.678934)97:0.226232)100:0.184317)100:0.178853)100:0.424850);

# Results of the Phylogenomics Analysis of TGBDP4 predicted protein sequence within Aplicomplexians 
