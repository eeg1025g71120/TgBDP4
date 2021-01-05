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

		(Toxoplasma_gondii_ME491:0.050005,Hammondia_hammondi1:0.058172,(Neospora_caninum1:0.289418,(Besnoitia_besnoiti1:0.476065,(Cystoisospora_suis1:0.648511,(Sarcocystis_neurona1:0.522131,((((Eimeria_mitis_Houghton1:0.217398,Eimeria_acervulina_Houghton1:0.108739)100:0.100404,Cyclospora_cayetanensis1:0.450960)90:0.068706,Eimeria_necatrix_Houghton1:0.182584)100:0.721755,(((Cryptosporidium_ubiquitum1:0.155990,(((Cryptosporidium_parvum1:0.012655,Cryptosporidium_tyzzeri1:0.018422)96:0.006612,Cryptosporidium_hominis1:0.020493)100:0.023975,Cryptosporidium_meleagridis1:0.047134)100:0.183798)100:0.320993,Cryptosporidium_andersoni1:0.591399)100:0.551159,(Plasmodium_falciparum1:0.000002,Plasmodium_praefalciparum1:0.002629)100:1.331543)100:0.737380)98:0.179617)54:0.078890)100:0.206511)100:0.187370)100:0.422114);

# Results of the Phylogenomics Analysis of TGBDP4 predicted protein sequence within Aplicomplexians  
# Results from NCBI BLAST Analysis
A secondary BLAST was run comparing TgBDP4’s predicted protein sequence to similar protein sequences within NCBI’s nr database. BLASTP was utilized to identify these similar sequences. This BLAST analysis yielded the top 50 organisms that appear to have a similar protein sequence within their genome. This list was further narrowed down through the elimination of multiple samples from the same species. These sequences are:
* Neospora caninum Liverpool: XP_003881437
* Besnoitia besnoiti: XP_029216533
* Cystoisospora suis: PHJ22926
* Eimeria praecox: CDI74223
* Eimeria mitis: XP_013350930
* Eimeria necatrix: XP_013436311
* Cyclospora cayetanensis: OEH79948
* Eimeria acervuline: XP_013250528
* Cryptosporidium felis: KAF7457151
* Plasmodium Falciparum: 4NXJ_A 
* Cryptosporidium ubiquitum: XP_028875204
* Cryptosporidium meleagridis: POM84118
* Cryptosporidium parvum: >XP_628299
* Cryptosporidium hominis: PPS94707
* Cryptosporidium tyzzeri: TRY50306
* Cryptosporidium andersoni: OII75837
* Cryptosporidium muris RN66: XP_002141275
* Plasmodium reichenowi: SOV83627
* Marchantia polymorpha:PTQ48795
* Plasmodium gaboni: XP_018639804
* Plasmodium sp. DRC-Itaito: SOV25474
* Plasmodium sp. gorilla clade G1: SOS81579
* Hepatocystis sp. ex Piliocolobus tephrosceles : VWU51739
* Ceratodon purpureus: KAG0617482
* Plasmodium vivax: XP_001616024
* Plasmodium ovale curtisi: SBS97265
* Trifolium subterraneum: GAU20458
* Chara braunii: GBG83407
* Plasmodium knowlesi strain H: XP_002260507
* Elaeis guineensis: XP_010940213
* Plasmodium fragile: XP_012336669
* Thalictrum thalictroides: KAF5207984
* Plasmodium relictum: XP_028534559
* Plasmodium malariae: SBS89083
* Musa balbisiana: THU69887
* Plasmodium inui San Antonio 1:XP_008819439
* Plasmodium gallinaceum: XP_028527377
* Ananas comosus: OAY69483
* Theileria equi strain WA: XP_004828861
* Populus euphratica: XP_011040801
* Musa acuminata subsp. malaccensis: XP_018686019
* Populus deltoides: KAF9856795
* Populus trichocarpa: XP_006384849
* Eragrostis curvula: TVU10293
* Trifolium pratense: PNY12559
* Coregonus sp. 'balchen': CAB1317005
* Phoenix dactylifera: XP_008812297
* Plasmodium gonderi: XP_028545308
*  Vanilla planifolia : KAG0490016
* Digitaria exilis: CAB3457785
# Phylogenomics Analysis of TGBDP4 predicted protein sequence with Uniprot
		#Toxoplasma gondii ME49
			 awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TGME49_306460.fa > Toxoplasma_gondii.fa
			 awk '/^>/{print ">Toxoplasma_gondii_ME49" ++i; next}{print}' Toxoplasma_gondii.fa > header_Toxoplasma_gondii.fa
		# Neospora caninum Liverpool
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_003881437.fa > Neospora_caninum.fa
			awk '/^>/{print ">Neospora_caninum " ++i; next}{print}' Neospora_caninum.fa > header_Neospora_caninum.fa
		# Besnoitia besnoiti: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_029216533.fa > Besnoitia_besnoiti.fa
			awk '/^>/{print ">Besnoitia_besnoiti" ++i; next}{print}' Besnoitia_besnoiti.fa > header_Besnoitia_besnoiti.fa
		# Cystoisospora suis
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PHJ22926.fa > Cystoisospora_suis.fa
			awk '/^>/{print ">Cystoisospora_suis " ++i; next}{print}' Cystoisospora_suis.fa > header_Cystoisospora_suis.fa
		# Eimeria praecox: CDI74223
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CDI74223.fa > Eimeria_praecox.fa
			awk '/^>/{print "> Eimeria_praecox " ++i; next}{print}' Eimeria_praecox.fa > header_Eimeria_praecox.fa
		# Eimeria mitis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_013350930.fa > Eimeria_mitis.fa
			awk '/^>/{print "> Eimeria_mitis " ++i; next}{print}' Eimeria_mitis.fa > header_Eimeria_mitis.fa
		# Eimeria necatrix: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_013436311.fa > Eimeria¬_necatrix.fa
			awk '/^>/{print "> Eimeria¬_necatrix " ++i; next}{print}' Eimeria¬_necatrix.fa > header_Eimeria¬_necatrix.fa
		# Cyclospora cayetanensis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' OEH79948.fa > Cyclospora_cayetanensis.fa
			awk '/^>/{print ">Cyclospora_cayetanensis " ++i; next}{print}' Cyclospora_cayetanensis.fa > header_Cyclospora_cayetanensis.fa
		# Eimeria acervuline: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_013250528.fa > Eimeria_acervuline.fa
			awk '/^>/{print ">Eimeria_acervuline " ++i; next}{print}' Eimeria_acervuline.fa > header_Eimeria_acervuline.fa
		# Cryptosporidium felis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' KAF7457151.fa > Cryptosporidium_felis.fa
			awk '/^>/{print ">Cryptosporidium_felis " ++i; next}{print}' Cryptosporidium_felis.fa > header_Cryptosporidium_felis.fa
		# Plasmodium Falciparum: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' 4NXJ_A.fa > Plasmodium_Falciparum.fa
			awk '/^>/{print ">Plasmodium_Falciparum " ++i; next}{print}' Plasmodium_Falciparum.fa > header_Plasmodium_Falciparum.fa
		# Cryptosporidium ubiquitum:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_028875204.fa > Cryptosporidium_ubiquitum.fa
			awk '/^>/{print ">Cryptosporidium_ubiquitum " ++i; next}{print}' Cryptosporidium_ubiquitum.fa > header_Cryptosporidium_ubiquitum.fa
		# Cryptosporidium meleagridis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' POM84118.fa > Cryptosporidium_meleagridis.fa
			awk '/^>/{print ">Cryptosporidium_meleagridis " ++i; next}{print}' Cryptosporidium_meleagridis.fa > header_Cryptosporidium_meleagridis.fa
		# Cryptosporidium parvum: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_628299.fa > Cryptosporidium_parvum.fa
			awk '/^>/{print ">Cryptosporidium_parvum " ++i; next}{print}' Cryptosporidium_parvum.fa > header_Cryptosporidium_parvum.fa
		# Cryptosporidium hominis:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PPS94707.fa > Cryptosporidium_hominis.fa
			awk '/^>/{print ">Cryptosporidium_hominis " ++i; next}{print}' Cryptosporidium_hominis.fa > header_Cryptosporidium_hominis.fa
		# Cryptosporidium tyzzeri: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TRY50306.fa > Cryptosporidium¬_tyzzeri.fa
			awk '/^>/{print ">Cryptosporidium¬_tyzzeri " ++i; next}{print}' Cryptosporidium¬_tyzzeri.fa > header_Cryptosporidium¬_tyzzeri.fa
		# Cryptosporidium andersoni: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' OII75837.fa > Cryptosporidium_andersoni.fa
			awk '/^>/{print "> Cryptosporidium_andersoni" ++i; next}{print}' Cryptosporidium_andersoni.fa > header_Cryptosporidium_andersoni.fa
		# Cryptosporidium muris RN66: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_002141275.fa > Cryptosporidium_muris.fa
			awk '/^>/{print ">Cryptosporidium_muris " ++i; next}{print}' Cryptosporidium_muris.fa > header_Cryptosporidium_muris.fa
		# Plasmodium reichenowi: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SOV83627.fa > Plasmodium_reichenowi.fa
			awk '/^>/{print ">Plasmodium_reichenowi " ++i; next}{print}' Plasmodium_reichenowi.fa > header_Plasmodium_reichenowi.fa
		# Marchantia polymorpha:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PTQ48795.fa > Marchantia_polymorpha.fa
			awk '/^>/{print "> Marchantia_polymorpha " ++i; next}{print}' Marchantia_polymorpha.fa > header_Marchantia_polymorpha.fa
		# Plasmodium gaboni:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_018639804.fa > Plasmodium_gaboni.fa
			awk '/^>/{print ">Plasmodium_gaboni " ++i; next}{print}' Plasmodium_gaboni.fa > header_Plasmodium_gaboni.fa
		# Plasmodium sp. DRC-Itaito:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SOV25474.fa > Plasmodium_sp_DRC.fa
			awk '/^>/{print "> Plasmodium_sp_DRC-Itaito" ++i; next}{print}' Plasmodium_sp_DRC.fa > header_Plasmodium_sp_DRC.fa
		# Plasmodium sp. gorilla clade G1:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SOS81579.fa > Plasmodium_sp_gorilla.fa
			awk '/^>/{print ">Plasmodium_sp_gorilla " ++i; next}{print}' Plasmodium_sp_gorilla.fa > header_Plasmodium_sp_gorilla.fa
		# Hepatocystis sp. ex Piliocolobus tephrosceles : 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' VWU51739.fa > Hepatocystis_sp.fa
			awk '/^>/{print ">Hepatocystis_sp" ++i; next}{print}' Hepatocystis_sp.fa > header_Hepatocystis_sp.fa
		# Ceratodon purpureus: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' KAG0617482.fa > Ceratodon_purpureus.fa
			awk '/^>/{print ">Ceratodon_purpureus" ++i; next}{print}' Ceratodon_purpureus.fa > header_Ceratodon_purpureus.fa


		
