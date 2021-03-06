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
			awk '/^>/{print ">Neospora_caninum" ++i; next}{print}' Neospora_caninum.fa > header_Neospora_caninum.fa
		# Besnoitia besnoiti: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_029216533.fa > Besnoitia_besnoiti.fa
			awk '/^>/{print ">Besnoitia_besnoiti" ++i; next}{print}' Besnoitia_besnoiti.fa > header_Besnoitia_besnoiti.fa
		# Cystoisospora suis
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PHJ22926.fa > Cystoisospora_suis.fa
			awk '/^>/{print ">Cystoisospora_suis" ++i; next}{print}' Cystoisospora_suis.fa > header_Cystoisospora_suis.fa
		# Eimeria praecox:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CDI74223.fa > Eimeria_praecox.fa
			awk '/^>/{print ">Eimeria_praecox" ++i; next}{print}' Eimeria_praecox.fa > header_Eimeria_praecox.fa
		# Eimeria mitis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_013350930.fa > Eimeria_mitis.fa
			awk '/^>/{print ">Eimeria_mitis" ++i; next}{print}' Eimeria_mitis.fa > header_Eimeria_mitis.fa
		# Eimeria necatrix: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_013436311.fa > Eimeria_necatrix.fa
			awk '/^>/{print ">Eimeria_necatrix" ++i; next}{print}' Eimeria_necatrix.fa > header_Eimeria_necatrix.fa
		# Cyclospora cayetanensis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' OEH79948.fa > Cyclospora_cayetanensis.fa
			awk '/^>/{print ">Cyclospora_cayetanensis" ++i; next}{print}' Cyclospora_cayetanensis.fa > header_Cyclospora_cayetanensis.fa
		# Eimeria acervuline: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_013250528.fa > Eimeria_acervuline.fa
			awk '/^>/{print ">Eimeria_acervuline" ++i; next}{print}' Eimeria_acervuline.fa > header_Eimeria_acervuline.fa
		# Cryptosporidium felis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' KAF7457151.fa > Cryptosporidium_felis.fa
			awk '/^>/{print ">Cryptosporidium_felis" ++i; next}{print}' Cryptosporidium_felis.fa > header_Cryptosporidium_felis.fa
		# Plasmodium Falciparum: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' 4NXJ_A.fa > Plasmodium_Falciparum.fa
			awk '/^>/{print ">Plasmodium_Falciparum " ++i; next}{print}' Plasmodium_Falciparum.fa > header_Plasmodium_Falciparum.fa
		# Cryptosporidium ubiquitum:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_028875204.fa > Cryptosporidium_ubiquitum.fa
			awk '/^>/{print ">Cryptosporidium_ubiquitum" ++i; next}{print}' Cryptosporidium_ubiquitum.fa > header_Cryptosporidium_ubiquitum.fa
		# Cryptosporidium meleagridis: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' POM84118.fa > Cryptosporidium_meleagridis.fa
			awk '/^>/{print ">Cryptosporidium_meleagridis" ++i; next}{print}' Cryptosporidium_meleagridis.fa > header_Cryptosporidium_meleagridis.fa
		# Cryptosporidium parvum: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_628299.fa > Cryptosporidium_parvum.fa
			awk '/^>/{print ">Cryptosporidium_parvum" ++i; next}{print}' Cryptosporidium_parvum.fa > header_Cryptosporidium_parvum.fa
		# Cryptosporidium hominis:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PPS94707.fa > Cryptosporidium_hominis.fa
			awk '/^>/{print ">Cryptosporidium_hominis" ++i; next}{print}' Cryptosporidium_hominis.fa > header_Cryptosporidium_hominis.fa
		# Cryptosporidium tyzzeri: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TRY50306.fa > Cryptosporidium_tyzzeri.fa
			awk '/^>/{print ">Cryptosporidium_tyzzeri" ++i; next}{print}' Cryptosporidium_tyzzeri.fa > header_Cryptosporidium_tyzzeri.fa
		# Cryptosporidium andersoni: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' OII75837.fa > Cryptosporidium_andersoni.fa
			awk '/^>/{print ">Cryptosporidium_andersoni" ++i; next}{print}' Cryptosporidium_andersoni.fa > header_Cryptosporidium_andersoni.fa
		# Cryptosporidium muris RN66: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_002141275.fa > Cryptosporidium_muris.fa
			awk '/^>/{print ">Cryptosporidium_muris" ++i; next}{print}' Cryptosporidium_muris.fa > header_Cryptosporidium_muris.fa
		# Plasmodium reichenowi: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SOV83627.fa > Plasmodium_reichenowi.fa
			awk '/^>/{print ">Plasmodium_reichenowi" ++i; next}{print}' Plasmodium_reichenowi.fa > header_Plasmodium_reichenowi.fa
		# Marchantia polymorpha:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PTQ48795.fa > Marchantia_polymorpha.fa
			awk '/^>/{print ">Marchantia_polymorpha" ++i; next}{print}' Marchantia_polymorpha.fa > header_Marchantia_polymorpha.fa
		# Plasmodium gaboni:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_018639804.fa > Plasmodium_gaboni.fa
			awk '/^>/{print ">Plasmodium_gaboni" ++i; next}{print}' Plasmodium_gaboni.fa > header_Plasmodium_gaboni.fa
		# Plasmodium sp. DRC-Itaito:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SOV25474.fa > Plasmodium_sp_DRC.fa
			awk '/^>/{print ">Plasmodium_sp_DRC-Itaito" ++i; next}{print}' Plasmodium_sp_DRC.fa > header_Plasmodium_sp_DRC.fa
		# Plasmodium sp. gorilla clade G1:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SOS81579.fa > Plasmodium_sp_gorilla.fa
			awk '/^>/{print ">Plasmodium_sp_gorilla" ++i; next}{print}' Plasmodium_sp_gorilla.fa > header_Plasmodium_sp_gorilla.fa
		# Hepatocystis sp. ex Piliocolobus tephrosceles : 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' VWU51739.fa > Hepatocystis_sp.fa
			awk '/^>/{print ">Hepatocystis_sp" ++i; next}{print}' Hepatocystis_sp.fa > header_Hepatocystis_sp.fa
		# Ceratodon purpureus: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' KAG0617482.fa > Ceratodon_purpureus.fa
			awk '/^>/{print ">Ceratodon_purpureus" ++i; next}{print}' Ceratodon_purpureus.fa > header_Ceratodon_purpureus.fa
		# Plasmodium vivax: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_001616024.fa > Plasmodium_vivax.fa
			awk '/^>/{print ">Plasmodium_vivax" ++i; next}{print}' Plasmodium_vivax.fa > header_Plasmodium_vivax.fa
		# Plasmodium ovale curtisi:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SBS97265.fa > Plasmodium_ovale.fa
			awk '/^>/{print ">Plasmodium_ovale" ++i; next}{print}' Plasmodium_ovale.fa > header_Plasmodium_ovale.fa
		# Trifolium subterraneum: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' GAU20458.fa > Trifolium_subterraneum.fa
			awk '/^>/{print ">Trifolium_subterraneum" ++i; next}{print}' Trifolium_subterraneum.fa > header_Trifolium_subterraneum.fa
		# Chara braunii: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' GBG83407.fa > Chara_braunii.fa
			awk '/^>/{print ">Chara_braunii" ++i; next}{print}' Chara_braunii.fa > header_Chara_braunii.fa
		# Plasmodium knowlesi strain H: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_002260507.fa > Plasmodium_knowlesi.fa
			awk '/^>/{print ">Plasmodium_knowlesi" ++i; next}{print}' Plasmodium_knowlesi.fa > header_Plasmodium_knowlesi.fa
		# Elaeis guineensis:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_010940213.fa > Elaeis_guineensis.fa
			awk '/^>/{print ">Elaeis_guineensis" ++i; next}{print}' Elaeis_guineensis.fa > header_Elaeis_guineensis.fa
		# Plasmodium fragile: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_012336669.fa > Plasmodium_fragile.fa
			awk '/^>/{print ">Plasmodium_fragile" ++i; next}{print}' Plasmodium_fragile.fa > header_Plasmodium_fragile.fa
		# Thalictrum thalictroides
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' KAF5207984.fa > Thalictrum_thalictroides.fa
			awk '/^>/{print "> Thalictrum_thalictroides" ++i; next}{print}' Thalictrum_thalictroides.fa  > header_Thalictrum_thalictroides.fa
		# Plasmodium relictum: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_028534559.fa > Plasmodium_relictum.fa
			awk '/^>/{print ">Plasmodium_relictum" ++i; next}{print}' Plasmodium_relictum.fa  > header_Plasmodium_relictum.fa
		# Plasmodium malariae:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SBS89083.fa > Plasmodium_malariae.fa
			awk '/^>/{print ">Plasmodium_malariae" ++i; next}{print}' Plasmodium_malariae.fa  > header_Plasmodium_malariae.fa
		# Musa balbisiana:
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' THU69887.fa > Musa_balbisiana.fa
			awk '/^>/{print ">Musa_balbisiana" ++i; next}{print}' Musa_balbisiana.fa  > header_Musa_balbisiana.fa
		# Plasmodium inui San Antonio 1
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_008819439.fa > Plasmodium_inui.fa
			awk '/^>/{print ">Plasmodium_inui" ++i; next}{print}' Plasmodium_inui.fa  > header_Plasmodium_inui.fa
		# Plasmodium gallinaceum
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_028527377.fa > Plasmodium_gallinaceum.fa
			awk '/^>/{print ">Plasmodium_gallinaceum" ++i; next}{print}' Plasmodium_gallinaceum.fa > header_Plasmodium_gallinaceum.fa
		# Ananas comosus
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' OAY69483.fa > Ananas_comosus.fa
			awk '/^>/{print ">Ananas_comosus" ++i; next}{print}' Ananas_comosus.fa > header_Ananas_comosus.fa
		# Theileria equi strain WA
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_004828861.fa > Theileria_equi.fa
			awk '/^>/{print ">Theileria_equi" ++i; next}{print}' Theileria_equi.fa > header_Theileria_equi.fa
		# Populus euphratica 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_011040801.fa > Populus_euphratica.fa
			awk '/^>/{print ">Populus_euphratica" ++i; next}{print}' Populus_euphratica.fa > header_Populus_euphratica.fa
		# Musa acuminata subsp. malaccensis
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_018686019.fa > Musa_acuminata.fa
			awk '/^>/{print ">Musa_acuminata" ++i; next}{print}' Musa_acuminata.fa > header_Musa_acuminata.fa
		# Populus deltoides
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' KAF9856795.fa > Populus_deltoides.fa
			awk '/^>/{print ">Populus_deltoides" ++i; next}{print}' Populus_deltoides.fa > header_Populus_deltoides.fa
		# Populus trichocarpa
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_006384849.fa > Populus_trichocarpa.fa
			awk '/^>/{print ">Populus_trichocarpa" ++i; next}{print}' Populus_trichocarpa.fa > header_Populus_trichocarpa.fa
		# Eragrostis curvula
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TVU10293.fa > Eragrostis_curvula.fa
			awk '/^>/{print ">Eragrostis_curvula" ++i; next}{print}' Eragrostis_curvula.fa > header_Eragrostis_curvula.fa
		# Trifolium pratense: 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PNY12559.fa > Trifolium_pratense.fa
			awk '/^>/{print ">Trifolium_pratense" ++i; next}{print}' Trifolium_pratense.fa > header_Trifolium_pratense.fa
		# Coregonus sp. 'balchen'
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CAB1317005.fa > Coregonus_sp.fa
			awk '/^>/{print ">Coregonus_sp" ++i; next}{print}' Coregonus_sp.fa > header_Coregonus_sp.fa
		# Phoenix dactylifera
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_008812297.fa > Phoenix_dactylifera.fa
			awk '/^>/{print ">Phoenix_dactylifera" ++i; next}{print}' Phoenix_dactylifera.fa > header_Phoenix_dactylifera.fa
		# Plasmodium gonderi
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' XP_028545308.fa > Plasmodium_gonderi.fa
			awk '/^>/{print ">Plasmodium_gonderi" ++i; next}{print}' Plasmodium_gonderi.fa > header_Plasmodium_gonderi.fa
		# Vanilla planifolia 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' KAG0490016.fa > Vanilla_planifolia.fa
			awk '/^>/{print ">Vanilla_planifolia" ++i; next}{print}' Vanilla_planifolia.fa > header_Vanilla_planifolia.fa
		# Digitaria exilis
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CAB3457785.fa > Digitaria_exilis.fa
			awk '/^>/{print ">Digitaria_exilis" ++i; next}{print}' Digitaria_exilis.fa > header_Digitaria_exilis.fa
Once all of the sequences where downloaded and formated, the files are then combined cat into a single file:
		
		cat header_Toxoplasma_gondii.fa header_Neospora_caninum.fa header_Besnoitia_besnoiti.fa header_Cystoisospora_suis.fa header_Eimeria_praecox.fa header_Eimeria_mitis.fa header_Eimeria_necatrix.fa header_Cyclospora_cayetanensis.fa header_Eimeria_acervuline.fa header_Cryptosporidium_felis.fa header_Plasmodium_Falciparum.fa header_Cryptosporidium_ubiquitum.fa header_Cryptosporidium_meleagridis.fa header_Cryptosporidium_parvum.fa header_Cryptosporidium_hominis.fa header_Cryptosporidium_tyzzeri.fa header_Cryptosporidium_andersoni.fa header_Cryptosporidium_muris.fa header_Plasmodium_reichenowi.fa header_Marchantia_polymorpha.fa header_Plasmodium_gaboni.fa header_Plasmodium_sp_DRC.fa header_Plasmodium_sp_gorilla.fa header_Hepatocystis_sp.fa header_Ceratodon_purpureus.fa header_Plasmodium_vivax.fa header_Plasmodium_ovale.fa header_Trifolium_subterraneum.fa header_Chara_braunii.fa header_Plasmodium_knowlesi.fa header_Elaeis_guineensis.fa header_Plasmodium_fragile.fa header_Thalictrum_thalictroides.fa header_Plasmodium_relictum.fa header_Plasmodium_malariae.fa header_Musa_balbisiana.fa header_Plasmodium_inui.fa header_Plasmodium_gallinaceum.fa header_Ananas_comosus.fa header_Theileria_equi.fa header_Populus_euphratica.fa header_Musa_acuminata.fa header_Populus_deltoides.fa header_Populus_trichocarpa.fa header_Eragrostis_curvula.fa header_Trifolium_pratense.fa header_Coregonus_sp.fa header_Phoenix_dactylifera.fa header_Plasmodium_gonderi.fa header_Vanilla_planifolia.fa header_Digitaria_exilis.fa > TGBDP4_NCBI_Phylogenetics_Tree.fa

A MAFFT analysis is then completed on the combined data file: 

	mafft --auto TGBDP4_NCBI_Phylogenetics_Tree.fa > Output_TGBDP4_NCBI_Phylogenetics_Tree.fa

Once the MAFFT analysis is complete, the phylogenomic tree was built utilizing iqtree.

	iqtree -s Output_TGBDP4_NCBI_Phylogenetics_Tree.fa -m LG -bb 1000 -pre output
This iqtree analysis outputs a Tiled called : `output.contree`

	cat output.contree
	
The output from the cat command will then be copy/pasted into Figtree to yield the phylogenetic trees.
The Output:
	
	(Toxoplasma_gondii_ME491:0.258297,(Neospora_caninum1:0.143118,Cystoisospora_suis1:0.313194)92:0.106610,(Besnoitia_besnoiti1:0.340847,((((Eimeria_praecox1:0.000002,(Eimeria_mitis1:0.000002,Eimeria_acervuline1:0.074491)100:0.049649)99:0.043755,Cyclospora_cayetanensis1:0.115521)99:0.111957,Eimeria_necatrix1:0.070021)100:0.308403,((((Eimeria_mitis2:0.000003,Eimeria_necatrix2:0.281481)80:0.016224,Eimeria_acervuline2:0.045379)100:0.340717,Cyclospora_cayetanensis2:0.000002)100:0.760381,((((Cryptosporidium_felis1:0.084187,(Cryptosporidium_ubiquitum1:0.000002,(Cryptosporidium_meleagridis1:0.009591,((Cryptosporidium_parvum1:0.000000,Cryptosporidium_tyzzeri1:0.000000):0.000002,Cryptosporidium_hominis1:0.000002)100:0.009742)99:0.029492)98:0.014971)96:0.092512,(Cryptosporidium_andersoni1:0.009782,Cryptosporidium_muris1:0.000002)100:0.109065)100:0.604914,((((Marchantia_polymorpha1:0.077114,Ceratodon_purpureus1:0.305611)98:0.107755,((((Trifolium_subterraneum1:0.000002,Trifolium_pratense1:0.019131)100:0.210519,(Musa_balbisiana1:0.015486,Musa_acuminata1:0.026083)100:0.351090)97:0.110763,Thalictrum_thalictroides1:0.274367)98:0.093828,(((Elaeis_guineensis1:0.017427,Phoenix_dactylifera1:0.061622)93:0.025104,((Ananas_comosus1:0.106175,Vanilla_planifolia1:0.153688)76:0.017652,(Eragrostis_curvula1:0.030986,Digitaria_exilis1:0.280184)98:0.162235)53:0.045414)55:0.088247,(Populus_euphratica1:0.006916,(Populus_deltoides1:0.006875,Populus_trichocarpa1:0.000002)87:0.006834)100:0.207880)78:0.169309)100:0.279607)97:0.124846,Chara_braunii1:0.188892)99:0.276637,Coregonus_sp1:0.704925)86:0.186629)86:0.272297,(((Cryptosporidium_felis2:0.109477,(Cryptosporidium_ubiquitum2:0.085285,(Cryptosporidium_meleagridis2:0.010152,(Cryptosporidium_parvum2:0.000002,(Cryptosporidium_hominis2:0.000002,Cryptosporidium_tyzzeri2:0.012595)35:0.000002)30:0.000002)95:0.067205)92:0.069059)88:0.165273,(Cryptosporidium_andersoni2:0.000003,Cryptosporidium_muris2:0.030957)99:0.188337)99:0.792215,((((((Plasmodium_Falciparum:0.000002,Plasmodium_sp_gorilla1:0.000003)97:0.015024,Plasmodium_reichenowi1:0.000002)87:0.021046,((((Hepatocystis_sp1:0.218234,Plasmodium_malariae1:0.004672)82:0.047435,((((Plasmodium_vivax1:0.009313,Plasmodium_inui1:0.038405)85:0.009568,Plasmodium_fragile1:0.028413)43:0.000003,Plasmodium_knowlesi1:0.048914)90:0.027777,Plasmodium_gonderi1:0.012070)99:0.095661)71:0.024488,Plasmodium_ovale1:0.047680)78:0.017973,(Plasmodium_relictum1:0.052895,Plasmodium_gallinaceum1:0.025699)90:0.090698)80:0.194702)54:0.031813,(Plasmodium_gaboni1:0.000002,Plasmodium_sp_DRC-Itaito1:0.000002)69:0.000002)72:0.365455,((Plasmodium_reichenowi2:0.018116,(Plasmodium_gaboni2:0.000002,Plasmodium_sp_DRC-Itaito2:0.000002)100:0.117953)84:0.084006,((((((Hepatocystis_sp2:0.700580,(Plasmodium_vivax2:0.000002,Plasmodium_inui2:0.062767)42:0.000002)43:0.011897,Plasmodium_knowlesi2:0.024823)41:0.008576,Plasmodium_fragile2:0.007281)42:0.039867,Plasmodium_gonderi2:0.042030)85:0.081688,(Plasmodium_relictum2:0.000002,Plasmodium_gallinaceum2:0.094589)100:0.187157)47:0.010960,Plasmodium_ovale2:0.101369)96:0.189296)100:2.005728)83:0.642211,Theileria_equi1:0.375639)75:0.244679)80:0.283121)77:0.048037)87:0.142786)87:0.131648)82:0.082490);


# Phylogenomics Analysis of TGBDP4 predicted protein sequence with Extended data
		#Toxoplasma gondii ME49
			 awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TGME49_306460.fa > Toxoplasma_gondii.fa
			 awk '/^>/{print ">Toxoplasma_gondii_ME49" ++i; next}{print}' Toxoplasma_gondii.fa > header_Toxoplasma_gondii.fa
		# Acanthamoeba castellanii str. Neff
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ACA1_215690.fa > Acanthamoeba_castellanii.fa
			awk '/^>/{print "> Acanthamoeba_castellanii" ++i; next}{print}' Acanthamoeba_castellanii.fa > header_Acanthamoeba_castellanii.fa
		# Phanerochaete chrysosporium
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' AGR57_3965.fa > Phanerochaete_chrysosporium.fa
			awk '/^>/{print "> Phanerochaete_chrysosporium" ++i; next}{print}' Phanerochaete_chrysosporium.fa > header_Phanerochaete_chrysosporium.fa
		# Allomyces macrogynus
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' AMAG_05489.fa > Allomyces_macrogynus.fa
			awk '/^>/{print "> Allomyces_macrogynus" ++i; next}{print}' Allomyces macrogynus.fa > header_Allomyces_macrogynus.fa
		# Batrachochytrium dendrobatidis
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BDEG_07408.fa > Batrachochytrium_dendrobatidis.fa
			awk '/^>/{print "> Batrachochytrium_dendrobatidis" ++i; next}{print}' Batrachochytrium_dendrobatidis.fa > header_ Batrachochytrium_dendrobatidis.fa
		# Besnoitia besnoiti strain Bb-Ger1 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BESB_018420.fa > Besnoitia_besnoiti.fa
			awk '/^>/{print "> Besnoitia_besnoiti" ++i; next}{print}' Besnoitia_besnoiti.fa > header_Besnoitia_besnoiti.fa
		# Bos taurus
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BOS_21376.fa > Bos_taurus.fa
			awk '/^>/{print "> Bos_taurus" ++i; next}{print}' Bos_taurus.fa > header_Bos_taurus.fa
		# Cryptosporidium hominis UdeA01
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CHUDEA7_1130.fa > Cryptosporidium_hominis.fa
			awk '/^>/{print "> Cryptosporidium_hominis " ++i; next}{print}' Cryptosporidium_hominis.fa > header_Cryptosporidium_hominis.fa
		# Cryptococcus neoformans var. grubii KN99 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CKF44_00375.fa > Cryptococcus_neoformans.fa
			awk '/^>/{print "> Cryptococcus_neoformans" ++i; next}{print}' Cryptococcus_neoformans.fa > header_Cryptococcus_neoformans.fa
		#  Cryptosporidium muris RN66 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CMU_033110.fa > Cryptosporidium_muris.fa
			awk '/^>/{print "> Cryptosporidium_muris" ++i; next}{print}' Cryptosporidium_muris.fa > header_Cryptosporidium_muris.fa
		#Cryptosporidium parvum IOWA-ATCC
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CPATCC_0006870.fa > Cryptosporidium_parvum.fa
			awk '/^>/{print ">Cryptosporidium_parvum" ++i; next}{print}' Cryptosporidium_parvum.fa > header_Cryptosporidium_parvum.fa
		# Cystoisospora suis strain Wien I 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CSUI_003223.fa > Cystoisospora_suis.fa
			awk '/^>/{print "> Cystoisospora_suis " ++i; next}{print}' Cystoisospora_suis.fa > header_Cystoisospora_suis.fa
		#Cryptosporidium tyzzeri isolate UGA55
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CTYZ_00000707.fa > Cryptosporidium_tyzzeri.fa
			awk '/^>/{print "> Cryptosporidium_tyzzeri" ++i; next}{print}' Cryptosporidium_tyzzeri.fa > header_Cryptosporidium_tyzzeri.fa
		# Cryptosporidium meleagridis strain UKMEL1 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CmeUKMEL1_10795.fa > Cryptosporidium_meleagridis.fa
			awk '/^>/{print "> Cryptosporidium_meleagridis" ++i; next}{print}' Cryptosporidium_meleagridis.fa > header_Cryptosporidium_meleagridis.fa
		# Eimeria acervulina Houghton
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EAH_00010620.fa > Eimeria_acervulina.fa
			awk '/^>/{print "> Eimeria_acervulina" ++i; next}{print}' Eimeria_acervulina.fa > header_Eimeria_acervulina.fa
		# Eimeria mitis Houghton
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EMH_0040220.fa > Eimeria_mitis.fa
			awk '/^>/{print "> Eimeria_mitis" ++i; next}{print}' Eimeria_mitis.fa > header_Eimeria_mitis.fa
		# Eimeria necatrix Houghton 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENH_00043160.fa > Eimeria_necatrix.fa
			awk '/^>/{print "> Eimeria_necatrix" ++i; next}{print}' Eimeria_necatrix.fa > header_Eimeria_necatrix.fa
		# Homo sapiens REF 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENSG0000009607019.fa > Homo_sapiens.fa
			awk '/^>/{print "> Homo_sapiens" ++i; next}{print}' Homo_sapiens.fa > header_Homo_sapiens.fa
		# Mus musculus C57BL6J 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENSMUSG0000006395215.fa > Mus_musculus.fa
			 awk '/^>/{print "> Mus_musculus" ++i; next}{print}' Mus_musculus.fa > header_Mus_musculus.fa
		# Hammondia hammondi strain H.H.34 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' HHA_306460.fa > Hammondia_hammondi.fa
			awk '/^>/{print "> Hammondia_hammondi" ++i; next}{print}' Hammondia_hammondi.fa > header_Hammondia_hammondi.fa
		# Mucor circinelloides 1006PhL 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' HMPREF1544_04515.fa > Mucor_circinelloides.fa
			awk '/^>/{print "> Mucor_circinelloides" ++i; next}{print}' Mucor_circinelloides.fa > header_Mucor_circinelloides.fa
		# Kwoniella bestiolae CBS 10118 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' I302_05158.fa > Kwoniella_bestiolae.fa
			awk '/^>/{print "> Kwoniella_bestiolae" ++i; next}{print}' Kwoniella_bestiolae.fa > header_Kwoniella_bestiolae.fa
		# Cryptococcus gattii VGIV IND107 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' I308_01860.fa > Cryptococcus_gattii.fa
			awk '/^>/{print ">Cryptococcus_gattii" ++i; next}{print}' Cryptococcus_gattii.fa > header_Cryptococcus_gattii.fa
		# Kwoniella heveanensis CBS 569 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' I317_06430.fa > Kwoniella_heveanensis.fa
			awk '/^>/{print "> Kwoniella_heveanensis" ++i; next}{print}' Kwoniella_heveanensis.fa > header_Kwoniella_heveanensis.fa
		# Cyclospora cayetanensis isolate NF1_C8
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' LOC34618577.fa > Cyclospora_cayetanensis.fa
			awk '/^>/{print "> Cyclospora_cayetanensis" ++i; next}{print}' Cyclospora_cayetanensis.fa > header_Cyclospora_cayetanensis.fa
		# Macaca mulatta isolate 17573
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' MACM_27154.fa > Macaca_mulatta.fa
			awk '/^>/{print "> Macaca_mulatta" ++i; next}{print}' Macaca_mulatta.fa > header_Macaca_mulatta.fa
		# Neospora caninum Liverpool 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NCLIV_044660.fa > Neospora_caninum.fa
			awk '/^>/{print ">Neospora_caninum " ++i; next}{print}' Neospora_caninum.fa > header_Neospora_caninum.fa
		# Phycomyces blakesleeanus
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PHYBL_128978.fa > Phycomyces_blakesleeanus.fa
			awk '/^>/{print ">Phycomyces_blakesleeanus" ++i; next}{print}' Phycomyces_blakesleeanus.fa > header_Phycomyces_blakesleeanus.fa
		# Globisporangium iwayamae DAOM BR242034
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' PIW_G002910.fa > Globisporangium_iwayamae.fa
			awk '/^>/{print ">Globisporangium_iwayamae" ++i; next}{print}' Globisporangium_iwayamae.fa > header_Globisporangium_iwayamae.fa
		# Rhizopus delemar RA 99-880 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' RO3G_02044.fa > Rhizopus_delemar.fa
			awk '/^>/{print ">Rhizopus_delemar" ++i; next}{print}' Rhizopus_delemar.fa > header_Rhizopus_delemar.fa
		# Sarcocystis neurona SN3
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SN3_01400140.fa > Sarcocystis_neurona.fa
			awk '/^>/{print ">Sarcocystis_neurona" ++i; next}{print}' Sarcocystis_neurona.fa > header_Sarcocystis_neurona.fa
		# Spizellomyces punctatus DAOM BR117
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SPPG_08745.fa > Spizellomyces_punctatus.fa
			awk '/^>/{print ">Spizellomyces_punctatus" ++i; next}{print}' Spizellomyces_punctatus.fa > header_Spizellomyces_punctatus.fa
		# Ustilago maydis 521
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' UMAG_06029.fa > Ustilago_maydis.fa
			awk '/^>/{print ">Ustilago_maydis" ++i; next}{print}' Ustilago_maydis.fa > header_Ustilago_maydis.fa
		# Vitrella brassicaformis CCMP3155
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Vbra_14604.fa > Vitrella_brassicaformis.fa
			awk '/^>/{print ">Vitrella_brassicaformis" ++i; next}{print}' Vitrella_brassicaformis.fa > header_Vitrella_brassicaformis.fa
		# Cryptosporidium andersoni isolate 30847
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cand_029210.fa > Cryptosporidium_andersoni.fa
			awk '/^>/{print ">Cryptosporidium_andersoni" ++i; next}{print}' Cryptosporidium_andersoni.fa > header_Cryptosporidium_andersoni.fa
		# Cryptosporidium ubiquitum isolate 39726
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cubi_02786.fa > Cryptosporidium_ubiquitum.fa
			awk '/^>/{print ">Cryptosporidium_ubiquitum" ++i; next}{print}' Cryptosporidium_ubiquitum.fa > header_Cryptosporidium_ubiquitum.fa
		# Sporisorium reilianum SRZ2 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' sr16504.fa > Sporisorium_reilianum.fa
			awk '/^>/{print ">Sporisorium_reilianum" ++i; next}{print}' Sporisorium_reilianum.fa > header_Sporisorium_reilianum.fa
		# Plasmodium falciparum 3D7 
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  PF3D7_1475600.fa > Plasmodium_falciparum.fa
			awk '/^>/{print ">Plasmodium_falciparum" ++i; next}{print}' Plasmodium_falciparum.fa > header_Plasmodium_falciparum.fa
		#Plasmodium reichenowi G01
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  PRG01_1475400.fa > Plasmodium_reichenowi.fa
			awk '/^>/{print ">Plasmodium_reichenowi" ++i; next}{print}' Plasmodium_reichenowi.fa > header_Plasmodium_reichenowi.fa
		#Plasmodium gaboni strain SY75
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  PGSY75_1475600.fa > Plasmodium_gaboni.fa
			awk '/^>/{print ">Plasmodium_gaboni " ++i; next}{print}' Plasmodium_gaboni.fa > header_Plasmodium_gaboni.fa
		#Plasmodium vivax Sal-1
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  PVX_118655.fa > Plasmodium_vivax.fa
			awk '/^>/{print ">Plasmodium_vivax" ++i; next}{print}' Plasmodium_vivax.fa > header_Plasmodium_vivax.fa
		#Plasmodium ovale curtisi GH01
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  PocGH01_12078300.fa > Plasmodium_ovale.fa
			awk '/^>/{print ">Plasmodium_ovale" ++i; next}{print}' Plasmodium_ovale.fa > header_Plasmodium_ovale.fa
		#Plasmodium knowlesi strain H
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  PKNH_1245900.fa > Plasmodium_knowlesi.fa
			awk '/^>/{print ">Plasmodium_knowlesi" ++i; next}{print}' Plasmodium_knowlesi.fa > header_Plasmodium_knowlesi.fa
		#Plasmodium malariae UG01
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  PmUG01_12080500.fa > Plasmodium_malariae.fa
			awk '/^>/{print ">Plasmodium_malariae " ++i; next}{print}' Plasmodium_malariae.fa > header_Plasmodium_malariae.fa
	cat header_Toxoplasma_gondii.fa header_Acanthamoeba_castellanii.fa header_Acanthamoeba_castellanii.fa header_Phanerochaete_chrysosporium.fa header_Allomyces_macrogynus.fa header_ Batrachochytrium_dendrobatidis.fa header_Besnoitia_besnoiti.fa header_Bos_taurus.fa header_Cryptosporidium_hominis.fa header_Cryptococcus_neoformans.fa header_Cryptosporidium_muris.fa header_Cryptosporidium_parvum.fa header_Cystoisospora_suis.fa header_Cryptosporidium_tyzzeri.fa header_Cryptosporidium_meleagridis.fa header_Eimeria_acervulina.fa header_Eimeria_mitis.fa header_Eimeria_necatrix.fa header_Homo_sapiens.fa header_Mus_musculus.fa header_Hammondia_hammondi.fa header_Mucor_circinelloides.fa header_Kwoniella_bestiolae.fa header_Cryptococcus_gattii.fa header_Kwoniella_heveanensis.fa header_Cyclospora_cayetanensis.fa header_Macaca_mulatta.fa header_Neospora_caninum.fa header_Phycomyces_blakesleeanus.fa header_Globisporangium_iwayamae.fa header_Rhizopus_delemar.fa header_Sarcocystis_neurona.fa header_Spizellomyces_punctatus.fa header_Ustilago_maydis.fa header_Vitrella_brassicaformis.fa header_Cryptosporidium_andersoni.fa header_Cryptosporidium_ubiquitum.fa header_Sporisorium_reilianum.fa header_Plasmodium_falciparum.fa header_Plasmodium_reichenowi.fa header_Plasmodium_gaboni.fa header_Plasmodium_vivax.fa header_Plasmodium_ovale.fa header_Plasmodium_knowlesi.fa header_Plasmodium_malariae.fa > TGBDP6_Phylogenetics_Trees_4.fa

	mafft --auto TGBDP6_Phylogenetics_Trees_4.fa > Output_TGBDP6_Phylogenetics_Tree_4.fa
	iqtree -s Output_TGBDP6_Phylogenetics_Tree_4.fa -m LG -bb 1000 -pre output
	cat output.contree

		(Toxoplasma_gondii_ME491:0.048237,((((Acanthamoeba_castellanii1:1.424495,((((((((((Phanerochaete_chrysosporium1:0.585166,Kwoniella_heveanensis1:0.482312)99:0.151661,Ustilago_maydis1:0.542989)100:0.230271,Phycomyces_blakesleeanus1:0.659800)65:0.087964,Batrachochytrium_dendrobatidis1:0.736081)98:0.165198,((Bos_taurus1:0.036371,(Homo_sapiens1:0.007325,Macaca_mulatta1:0.003477)100:0.015829)96:0.036956,Mus_musculus1:0.023660)100:0.981044)96:0.179396,Allomyces_macrogynus1:0.932369)100:0.604598,(((Mucor_circinelloides1:0.251343,Rhizopus_delemar1:0.272340)100:0.741317,Spizellomyces_punctatus1:1.039762)97:0.233074,Sporisorium_reilianum1:1.323095)75:0.171247)26:0.091621,(((Cryptococcus_neoformans1:0.056369,Cryptococcus_gattii1:0.042179)100:0.337033,Kwoniella_bestiolae1:0.125625)100:1.148825,Globisporangium_iwayamae1:1.275564)36:0.087171)31:0.068453,(((Plasmodium_falciparum1:0.030968,Plasmodium_reichenowi1:0.030172)100:0.052205,Plasmodium_gaboni:0.107350)100:0.342610,(((Plasmodium_vivax1:0.082937,Plasmodium_knowlesi1:0.103766)100:0.398339,Plasmodium_ovale1:0.298334)67:0.066485,Plasmodium_malariae:0.303978)89:0.260128)100:1.227908)40:0.057370,((((Cryptosporidium_hominis:0.016954,(Cryptosporidium_parvum1:0.012797,Cryptosporidium_tyzzeri1:0.018542)100:0.009237)100:0.194296,Cryptosporidium_ubiquitum1:0.162127)100:0.279265,(Cryptosporidium_muris1:0.018275,Cryptosporidium_andersoni1:0.014128)100:0.623483)100:0.559978,Vitrella_brassicaformis1:1.147133)91:0.327593)69:0.177330)100:0.500157,(Cystoisospora_suis:0.609537,((((Eimeria_acervulina1:0.120587,Eimeria_mitis1:0.196804)100:0.095205,Cyclospora_cayetanensis1:0.289500)99:0.075362,Eimeria_necatrix1:0.181139)100:0.719916,Sarcocystis_neurona1:0.588182)96:0.228507)57:0.066309)100:0.222093,Besnoitia_besnoiti1:0.514156)100:0.140504,Neospora_caninum:0.307208)100:0.418304,Hammondia_hammondi1:0.061148);
