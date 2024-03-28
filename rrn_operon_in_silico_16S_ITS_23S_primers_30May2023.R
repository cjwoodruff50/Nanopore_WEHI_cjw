# Script to hold rrn operon in-silico primers for segmenting 16S-ITS-23S regions.
# Extended 28 August 2020 to give an 8F primer
# Extended 04 August 2021 to include Zhengming Zhang's Zregions for the 23S lsu
#       16S Variable Regions:                  23S Variable Regions 
#                                                Z0:  0 - 83
#         V1: 69 – 99 (1-89)                     Z1: 84 - 135
#         V2: 137 – 242 (92-336)                 Z2: 136 - 382
#         V3: 433 – 497 (336-516)                Z3: 383 - 727
#         V4: 576 – 682 (562-784)                Z4: 728 - 985
#         V5: 822 – 879 (785-909)                Z5: 986 - 1587
#         V6: 986 – 1043  (961-1062)             Z6: 1588 - 1722
#         V7: 1117 – 1173 (1100-1219)            Z7: 1723 - 1820
#         V8: 1243 – 1294 (1219-1390)            Z8: 1821 - 2143
#         V9: 1435 - 1465 (1390-1492)            Z9: 2144 - 2488
#
# Update 22 April 2020 from version of 17 July 2019 by adding primer_sets 5,6,7,8 and renumbering 
# the following primer sets.  This renumbering means code changes elsewhere.                   [cjw]
# Update 30 August 2021 to include Z0, defining boundaries in terms of already-entered primer sets.
# Update 29 May 2023 to remove primers 15, 16 and 17 from primer_sets 4 and 40 and add them to a new primersset 41.

primer_type = c("16S_5p","16S_3p","23S_5p","23S_3p","V1_5p","V1_3p","V2_5p","V2_3p","V3_5p",
                "V3_3p","V4_5p","V4_3p","V5_5p","V5_3p","V6_5p","V6_3p","V7_5p","V7_3p","V9_5p","V9_3p",
                "Z0_5p","Z0_3p","Z1_5p","Z1_3p","Z2_5p","Z2_3p","Z3_5p","Z3_3p","Z4_5p","Z4_3p","Z5_5p","Z5_3p",
                "Z6_5p","Z6_3p","Z7_5p","Z7_3p","Z8_5p","Z8_3p","Z9_5p","Z9_3p")
primer_sets = vector("list",41)
primer_sets[[1]] = 1:2                # 16S_5p.  Note that primer 61 is the same as primer 2
primer_sets[[2]] = 3:4                # 16S_3p
primer_sets[[3]] = 5:10               # 23S_5p
primer_sets[[4]] = 11:14              # 23S_3p
primer_sets[[5]] = 61                 # V1-5p;   Identical to 67. 61 labelled as 27F, 67 as 8F
primer_sets[[6]] = 66                 # V1-3p
primer_sets[[7]] = 63:65              # V2-5p
primer_sets[[8]] = 62                 # V2-3p
primer_sets[[9]] = c(18,22,23)        # V3_5p
primer_sets[[10]] = c(19,24)          # V3_3p
primer_sets[[11]] = 27:36             # V4_5p;  Note that 25:26 targets a supposedly slightly different location - but primer sequences don't match.
primer_sets[[12]] = 37:48             # V4_3p
primer_sets[[13]] = 54                # V5_5p
primer_sets[[14]] = 55:60             # V5_3p
primer_sets[[15]] = 20                # V6_5p
primer_sets[[16]] = 21                # V6_3p
primer_sets[[17]] = 49                # V7_5p
primer_sets[[18]] = 50:53             # V7_3p
primer_sets[[19]] = 68:69             # V9_5p
primer_sets[[20]] = 70                # V9_3p
primer_sets[[21]] = 5:10              # Z0_5p;  Identical to 23S_5'
primer_sets[[22]] = 71:79             # Z0_3p
primer_sets[[23]] = 71:79             # Z1_5p
primer_sets[[24]] = 80:81             # Z1_3p
primer_sets[[25]] = 80:81             # Z2_5p
primer_sets[[26]] = 82:85             # Z2_3p
primer_sets[[27]] = 82:85             # Z3_5p
primer_sets[[28]] = 86:101            # Z3_3p
primer_sets[[29]] = 86:101            # Z4_5p
primer_sets[[30]] = 102:109           # Z4_3p
primer_sets[[31]] = 102:109           # Z5_5p
primer_sets[[32]] = 110:117           # Z5_3p
primer_sets[[33]] = 110:117           # Z6_5p
primer_sets[[34]] = 118:125           # Z6_3p
primer_sets[[35]] = 118:125           # Z7_5p
primer_sets[[36]] = 126               # Z7_3p
primer_sets[[37]] = 126               # Z8_5p
primer_sets[[38]] = 127:134           # Z8_3p
primer_sets[[39]] = 127:134           # Z9_5p
primer_sets[[40]] = 11:14             # Z9_3p    same as 23S_3p
primer_sets[[41]] = 15:17             # 2241R    Unclear 

primers = vector("list",134)   # Note that this length of the list primers must be exact
#primers[[1:2]] = list(type="16S_5p",id="27F_A",sequence="AGAGTTTGATCMTGGCTCAG") 
primers[[1]] = list(type="16S_5p",id="27F_A",sequence="AGAGTTTGATCATGGCTCAG")
primers[[2]] = list(type="16S_5p",id="27F_C",sequence="AGAGTTTGATCCTGGCTCAG") 
#primers[[3:4]] = list(type="16S_3p",id="1492R_C",sequence="TACGGYTACCTTGTTACGACTT")
primers[[3]] = list(type="16S_3p",id="1492R_C",sequence="TACGGCTACCTTGTTACGACTT")   # "GGTTACCTTGTTACGACTT"
primers[[4]] = list(type="16S_3p",id="1492R_T",sequence="TACGGTTACCTTGTTACGACTT")
#primers[[5:10]] = list(type="23S_5p",id="129F_YV",sequence="CYGAATGGGGVAACC")
primers[[5]] = list(type="23S_5p",id="129F_CA",sequence="CCGAATGGGGAAACC")
primers[[6]] = list(type="23S_5p",id="129F_CC",sequence="CCGAATGGGGCAACC")
primers[[7]] = list(type="23S_5p",id="129F_CG",sequence="CCGAATGGGGGAACC")
primers[[8]] = list(type="23S_5p",id="129F_TA",sequence="CTGAATGGGGAAACC")
primers[[9]] = list(type="23S_5p",id="129F_TC",sequence="CTGAATGGGGCAACC")
primers[[10]] = list(type="23S_5p",id="129F_TG",sequence="CTGAATGGGGGAACC")
#primers[[11:14]] = list(type="23S_3p",id="U2482R_AA",sequence="CCRAMCTGTCTCACGACG")
primers[[11]] = list(type="23S_3p",id="U2482R_AA",sequence="CCAAACTGTCTCACGACG")
primers[[12]] = list(type="23S_3p",id="U2482R_AC",sequence="CCAACCTGTCTCACGACG")
primers[[13]] = list(type="23S_3p",id="U2482R_GA",sequence="CCGAACTGTCTCACGACG")
primers[[14]] = list(type="23S_3p",id="U2482R_GC",sequence="CCGACCTGTCTCACGACG")
#primers[[15:17]] = list(type="23S_3p",id="2241R",sequence="ACCGCCCCAGTHAAACT")
primers[[15]] = list(type="23S_xx",id="2241R_A",sequence="ACCGCCCCAGTAAAACT")
primers[[16]] = list(type="23S_xx",id="2241R_C",sequence="ACCGCCCCAGTCAAACT")
primers[[17]] = list(type="23S_xx",id="2241R_T",sequence="ACCGCCCCAGTTAAACT")
#primers[18:21] From Chakravorty S et al., J Microbiol Methods. 2007 May; 69(2): 330–339. 
primers[[18]] = list(type="V3_5p",id="V3F",sequence="CCAGACTCCTACGGGAGGCAG")
primers[[19]] = list(type="V3_3p",id="V3R",sequence="CGTATTACCGCGGCTGCTG")
primers[[20]] = list(type="V6_5p",id="V6F",sequence="TCGATGCAACGCGAAGAA")
primers[[21]] = list(type="V6_3p",id="V6R",sequence="ACATTTCACAACACGAGCTGACGA")
# primers[22:40+] From Zhang J etal Science of the Total Environment 618 (2018) 1254–1267
#  probeBase:  338R  "TGCTGCCTCCCGTAGGAGT" for V3_5p
#  probeBase:  P338F  "ACTCCTACGGGAGGCAGCAG" for V3_5p
#  probeBase:  533R  "TTACCGCGGCTGCTGGCAC" for V3_3p
#  probeBase:  520F  "AYTGGGYDTAAAGNG" for V4_5p
#  probeBase:  802R  "TACNVGGGTATCTAATCC" for V4_3p
# Note that "N" has been replaced only by "A" in 520F primers
primers[[22]] = list(type="V3_5p",id="P338F",sequence="ACTCCTACGGGAGGCAGCAG")
primers[[23]] = list(type="V3_5p",id="338R",sequence="TGCTGCCTCCCGTAGGAGT")
primers[[24]] = list(type="V3_3p",id="533R",sequence="TTACCGCGGCTGCTGGCAC")
# The V4 primers are 515F (5′-GTGCCAGCMGCCGCGGTAA-3′) and 806R (5′-GGACTACHVGGGTWTCTAAT-3′)
#  from Caporaso, J. G., Lauber, C. L., Walters, W. A., Berg-Lyons, D., Lozupone, C. A., Turnbaugh, P. J., et al. (2011). 
#       Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. Proc. Natl. Acad. Sci. U.S.A. 108, 4516–4522. doi: 10.1073/pnas.1000080107
# GTGYCAGCMGCCGCGGTAA   - 515F V4-5p  Parada   (cf proBase GTGCCAGCMGCCGCGGTAA) where M=A or C) 
#  probeBase:  Don't use this one 520F  "AYTGGGYDTAAAGNG" for V4_5p
# Note that "N" has been replaced only by "A" in 520F primers - this is sub-optimal
# 5′-CAGCMGCCGCGGTAATAC-3′   - 519F  : Burggraf S, Huber H, Stetter KO. Reclassification of the crenarchaeal orders and families in accordance with 16S rRNA sequence data. 
#                                      Int J Syst Bacteriol. 1997;47:657–660. doi: 10.1099/00207713-47-3-657
primers[[25]] = list(type="V4_5p",id="515F_CCA",sequence="GTGCCAGCAGCCGCGGTAA")
primers[[26]] = list(type="V4_5p",id="515F_CTA",sequence="GTGCCAGCCGCCGCGGTAA")
primers[[27]] = list(type="V4_5p",id="520F_TCA",sequence="ATTGGGCATAAAGAG")
primers[[28]] = list(type="V4_5p",id="520F_TTA",sequence="ATTGGGTATAAAGAG")
primers[[29]] = list(type="V4_5p",id="520F_CCG",sequence="ACTGGGCGTAAAGAG")
primers[[30]] = list(type="V4_5p",id="520F_CTG",sequence="ACTGGGTGTAAAGAG")
primers[[31]] = list(type="V4_5p",id="520F_TCG",sequence="ATTGGGCGTAAAGAG")
primers[[32]] = list(type="V4_5p",id="520F_TTG",sequence="ATTGGGTGTAAAGAG")
primers[[33]] = list(type="V4_5p",id="520F_CCT",sequence="ACTGGGCTTAAAGAG")
primers[[34]] = list(type="V4_5p",id="520F_CTT",sequence="ACTGGGTTTAAAGAG")
primers[[35]] = list(type="V4_5p",id="520F_TCT",sequence="ATTGGGCTTAAAGAG")
primers[[36]] = list(type="V4_5p",id="520F_TTT",sequence="ATTGGGTTTAAAGAG")
#  probeBase:  802R  "TACNVGGGTATCTAATCC" for V4_3p
# Aprill etal  806R  "GGACTACNVGGGTWTCTAAT" for V4_3p
primers[[37]] = list(type="V4_3p",id="802R_AA",sequence="TACAAGGGTATCTAATCC")
primers[[38]] = list(type="V4_3p",id="802R_AC",sequence="TACACGGGTATCTAATCC")
primers[[39]] = list(type="V4_3p",id="802R_AG",sequence="TACAGGGGTATCTAATCC")
primers[[40]] = list(type="V4_3p",id="802R_CA",sequence="TACCAGGGTATCTAATCC")
primers[[41]] = list(type="V4_3p",id="802R_CC",sequence="TACCCGGGTATCTAATCC")
primers[[42]] = list(type="V4_3p",id="802R_CG",sequence="TACCGGGGTATCTAATCC")
primers[[43]] = list(type="V4_3p",id="802R_GA",sequence="TACGAGGGTATCTAATCC")
primers[[44]] = list(type="V4_3p",id="802R_GC",sequence="TACGCGGGTATCTAATCC")
primers[[45]] = list(type="V4_3p",id="802R_GG",sequence="TACGGGGGTATCTAATCC")
primers[[46]] = list(type="V4_3p",id="802R_TA",sequence="TACTAGGGTATCTAATCC")
primers[[47]] = list(type="V4_3p",id="802R_TC",sequence="TACTCGGGTATCTAATCC")
primers[[48]] = list(type="V4_3p",id="802R_TG",sequence="TACTGGGGTATCTAATCC")
# microbiome: V7_5p with 16S.1100.F16; V7_3p with 1237F
primers[[49]] = list(type="V7_5p",id="16S.1100.F16",sequence="CAACGAGCGCAACCCT")
#primers[[50:53]] = list(type="V7_3p",id="1237F",sequence="GGGCTACACACGYGCWAC"), W=A or T; Y = C or T
primers[[50]] = list(type="V7_3p",id="1237F_CA",sequence="GGGCTACACACGCGCAAC")
primers[[51]] = list(type="V7_3p",id="1237F_CT",sequence="GGGCTACACACGCGCAAC")
primers[[52]] = list(type="V7_3p",id="1237F_TA",sequence="GGGCTACACACGTGCTAC")
primers[[53]] = list(type="V7_3p",id="1237F_TT",sequence="GGGCTACACACGTGCTAC")
# probeBase: V5_5p with 784F;  V5_3p with 908F, 908R
primers[[54]] = list(type="V5_5p",id="784F",sequence="AGGATTAGATACCCT")
#primers[[55:58]] = list(type="V5_3p",id="909F",sequence="ACTCAAAKGAATWGACGG")
primers[[55]] = list(type="V5_3p",id="909F_GA",sequence="ACTCAAAGGAATAGACGG")
primers[[56]] = list(type="V5_3p",id="909F_GT",sequence="ACTCAAAGGAATTGACGG")
primers[[57]] = list(type="V5_3p",id="909F_TA",sequence="ACTCAAATGAATAGACGG")
primers[[58]] = list(type="V5_3p",id="909F_TT",sequence="ACTCAAATGAATTGACGG")
#primers[[59:60]] = list(type="V5_3p",id="908R",sequence="CGTCAATTCMTTTGAGTT"), M=A or C
primers[[59]] = list(type="V5_3p",id="908R_A",sequence="CGTCAATTCATTTGAGTT")
primers[[60]] = list(type="V5_3p",id="908R_C",sequence="CGTCAATTCCTTTGAGTT")
primers[[61]] = list(type="16S-5p",id="27F_C",sequence="AGAGTTTGATCCTGGCTCAG")    # From Graspeuntner S etal "Selection of validated hypervariable 
primers[[62]] = list(type="V2_3p",id="338R",sequence="TGCTGCCTCCCGTAGGAGT")      #  regions is crucial in 16S-based microbiota studies of the female 
                                                                                 #   genital tract" Scientific Reports 26 June 2018
# primers[[63:65]] = list(type="V2_5p",id="97F",sequence="GGCGVACGGGTGAGTAA")    # This is a tentative entry, based on Yang B Wang Qian P-Y Bioinformatics 2016        
primers[[63]] = list(type="V2_5p",id="97F",sequence="GGCGAACGGGTGAGTAA")         # This is a tentative entry, based on Yang B Wang Qian P-Y Bioinformatics 2016        
primers[[64]] = list(type="V2_5p",id="97F",sequence="GGCGCACGGGTGAGTAA")         # This is a tentative entry, based on Yang B Wang Qian P-Y Bioinformatics 2016        
primers[[65]] = list(type="V2_5p",id="97F",sequence="GGCGGACGGGTGAGTAA")         # This is a tentative entry, based on Yang B Wang Qian P-Y Bioinformatics 2016        
#      id="64F",5'- BGY CTW ANR CAT GCA AGT YG")                                 # For V1-5p primer: Probase - multiple ambiguity codes to be resolved
#      id="P63F",sequence="CAGGCCTAACACATGCAAGTC")                               # For V1-5p primer:  Probase - matches one version of 64F
#      id="pBR-V1.ASR",sequence="TTA CTC ACC CGT CCG CCA CT")                    # For V1-3p reverse primer  " Probase
primers[[66]] = list(type="V1_3p",id="pBR-V1.ASR",sequence="TTACTCACCCGTCCGCCACT")
primers[[67]] = list(type="V1_5p",id="P63F",sequence="CAGGCCTAACACATGCAAGTC")      
primers[[68]] = list(type="V9_5p",id="1381R",sequence="GACGGGCGGTGTGTACA")
primers[[69]] = list(type="V9_5p",id="1381R",sequence="GACGGGCGGTGTGTGCA")
primers[[70]] = list(type="V9_3p",id="1492R",sequence="GGTTACCTTGTTACGACTT")
# Now add Zhengming Zhang's Zregion primers.
#primers[[71:79]] = list(type="Z1_5p",id="Z1_84",sequence="GAABTGAAACATCTHAGTA")
primers[[71]] = list(type="Z0_3p",id="Z01F_84",sequence="GAACTGAAACATCTAAGTA")
primers[[72]] = list(type="Z0_3p",id="Z02F_84",sequence="GAAGTGAAACATCTAAGTA")
primers[[73]] = list(type="Z0_3p",id="Z03F_84",sequence="GAATTGAAACATCTAAGTA")
primers[[74]] = list(type="Z0_3p",id="Z04F_84",sequence="GAACTGAAACATCTCAGTA")
primers[[75]] = list(type="Z0_3p",id="Z05F_84",sequence="GAAGTGAAACATCTCAGTA")
primers[[76]] = list(type="Z0_3p",id="Z06F_84",sequence="GAATTGAAACATCTCAGTA")
primers[[77]] = list(type="Z0_3p",id="Z07F_84",sequence="GAACTGAAACATCTTAGTA")
primers[[78]] = list(type="Z0_3p",id="Z08F_84",sequence="GAAGTGAAACATCTTAGTA")
primers[[79]] = list(type="Z0_3p",id="Z09F_84",sequence="GAATTGAAACATCTTAGTA")
#primers[[80:81]] = list(type="Z1_3p",id="Z2_136",sequence="AGTAGYGGCGAGCGAA")
primers[[80]] = list(type="Z1_3p",id="Z01F_136",sequence="AGTAGCGGCGAGCGAA")
primers[[81]] = list(type="Z1_3p",id="Z02F_136",sequence="AGTAGTGGCGAGCGAA")
#primers[[82:85]] = list(type="Z3_5p",id="Z3_384",sequence="AGTACYGTGARGGAA")
primers[[82]] = list(type="Z2_3p",id="Z01F_384",sequence="AGTACCGTGAAGGAA")
primers[[83]] = list(type="Z2_3p",id="Z02F_384",sequence="AGTACTGTGAAGGAA")
primers[[84]] = list(type="Z2_3p",id="Z03F_384",sequence="AGTACCGTGAGGGAA")
primers[[85]] = list(type="Z2_3p",id="Z04F_384",sequence="AGTACTGTGATGGAA")
#primers[[86:101]] = list(type="Z4_5p",id="Z4_728",sequence="WRATAGCTSGTWCTC")
primers[[86]] = list(type="Z3_3p",id="Z01F_728",sequence="AAATAGCTCGTACTC")
primers[[87]] = list(type="Z3_3p",id="Z02F_728",sequence="AAATAGCTCGTTCTC")
primers[[88]] = list(type="Z3_3p",id="Z03F_728",sequence="AAATAGCTGGTACTC")
primers[[89]] = list(type="Z3_3p",id="Z04F_728",sequence="AAATAGCTGGTTCTC")
primers[[90]] = list(type="Z3_3p",id="Z05F_728",sequence="AGATAGCTCGTACTC")
primers[[91]] = list(type="Z3_3p",id="Z06F_728",sequence="AGATAGCTCGTTCTC")
primers[[92]] = list(type="Z3_3p",id="Z07F_728",sequence="AGATAGCTGGTACTC")
primers[[93]] = list(type="Z3_3p",id="Z08F_728",sequence="AGATAGCTGGTTCTC")
primers[[94]] = list(type="Z3_3p",id="Z09F_728",sequence="TAATAGCTCGTACTC")
primers[[95]] = list(type="Z3_3p",id="Z10F_728",sequence="TAATAGCTCGTTCTC")
primers[[96]] = list(type="Z3_3p",id="Z11F_728",sequence="TAATAGCTGGTACTC")
primers[[97]] = list(type="Z3_3p",id="Z12F_728",sequence="TAATAGCTGGTTCTC")
primers[[98]] = list(type="Z3_3p",id="Z13F_728",sequence="TGATAGCTCGTACTC")
primers[[99]] = list(type="Z3_3p",id="Z14F_728",sequence="TGATAGCTCGTTCTC")
primers[[100]] = list(type="Z3_3p",id="Z15F_728",sequence="TGATAGCTGGTACTC")
primers[[101]] = list(type="Z3_3p",id="Z16F_728",sequence="TGATAGCTGGTTCTC")
#primers[[102:109]] = list(type="Z5_5p",id="Z5_986",sequence="GTTRGCTYRGAAGCAGC")
primers[[102]] = list(type="Z4_3p",id="Z01F_986",sequence="GTTAGCTCAGAAGCAGC")
primers[[103]] = list(type="Z4_3p",id="Z02F_986",sequence="GTTAGCTCGGAAGCAGC")
primers[[104]] = list(type="Z4_3p",id="Z03F_986",sequence="GTTAGCTTAGAAGCAGC")
primers[[105]] = list(type="Z4_3p",id="Z04F_986",sequence="GTTAGCTTGGAAGCAGC")
primers[[106]] = list(type="Z4_3p",id="Z05F_986",sequence="GTTGGCTCAGAAGCAGC")
primers[[107]] = list(type="Z4_3p",id="Z06F_986",sequence="GTTGGCTCGGAAGCAGC")
primers[[108]] = list(type="Z4_3p",id="Z07F_986",sequence="GTTGGCTTAGAAGCAGC")
primers[[109]] = list(type="Z4_3p",id="Z08F_986",sequence="GTTGGCTTGGAAGCAGC")
#primers[[110:117]] = list(type="Z6_5p",id="Z6_1588",sequence="AGGAAYTMKGCAA")
primers[[110]] = list(type="Z5_3p",id="Z01F_1588",sequence="AGGAACTAGGCAA")
primers[[111]] = list(type="Z5_3p",id="Z02F_1588",sequence="AGGAACTATGCAA")
primers[[112]] = list(type="Z5_3p",id="Z03F_1588",sequence="AGGAACTCGGCAA")
primers[[113]] = list(type="Z5_3p",id="Z04F_1588",sequence="AGGAACTCTGCAA")
primers[[114]] = list(type="Z5_3p",id="Z05F_1588",sequence="AGGAATTAGGCAA")
primers[[115]] = list(type="Z5_3p",id="Z06F_1588",sequence="AGGAATTATGCAA")
primers[[116]] = list(type="Z5_3p",id="Z07F_1588",sequence="AGGAATTCGGCAA")
primers[[117]] = list(type="Z5_3p",id="Z08F_1588",sequence="AGGAATTCTGCAA")
#primers[[118:125]] = list(type="Z7_5p",id="Z7_1723",sequence="TGAYRCCTGCCCRGTGC")
primers[[118]] = list(type="Z6_3p",id="Z01F_1723",sequence="TGACACCTGCCCAGTGC")
primers[[119]] = list(type="Z6_3p",id="Z02F_1723",sequence="TGACACCTGCCCGGTGC")
primers[[120]] = list(type="Z6_3p",id="Z03F_1723",sequence="TGACGCCTGCCCAGTGC")
primers[[121]] = list(type="Z6_3p",id="Z04F_1723",sequence="TGACGCCTGCCCGGTGC")
primers[[122]] = list(type="Z6_3p",id="Z05F_1723",sequence="TGATACCTGCCCAGTGC")
primers[[123]] = list(type="Z6_3p",id="Z06F_1723",sequence="TGATACCTGCCCGGTGC")
primers[[124]] = list(type="Z6_3p",id="Z07F_1723",sequence="TGATGCCTGCCCAGTGC")
primers[[125]] = list(type="Z6_3p",id="Z08F_1723",sequence="TGATGCCTGCCCGGTGC")
primers[[126]] = list(type="Z7_3p",id="Z01F_1821",sequence="TCCTAAGGTAGCGAAATTCCTTG")
#primers[[127:134]] = list(type="Z9_5p",id="Z9_2144",sequence="ACTGGGGYGGTYKCCTCC")
primers[[127]] = list(type="Z8_3p",id="Z01F_2144",sequence="ACTGGGGCGGTCGCCTCC")
primers[[128]] = list(type="Z8_3p",id="Z02F_2144",sequence="ACTGGGGCGGTCTCCTCC")
primers[[129]] = list(type="Z8_3p",id="Z03F_2144",sequence="ACTGGGGCGGTTGCCTCC")
primers[[130]] = list(type="Z8_3p",id="Z04F_2144",sequence="ACTGGGGCGGTTTCCTCC")
primers[[131]] = list(type="Z8_3p",id="Z05F_2144",sequence="ACTGGGGTGGTCGCCTCC")
primers[[132]] = list(type="Z8_3p",id="Z06F_2144",sequence="ACTGGGGTGGTCTCCTCC")
primers[[133]] = list(type="Z8_3p",id="Z07F_2144",sequence="ACTGGGGTGGTTGCCTCC")
primers[[134]] = list(type="Z8_3p",id="Z08F_2144",sequence="ACTGGGGTGGTTTCCTCC")

