Dataset Title: Diet DNA metabarcoding data from spider individuals (Heteropoda venatoria) either surface sterilized with bleach or not from Palmyra Atoll 2015-2017

Abstract: These are data and code from a study examining the potential for surface contamination to influence diet DNA metabarcoding datasets when DNA is sequenced from full body parts (in this case, the opisthosomas of spider individuals). These datasets include the raw sequencing data, all downstream datasets, and taxonomic assignments collected from database searches on BOLD and GenBank (accessed 2019). The code includes code to reproduce all bioinformatics (merge, filter, match to taxonomies, rarefy, sort) as well as all statistics and figures generated from analyses. Raw data are from DNA extractions of predator gut regions (opisthosomas) and amplification of the CO1 gene using PCR. The predator species is Heteropoda venatoria collected individually with sterilized implements and either the diet sequences from their natural diets were extracted or their diets following feeding spiders in a feeding trial. 

Creators: Miller-ter Kuile, Ana, Austen Apigo, Hillary Young

Other personnel names and roles:
Field technicians: Jasmine Childress, Katherine Plummer, Magalay Espinoza, John McLaughlin, Michelle Lee, Carina Motta
Lab technicians: Tessa Chou, Emily Lutz

Contact: Ana Miller-ter Kuile, ana.miller.ter.kuile@gmail.com

Keywords: consumptive interactions, invertebrates, contamination, food web, predator prey interactions, diet analysis

Funding of this work:
Hillary S. Young, National Science Foundation, DEB-1457371
Rodolfo Dirzo, National Geographic Society
Hillary S. Young, UC Santa Barbara Faculty Senate Grant

Timeframe:
Being date: July 2015
End date: Dec 2021

Geographic Location:
Verbal description: Palmyra Atoll National Wildlife Refuge, Northern Line Islands
North bounding coordinate: 5.883333
South bounding coordinate: -162.08333

Taxonomic species or groups:
Phylum: Arthropoda, Class: Arachnida, Order: Araneae, Family: Sparassidae, Genus: Heteropoda, Species: venatoria

Methods:
Field site and collections
We conducted fieldwork on Palmyra Atoll National Wildlife Refuge, Northern Line Islands (5º53’ N, 162º05’W). Palmyra Atoll has a well-characterized species list and is relatively species poor, allowing for relatively complete characterization of consumer and diet items (Handler, Gruner, Haines, Lange, & Kaneshiro, 2007). We targeted a generalist, active hunting spider species (Heteropoda venatoria) because a) it occurs in high abundance on the atoll and is easy to collect, b) it is a generalist species that feeds on a wide suite of organisms (including spiders, other invertebrates, and two geckos in the genus Lepidodactylus), c) it is the only species in its family on the atoll, meaning consumer DNA can be differentiated from potential diet DNA. All individuals were stored individually in sterilized containers (Greenstone et al., 2011). 

Natural environment consumer collection
In 2015, we collected consumers (n = 47) from natural environments, which had fed on available diet items and come into contact with environmental surfaces, to test if DNA metabarcoding detects diet DNA effectively. We froze all individuals at -80C immediately following collection until surface sterilization and DNA extraction in 2019.

Feeding trial consumer set-up and feeding
In 2017, we conducted laboratory trials (n = 26) to test if DNA metabarcoding detects DNA from diet items offered in a contained environment. We created feeding environments from one-liter yogurt containers with holes for air transfer, and placed one H. venatoria in each container. After 12 hours, we placed one large grasshopper (Oxya japonica, a likely diet item (Handler et al., 2007)) in each container and left all containers for 24 hours. We then froze (-20C) each H. venatoria that had killed the grasshopper (n = 25, consumption was not easily detectable and thus not considered in analyses). We cleaned all containers between trials with 10% bleach solution.

To test surface sterilization’s efficacy at removing possible contaminants, we used a surface sterilization treatment (Schulz, Wanke, Draeger, & Aust, 1993; Burgdorf et al., 2014) on ~half the consumers for each set: those collected from the natural environment, and those subjected to controlled feeding trials.  We submerged and stirred each (whole) consumer in 10% commercial bleach by volume for two minutes and washed each in deionized water for two minutes. Similar bleach submersion leads to undetectable DNA degradation in similar soft-exoskeleton consumers (Greenstone et al., 2012; Linville & Wells, 2002). Natural environment consumers (2015) had been frozen at -80C since collection; we surface sterilized these consumers in a sterilized laminar flow hood in 2019 just before DNA extraction (n = 22 surface sterilized, n = 25 not surface sterilized; Table 1). We surface sterilized feeding trial consumers (2017) in the lab on the atoll in 2017 following freezing at -20C, then stored each in individual vials of 95% ethanol in a -20ºC freezer until DNA extraction (no -80C freezer was available at the field station that year) (n = 10 surface sterilized; n = 14 not surface sterilized). Prior to DNA extraction, all samples dried for 1-3 hours in a sterilized laminar flow hood and the opisthosoma (hind gut region) was removed using a sterilized scalpel.  Between all steps, tools were sterilized with either ethanol and flame (scalpels and forceps) or 10% bleach (surfaces) between handling each individual. 

DNA extraction and removal of consumer DNA with Ampure XP beads
We extracted DNA from each consumer following a modified CTAB extraction protocol (Fulton, Chunwongse, & Tanksley, 1995). We quantified DNA using a Qubit (Invitrogen) fluorometer with the high sensitivity double-stranded DNA quantification kit. We followed Krehenwinkel et al. (2017) to isolate a proportion of lower molecular weight DNA with Ampure XP beads prior to PCR (Appendix E, Figure 1). We diluted each DNA sample to 20ng/microliter (creating a total sample volume of 40microLiter), mixed each sample using Ampure XP beads (0.75x bead-to-DNA ratio), and kept the supernatant. With the supernatant, we precipitating the DNA pellets with isopropanol and 5M potassium acetate, and washed DNA pellets with ethanol (Appendix F). We quantified this cleaned DNA again using a Qubit fluorometer and diluted all samples to 10ng/microLiter prior to PCR steps. All DNA pellets were stored in and diluted with TE buffer.

PCR amplification, library preparation, and sequencing
We amplified the CO1 gene with general metazoan primers (Yu et al., 2012; Leray et al., 2013; Krehenwinkel et al., 2017, Table 2). We performed all PCR preparation steps in a UV-sterilized biosafety cabinet. We used PCR reaction volumes of 25microLiter (9microLiter nuclease free water, 12.5microLiter GoTaq Green Master Mix (Promega Corp.), 1.25 microLiter of each of the primers (at 10mM), and 1 microLiter of DNA template (at 10ng/microLiter)). We ran each sample in duplicate along with duplicated negative samples each PCR run. PCR reactions are as follows: initial denaturation step at 95C for 3 minutes, then 35 cycles of: 1) 95C for 30 seconds, 2) 46C for 30 seconds, and 3) 72C for one minute, followed by a final 5 minutes at 72C. We cleaned PCR products with Ampure XP beads at a 0.8x bead-to-DNA ratio and resuspended from beads using a 10mM TRIS buffer. 

We attached Illumina index primers with an additional PCR step following standard protocols (Nextera XT Index Kit v2, Illumina, 2019). We combined duplicate samples for which both duplicates successfully amplified and diluted to a concentration of 5nM. We multiplexed all samples with one negative control and two fungal clone positive controls (GenBank accession numbers: MG840195 and MG840196;  Apigo & Oono, 2018; Clark et al., 2016; Toju et al., 2012). We submitted multiplexed samples for sequencing at the University of California, Santa Barbara Biological Nanostructures Laboratory Genetics Core. Samples were run on an Illumina MiSeq platform (v2 chemistry, 500 cycles, paired-end reads) with a 15% spike-in of PhiX. Following sequencing, samples were demultiplexed using Illumina’s bcl2fastq conversion software (v2.20) at the Core facility. Our full protocol from DNA extraction through submission for Illumina sequencing can be found in Appendix F. 

Sequence merging, filtering, and clustering with UNOISE3
We merged, filtered (max ee  = 1.0), and denoised (clustered) our sequences around amplicon sequence variants (ASVs) using the UNOISE3 algorithm (unoise3 command in the open-source USEARCH 32-bit version 11.0.667; Edgar, 2016, Appendix E, Figure 3). Prior to denoising with UNOISE3, we used cutadapt (version 1.18, Martin, 2011) to remove primers from each sequence. We also repeated analyses with the DADA2 algorithm run through R (dada2 package version 1.1.14.0; Callahan et al., 2016) and with a data cleaning step run through BBSplit (Bushnell, 2019) to remove consumer DNA prior to ASV assignment  (because ASV assignment is abundance-sensitive). We considered analyses from the UNOISE3 algorithm only because UNOISE3 assigned more sequence reads to positive controls than DADA2 (on average, 3x as many reads per positive control) and the cleaning step paired with either DADA2 or UNOISE3 did not increase potential diet DNA detection (summary and comparisons in Appendices A and B). 

We created a list of unique ASVs and a matrix of ASV abundances across samples. We matched ASVs to taxonomies in the GenBank and BOLD databases. For GenBank, we used BLAST (version 2.7.1) with the blastn command for taxonomic assignment of each ASV using the computing cluster at UC Santa Barbara, comparing against the GenBank nucleotide database with an evalue of 0.01 (downloaded on November 20, 2019). We visualized and exported taxonomic alignment using MEGAN Community Edition (version 6.18.0, Huson et al., 2016), using default settings (LCA=naïve, MinScore = 50.0, MaxExpected  = 0.01, TopPercent = 10.0, MinSupportPercent = 0.05) and selecting the subtree with all possible diet items for this species (Kingdom: Animalia, Clade: Bilateria). For taxonomies which were not assigned below the order level (n =24), we submitted each ASV individually to the BLAST Basic Local Alignment Search Tool and assigned them a family based on the best sequence match in the database, given that the top ten database matches were from the same family. For BOLD taxonomic assignment, we used the BOLD IDEngine of the CO1 gene with Species Level Barcode Records (accessed February 5-16, 2020; 3,825,490 Sequences, 216,704 Species, and 95,537 Interim Species in database) to match each ASV list to taxonomies. We combined taxonomic assignments from both programs and discarded taxonomic assignments that were mismatched at the family level or higher (Elbrecht, Peinert, & Leese, 2017). 

Detection of potential diet items
For consumers from both the natural environment and feeding trials, we asked whether surface sterilization altered detection of potential diet items for each consumer. For natural environment consumers, we examined all potential diet items (which could represent either diet or surface contaminants). For feeding trial consumers, we focused our detection analysis on the offered diet item we provided the consumers in the feeding trial environment (O. japonica, which all consumers were observed to have killed, but not necessarily ingested). We rarefied (McKnight et al., 2019, Appendix E, Figure 4) based on the sample with the lowest sequencing depth which had been sequenced with 95%+ sampling completeness based on iNEXT (version 2.0.20) interpolation and extrapolation methods (Hsieh & Chao, 2017, 16,004 reads for natural environment and 55,205 reads for feeding trial consumers). We rarefied using the rrarefy() function in the vegan (version 2.5.6) package in R and rarefied the field and lab consumers separately. 

We then selected all ASVs that matched potential diet items for the natural environment consumers (Kingdom: Animalia; Clade: Bilateria, excluding consumer DNA) and just the offered diet item for the feeding trial consumers (including species: Oxya japonica, genus: Oxya, family: Acrididae, excluding those which only matched to order). Because the consumer species H. venatoria is the only species in the family Sparassidae on Palmyra Atoll, removing consumer DNA meant excluding all ASVs that received a family-level taxonomic assignment of “Sparassidae”. As all ASVs received family-level taxonomic assignment, we pooled ASVs that matched at the family level into one taxonomic unit using cumulative read abundance (i.e. all ASVs matched to diet family A were pooled into diet family A taxonomic unit), a practice common in diet metabarcoding (Kartzinel et al., 2015) and predator-prey interaction (Brose et al., 2019) studies.

Statistical analyses
For potential diet detection and rarefied abundance in both sets of consumers (natural environment and feeding trial) we used generalized linear models to assess the effect of surface sterilization treatment. For prey detection, we used all potential (natural environment) or offered (feeding trial) diet item detection (presence-absence per sample) as the response variable in the full model with surface sterilization as a fixed effect and a binomial distribution. For rarefied diet abundance, we only assessed consumers for which we had detected diet and not those with no diet detection (n = 33 of 37 for natural environment; n = 14 of 19 for feeding trials). For this model, we treated the number of all potential (natural environment) or offered (feeding trials) diet DNA reads per sample as the response variable, surface sterilization treatment as a fixed effect, total read abundance of the sample (constant across all) as an offset term, and a Poisson or negative binomial distribution (to correct for overdispersion when needed). We assessed differences in per sample potential diet richness among sterilization treatments for the natural environment consumers using generalized linear models with the number of potential diet items per sample as the response variable (both family-level taxonomic units or ASVs), surface sterilization treatment as the fixed effect and a Poisson or negative binomial distribution (to correct for overdispersion when needed). We assessed differences in potential diet item composition with family-level taxonomic units between surface sterilized and unsterilized consumers using a presence-absence PERMANOVA model fit with a binomial mixed effects model with surface sterilization treatment as a fixed effect, a random intercept term for potential diet item, and a random slope term for surface sterilization treatment. We also assessed ASV composition as a representation of potential prey composition using a canonical correspondence analysis (CCA) with surface sterilization as a predictor variable. We performed these analyses along with multiple other supplementary analyses and approaches, which can be found in the Supplementary Information (Appendix D and E).  

For all generalized linear models and mixed models, we performed model selection by comparing the full model (including the fixed effect of surface sterilization treatment) to a null model without this effect. All models were called in the glmmTMB package (version 1.0.0, Brooks et al., 2017) in R (version 3.6.1) We chose the best fitting model based on size corrected AIC values (MuMIn package version 1.43.15). For responses for which the best model included the surface sterilization treatment term, we examined the model summary to determine the standardized coefficients (β) and p-value of the significance between marginal means of the levels of the surface sterilization fixed effect. We assessed model fit using diagnostics in the DHARMa package (version 0.2.7), including tests for heteroskedasticity, and for count models (Poisson or negative binomial), zero inflation and overdispersion (Bolker et al., 2009; Zuur, Ieno, Walker, Saveliev, Anatoly, & Smith, 2009). We performed the CCA analysis using the vegan package in R, comparing a model with surface sterilization as a fixed effect to a null model using an ANOVA. All raw data, data cleaning, and data analyses are available online (Miller-ter Kuile, 2020b, 2020a), and model outputs for primary and supplemental models can be found in Appendices C and D.

### 

Data Tables
*Note: We did not include summaries of data tables in the “outputs” folder of the data folder. These are all intermediate steps in the data preparation process and the steps to generate them are found in the scripts in the ‘2_data_prep’ folder.

# 
Table name: dada2_c_pred_asvs.fasta, dada2_c_prey_asvs.fasta, dada2_c_um_asvs.fasta, dada2_uc_asvs.fasta

Table description: These are all lists of sequences generated from the DADA2 pipeline. Any file with a "_c_" indicates these data were cleaned with BBSplit prior to DADA2, the one with "_uc_" was not run through BBSplit. (c_pred = predator DNA split during BBSplit, c_prey = prey DNA split during BBSplit, c_um = DNA not matched to either predator or prey during BBSplit)

Column name | Description | Unit or code explanation or date format | Missing value code

ASV | A CO1 sequence merged and filtered via DADA2 | NA | NA

# 
Table name: unoise_c_pred_zotus.fasta, unoise_c_prey_zotus.fasta, unoise_c_um_zotus.fasta, unoise_uc_zotus.fasta

Table description: These are all lists of sequences generated from the UNOISE3 pipeline. Any file with a "_c_" indicates these data were cleaned with BBSplit prior to UNOISE3, the one with "_uc_" was not run through BBSplit. (c_pred = predator DNA split during BBSplit, c_prey = prey DNA split during BBSplit, c_um = DNA not matched to either predator or prey during BBSplit)

Column name | Description | Unit or code explanation or date format | Missing value code

ZOTU | A CO1 sequence merged and filtered via UNOISE3 | NA | NA

# 
Table name: dada2_c_pred_asv_tab.tsv, dada2_c_prey_asv_tab.tsv, dada2_c_um_asv_tab.tsv, dada2_uc_asv_tab.tsv

Table description: A sample by ASV matrix of sequences merged through the DADA2 pipeline. Any file with a "_c_" indicates these data were cleaned with BBSplit prior to DADA2, the one with "_uc_" was not run through BBSplit. (c_pred = predator DNA split during BBSplit, c_prey = prey DNA split during BBSplit, c_um = DNA not matched to either predator or prey during BBSplit)

Column name | Description | Unit or code explanation or date format | Missing value code

ASV | An ASV created through DADA2 | NA | NA

columns with some version of a name "*HEV-[number]_S[number] | a sample that was sequenced through DADA2 | A count of the number of sequences matched to that ASV in that sample | NA

# 
Table name: unoise_c_pred_zotu_tab.txt, unoise_c_prey_zotu_tab.txt, unoise_c_um_zotu_tab.txt, unoise_uc_zotu_tab.txt

Table description: A sample by ASV matrix of sequences merged through the UNOISE3 pipeline. Any file with a "_c_" indicates these data were cleaned with BBSplit prior to UNOISE3, the one with "_uc_" was not run through BBSplit. (c_pred = predator DNA split during BBSplit, c_prey = prey DNA split during BBSplit, c_um = DNA not matched to either predator or prey during BBSplit)

Column name | Description | Unit or code explanation or date format | Missing value code

#OTU ID | An ZOTU created through UNOISE3 | NA | NA

columns with some version of a name "*HEV[number] | a sample that was sequenced through UNOISE3 | A count of the number of sequences matched to that ZOTU in that sample | NA


# 
Table name:Sample_Metadata.csv

Table description: Information on where and when each sample was collected as well as the spider's size, whether it was surface sterilized prior to DNA extraction, and whether it was part of the lab feeding trial or field-collected

Column name | Description | Unit or code explanation or date format | Missing value code

Method | A method of collection for each sample | Hand = collected via individual hand collection methods with sterilized implements | NA

Island | The islet site name where the sample was collected | NA | NA

Year | The year in which the sample was collected | YYYY | NA

ID | The species ID for each sample that was collected | NA | NA

Date Collected | The date the sample was collected in the field | MM/DD/YY | NA

sample | The sample ID that corresponds to the sample's name in DNA sequencing pipelines | NA | NA

Size | The size of each spider individual, measured as a length in millimeters from the front of the head region to the end of the opisthosoma | millimeters | NA

Extract.Date | The date that the sample DNA was extracted | MM/DD/YY |NA | NA

Sterilized | Whether the sample was surface sterilized with bleach washes prior to DNA extraction | SS = surface sterilized, NS = not surface sterilized | NA

Source | Whether the sample was part of the lab feeding trial or the field collection | LAB = lab feeding trial, FIELD = field-collected | NA

# 
Table name: all tables of: IDEngine_Results_Summary ([NUMBER]).xls

Table description: Sequence taxonomic matches from BOLD sequence database

Column name | Description | Unit or code explanation or date format | Missing value code

Query ID | The ASV or ZOTU matched to sequences on BOLD | NA | NA

Best ID | The best match from the BOLD database | NA | No match = no matches in database

Search DB | The database on BOLD that was searched | NA | NA

Top % | The top percent sequence match to the sequence | NA | blank cells

Low % | The lower bound of the percent sequence match for this sequence | NA | blank cells

# 
Table name: dada2_c_pred_taxa_bold.csv, dada2_c_prey_taxa_bold.csv, dada2_c_um_taxa_bold.csv, dada2_uc_taxa_bold.csv, unoise3_c_pred_bold_taxa.csv, unoise3_c_prey_bold_taxa.csv, unoise3_c_um_bold_taxa.csv, unoise_uc_bold_taxa.csv

Table description: Tables of the taxonomic identifications at multiple levels of the ASVs or ZOTUs from UNOISE3 or DADA2, manually-entered after the initial BOLD sequence match process

Column name | Description | Unit or code explanation or date format | Missing value code

Query ID | An ASV or ZOTU input into the BOLD database for taxonomic matches | NA | NA

ID | The identity of the sequence from the BOLD database | NA | NA

Level | The taxonomic level the sequence was matched to in the BOLD database | NA | NA

Kingdom | The kingdom of the taxonomic match | NA | NA

Phylum | The phylum of the taxonomic match | NA | NA

Class | The class of the taxonomic match | NA | NA

Order | The order of the taxonomic match | NA | NA

Family | The family of the taxonomic match | NA | NA

Genus | The genus of the taxonomic match | NA | NA

Species | The species of the taxonomic match | NA | NA


# 
Table name: [pipeline name].rma6

Table description: A file that can be opened in MEGAN community edition to visualize and export all or portions of taxonomic trees of those taxonomic assignments

Column name | Description | Unit or code explanation or date format | Missing value code

No columns

# 
Table name: dada2_c_pred_ALL_ncbi_taxa.csv, dada2_c_prey_ALL.csv, dada2_c_prey_ncbi_taxa.csv, dada2_c_um_ALL.csv, dada2_c_um_ncbi_taxa.csv, dada2_uc_all.csv, dada2_uc_ncbi_taxa.csv, unoise_c_pred_ncbi_ALL_taxa.csv, unoise_c_prey_ALL_ncbi_taxa.csv, unoise_c_um_ncbi_ALL_taxa.csv, unoise_uc_ncbi_taxa.csv

Table description: Tables of taxonomic assignments through either pipeline matched to taxonomies through the NCBI GenBank nucleotide database.

Column name | Description | Unit or code explanation or date format | Missing value code

ASV | An ASV/ZOTU submitted to GenBank | NA | NA

ID | A taxonomic identity drawn from the blastN command | NA | NA

Category | a category of type of DNA assigned to this sequence | non-diet = a sequence that is not diet, predator = predator DNA, prey = prey/diet DNA | NA | NA

Level | The taxonomic level at which the taxonomic assignment was returned | NA | NA

Kingdom | The kingdom of the taxonomic assignment | NA |NA

Phylum | The phylum of the taxonomic assignment | NA | NA

Class | THE class of the taxonomic assignment | NA |NA 

Order | The order of the taxonomic assignment | NA | NA

Family | The family of the taxonomic assignment | NA | NA

Genus | The genus of the taxonomic assignment | NA | NA

Species | The species of the taxonomic assignment | NA | NA

*sometimes, some version of these (or b_[other letters]):

BLAST_ID | An identifier from a single ASV BLAST webpage search of the particular ASV | NA | NA

BLAST_Category | The taxonomic category of the single ASV BLAST webpage search | NA | NA

BLAST_Level | The taxonomic level of the single ASV BLAST webpage search | NA |NA

BLAST_Kingdom | The kingdom of the single ASV BLAST webpage search | NA | NA 

BLAST_Phylum | The phylum of the single ASV BLAST webpage search | NA | NA 

BLAST_Class | The class of the single ASV BLAST webpage search | NA | NA 

BLAST_Order | The order of the single ASV BLAST webpage search | NA | NA 

BLAST_Family | The family of the single ASV BLAST webpage search | NA | NA 

BLAST_Genus | The genus of the single ASV BLAST webpage search | NA | NA 

BLAST_Species | The species of the single ASV BLAST webpage search | NA | NA 


# 
Table name: [name].tre

Table description: Taxonomic trees exported from MEGAN community edition which can be imported into R or displayed in MEGAN.

Column name | Description | Unit or code explanation or date format | Missing value code

No columns

###

Scripts/code

#
File name: cutadapt_primer_trimming.txt
Description: Script to remove primers from DNA sequences prior to merging and filtering
Category: Bioinformatics
Scripting language: Bash/command line

#
File name: bbsplit_loop.txt
Description: A script which sorts samples prior to merging and filtering ("cleaning" step)
Category: Bioinformatics
Scripting language: bbsplit in bash/command line

#
File name: bbsplit_ref_creation.R
Description: A script that creates the prey and predator reference categories for bbpsplit and is based on unclean (no bbsplit) taxonomic assignments
Category: Bioinformatics
Scripting language: R

#
File name: dada2_unclean.R
Description: A script that runs the dada2 algorithm on the unclean (no bbsplit) DNA sequences
Category: Bioinformatics
Scripting language: R

#
File name: usearch_script.txt
Description: The script used to run the USEARCH3 algorithm. Both versions of this script are the same
Category: Bioinformatics
Scripting language: the free version of USEARCH (32-bit) in bash/command line 

#
File name: dada2_cleaned.R
Description: The dada2 algorithm run on the cleaned (run through bbsplit) DNA sequences
Category: Bioinformatics
Scripting language: R

#
File name: 2a_taxonomic_assignments.R
Description: A script that attaches the taxonomic assignments from GenBank and BOLD to the per sample ASV counts 
Category: Data preparation
Scripting language: R

#
File name: 2b_divide_rarefy.R
Description: A script that divides out the field-collected from the lab-fed samples and rarefies samples to similar sampling depths within each experiment
Category: Data preparation
Scripting language: R

#
File name: 2c_rarefied_taxonomically_sorted.R
Description: A script that divides the predator, prey, and unmapped DNA from each other for each set of samples
Category: Data preparation
Scripting language: R

#
File name: 3a_asv_table_cleaning.R
Description: This script checks the "error" in sequencing suggested by the number of reads and ASVs assigned to positive and negative control samples. Both matched really well so we did not adjust our reads
Category: QC
Scripting language: R

#
File name: 3b_sampling_depth.R
Description: This code looks at sampling depth across samples for the raw sequence reads to ensure that samples were sequenced sufficiently and to decide if there should be a cutoff to remove low-sequenced samples
Category: Data preparation
Scripting language: R

#
File name: 4a_prey_detection.R
Description: This script analyses prey detection in the lab and field sterilization experiments
Category: Statistical analyses and figures
Scripting language: R

#
File name: 4b_abundance_analyses.R
Description: This script looks at the abundance of prey DNA in the lab and field sterilization experiments
Category: Statistical analyses and figures
Scripting language: R

#
File name: 4c_field_composition.R
Description: This script looks at the species richness and composition of the field prey DNA with surface sterilization treatment
Category: Statistical analyses and figures
Scripting language: R

#
File name: 4d_field_richness.R
Description: This script un-aggregates the above analyses, instead looking at ASVs as opposed to combined taxonomic groups
Category: Statistical analyses and figures
Scripting language: R

#
File name: 4e_supp_lab_composition.R
Description: This script analyzes the prey composition of the lab-fed spiders
Category: Statistical analyses and figures
Scripting language: R

#
File name: 4f_amplification_success.R
Description: Summary stats of amplification success across samples in our study
Category: Statistical analyses and figures
Scripting language: R

#
File name: 4g_sample_size_stuff.R
Description: a power test for sample size 
Category: Statistical analyses and figures 
Scripting language: R

###

References

Apigo, A., & Oono, R. (2018). MG840195 and MG840196.

Bolker, B. M., Brooks, M. E., Clark, C. J., Geange, S. W., Poulsen, J. R., Stevens, M. H. H., & White, J. S. S. (2009). Generalized linear mixed models: a practical guide for ecology and evolution. Trends in Ecology and Evolution, 24(3), 127–135. doi: 10.1016/j.tree.2008.10.008

Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., … Bolker, B. M. (2017). Modeling Zero-Inflated Count Data With glmmTMB. BioRxiv, 132753. doi: 10.1101/132753

Burgdorf, R. J., Laing, M. D., Morris, C. D., & Jamal-Ally, S. F. (2014). A procedure to evaluate the efficiency of surface sterilization methods in culture-independent fungal endophyte studies. Brazilian Journal of Microbiology, 45(3), 977–983. doi: 10.1590/S1517-83822014000300030

Bushnell, B. (2019). BBMap.

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581–583. doi: 10.1038/nmeth.3869

Clark, K., Karsch-Mizrachi, I., Lipman, D. J., Ostell, J., & Sayers, E. W. (2016). GenBank. Nucleic Acids Research, 44(D1), D67–D72. doi: 10.1093/nar/gkv1276
Edgar, R. C. (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. BioRxiv, 081257. doi: 10.1101/081257

Elbrecht, V., Peinert, B., & Leese, F. (2017). Sorting things out: Assessing effects of unequal specimen biomass on DNA metabarcoding. Ecology and Evolution, 7(17), 6918–6926. doi: 10.1002/ece3.3192

Fulton, T. M., Chunwongse, J., & Tanksley, S. D. (1995). Microprep protocol for extraction of DNA from tomato and other herbaceous plants. Plant Molecular Biology Reporter, 13(3), 207–209. doi: 10.1007/BF02670897

Greenstone, M. H., Weber, D. C., Coudron, T. A., Payton, M. E., & Hu, J. S. (2012). Removing external DNA contamination from arthropod predators destined for molecular gut-content analysis. Molecular Ecology Resources, 12(3), 464–469. doi: 10.1111/j.1755-0998.2012.03112.x

Greenstone, M. H., Weber, D. C., Coudron, T. C., & Payton, M. E. (2011). Unnecessary roughness? Testing the hypothesis that predators destined for molecular gut-content analysis must be hand-collected to avoid cross-contamination. Molecular Ecology Resources, 11(2), 286–293. doi: 10.1111/j.1755-0998.2010.02922.x

Handler, A., Gruner, D., Haines, W., Lange, M., & Kaneshiro, K. (2007). Arthropod surveys on Palmyra Atoll, Line Islands, and insights into the decline of the native tree Pisonia grandis (Nyctaginaceae). Pacific Science, 61(4), 485–502. doi: 10.2984/1534-6188(2007)61
Hsieh, T. C., & Chao, A. (2017). Rarefaction and extrapolation: Making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. Systematic Biology, 66(1), 100–111. doi: 10.1093/sysbio/syw073

Huson, D. H., Beier, S., Flade, I., Górska, A., El-Hadidi, M., Mitra, S., … Tappu, R. (2016). MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data. PLoS Computational Biology, 12(6), 1–12. doi: 10.1371/journal.pcbi.1004957

Illumina. (2019). Illumina. Nextera XT DNA Library Prep Reference Guide., (May), Document # 15031942 v05.

Krehenwinkel, H., Kennedy, S., Pekár, S., & Gillespie, R. G. (2017). A cost‐efficient and simple protocol to enrich prey DNA from extractions of predatory arthropods for large‐scale gut content analysis by Illumina sequencing. Methods in Ecology and Evolution, 8(1), 126–134. doi: 10.1111/2041-210X.12647

Leray, M., Yang, J. Y., Meyer, C. P., Mills, S. C., Agudelo, N., Ranwez, V., … Machida, R. J. (2013). A new versatile primer set targeting a short fragment of the mitochondrial COI region for metabarcoding metazoan diversity: Application for characterizing coral reef fish gut contents. Frontiers in Zoology, 10(34), 1–14. doi: 10.1186/1742-9994-10-34

Linville, J. G., & Wells, J. D. (2002). Surface sterilization of a maggot using bleach does not interfere with mitochondrial DNA analysis of crop contents. Journal of Forensic Sciences, 47(5), 15532J. doi: 10.1520/jfs15532j

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBNet Journal, 17(1), 10–12. doi: doi:10.14806/ej.17.1.200

McKnight, D. T., Huerlimann, R., Bower, D. S., Schwarzkopf, L., Alford, R. A., & Zenger, K. R. (2019). Methods for normalizing microbiome data: An ecological perspective. Methods in Ecology and Evolution, 10(3), 389–400. doi: 10.1111/2041-210X.13115

Miller-ter Kuile, A. (2020a). BioProject: PRJNA639981.

Miller-ter Kuile, A. (2020b). DNA_Diet_Methods.

Schulz, B., Wanke, U., Draeger, S., & Aust, H. J. (1993). Endophytes from herbaceous plants and shrubs: effectiveness of surface sterilization methods. Mycological Research, 97(12), 1447–1450. doi: 10.1016/S0953-7562(09)80215-3

Toju, H., Tanabe, A. S., Yamamoto, S., & Sato, H. (2012). High-coverage ITS primers for the DNA-based identification of ascomycetes and basidiomycetes in environmental samples. PLoS ONE, 7(7). doi: 10.1371/journal.pone.0040863

Yu, D. W., Ji, Y., Emerson, B. C., Wang, X., Ye, C., Yang, C., & Ding, Z. (2012). Biodiversity soup: Metabarcoding of arthropods for rapid biodiversity assessment and biomonitoring. Methods in Ecology and Evolution, 3(4), 613–623. doi: 10.1111/j.2041-210X.2012.00198.x

Zuur, A. F., Ieno, E. N., Walker, N. J., Saveliev, Anatoly, A., & Smith, G. M. (2009). Mixed Effects Models and Extensions in Ecology with R. In Mixed Effects Models and Extensions in Ecology with R (Vol. 53). doi: 10.1017/CBO9781107415324.004



