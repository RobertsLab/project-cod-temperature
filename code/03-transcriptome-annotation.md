Oh What a Blast!
================
Steven Roberts
17 April, 2024

- <a href="#1-getting-the-query-fasta-file"
  id="toc-1-getting-the-query-fasta-file">1 Getting the query fasta
  file</a>
- <a href="#2-database-creation" id="toc-2-database-creation">2 Database
  Creation</a>
  - <a href="#21-obtain-fasta-uniprotswiss-prot"
    id="toc-21-obtain-fasta-uniprotswiss-prot">2.1 Obtain Fasta
    (UniProt/Swiss-Prot)</a>
  - <a href="#22-making-the-database" id="toc-22-making-the-database">2.2
    Making the database</a>
- <a href="#3-running-blastx" id="toc-3-running-blastx">3 Running
  Blastx</a>
- <a href="#4-joining-blast-table-with-annoations"
  id="toc-4-joining-blast-table-with-annoations">4 Joining Blast table
  with annoations.</a>
  - <a href="#41-prepping-blast-table-for-easy-join"
    id="toc-41-prepping-blast-table-for-easy-join">4.1 Prepping Blast table
    for easy join</a>
  - <a href="#42-could-do-some-cool-stuff-in-r-here-reading-in-table"
    id="toc-42-could-do-some-cool-stuff-in-r-here-reading-in-table">4.2
    Could do some cool stuff in R here reading in table</a>

# 1 Getting the query fasta file

``` bash
curl https://gannet.fish.washington.edu/seashell/snaps/Gadus_macrocephalus.coding.gene.V1.cds \
-k \
> ../data/Gadus_macrocephalus.coding.gene.V1.cds
```

Exploring what fasta file

``` bash
head -3 ../data/Gadus_macrocephalus.coding.gene.V1.cds
```

    ## >Gma_1G0000010.1 locus=chr1:81612:97483:+    len:2343
    ## ATGCCTGTGAACGCGCGGGACCGGACAGTGCTGGGGCGTTTCCCCGGGGTCACGCTGGAA
    ## CCGGTGGAGGAGGAGGTGGAGGAGGAGGAGGAGGTGGAAGAGGACCAGGTGGAGCGAGGC

``` bash
echo "How many sequences are there?"
grep -c ">" ../data/Gadus_macrocephalus.coding.gene.V1.cds
```

    ## How many sequences are there?
    ## 23843

``` r
# Read FASTA file
fasta_file <- "../data/Gadus_macrocephalus.coding.gene.V1.cds"  # Replace with the name of your FASTA file
sequences <- readDNAStringSet(fasta_file)

# Calculate sequence lengths
sequence_lengths <- width(sequences)

# Create a data frame
sequence_lengths_df <- data.frame(Length = sequence_lengths)

# Plot histogram using ggplot2
ggplot(sequence_lengths_df, aes(x = Length)) +
  geom_histogram(binwidth = 1, color = "grey", fill = "blue", alpha = 0.75) +
  labs(title = "Histogram of Sequence Lengths",
       x = "Sequence Length",
       y = "Frequency") +
  theme_minimal()
```

<img src="03-transcriptome-annotation_files/figure-gfm/histogram-1.png" style="display: block; margin: auto;" />

``` r
# Read FASTA file
fasta_file <- "../data/Gadus_macrocephalus.coding.gene.V1.cds"
sequences <- readDNAStringSet(fasta_file)

# Calculate base composition
base_composition <- alphabetFrequency(sequences, baseOnly = TRUE)

# Convert to data frame and reshape for ggplot2
base_composition_df <- as.data.frame(base_composition)
base_composition_df$ID <- rownames(base_composition_df)
base_composition_melted <- reshape2::melt(base_composition_df, id.vars = "ID", variable.name = "Base", value.name = "Count")

# Plot base composition bar chart using ggplot2
ggplot(base_composition_melted, aes(x = Base, y = Count, fill = Base)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Base Composition",
       x = "Base",
       y = "Count") +
  theme_minimal() +
  scale_fill_manual(values = c("A" = "green", "C" = "blue", "G" = "yellow", "T" = "red"))
```

<img src="03-transcriptome-annotation_files/figure-gfm/ACGT-1.png" style="display: block; margin: auto;" />

``` r
# Read FASTA file
fasta_file <- "../data/Gadus_macrocephalus.coding.gene.V1.cds"
sequences <- readDNAStringSet(fasta_file)

# Count CG motifs in each sequence
count_cg_motifs <- function(sequence) {
  cg_motif <- "CG"
  return(length(gregexpr(cg_motif, sequence, fixed = TRUE)[[1]]))
}

cg_motifs_counts <- sapply(sequences, count_cg_motifs)

# Create a data frame
cg_motifs_counts_df <- data.frame(CG_Count = cg_motifs_counts)

# Plot CG motifs distribution using ggplot2
ggplot(cg_motifs_counts_df, aes(x = CG_Count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "blue", alpha = 0.75) +
  labs(title = "Distribution of CG Motifs",
       x = "Number of CG Motifs",
       y = "Frequency") +
  theme_minimal()
```

<img src="03-transcriptome-annotation_files/figure-gfm/cg-1.png" style="display: block; margin: auto;" />

# 2 Database Creation

## 2.1 Obtain Fasta (UniProt/Swiss-Prot)

This is from here picur reviewe sequences I named based on the identify
of the database given

``` bash
cd ../data
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta.gz uniprot_sprot_r2023_04.fasta.gz
gunzip -k uniprot_sprot_r2023_04.fasta.gz
```

## 2.2 Making the database

``` bash
mkdir ../blastdb
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/uniprot_sprot_r2023_04.fasta \
-dbtype prot \
-out ../blastdb/uniprot_sprot_r2023_04
```

# 3 Running Blastx

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastx \
-query ../data/Gadus_macrocephalus.coding.gene.V1.cds \
-db ../blastdb/uniprot_sprot_r2023_04 \
-out ../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx.tab \
-evalue 1E-20 \
-num_threads 20 \
-max_target_seqs 1 \
-outfmt 6
```

``` bash
head -2 ../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx.tab
```

    ## Gma_1G0000010.1  sp|P22735|TGM1_HUMAN    50.659  683 318 7   328 2334    109 786 0.0 688
    ## Gma_1G0000020.1  sp|Q9JI35|HRH3_CAVPO    54.684  395 160 4   136 1266    50  443 1.98e-140   411

``` bash
echo "Number of lines in output"
wc -l ../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx.tab
```

    ## Number of lines in output
    ## 11575 ../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx.tab

# 4 Joining Blast table with annoations.

## 4.1 Prepping Blast table for easy join

``` bash
tr '|' '\t' < ../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx.tab \
> ../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx_sep.tab

head -1 ../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx_sep.tab
```

    ## Gma_1G0000010.1  sp  P22735  TGM1_HUMAN  50.659  683 318 7   328 2334    109 786 0.0 688

## 4.2 Could do some cool stuff in R here reading in table

``` r
bltabl <- read.csv("../output/03-transcriptome-annotation/Gm.cds-uniprot_blastx_sep.tab", sep = '\t', header = FALSE)

spgo <- read.csv("https://gannet.fish.washington.edu/seashell/snaps/uniprot_table_r2023_01.tab", sep = '\t', header = TRUE)
```

``` r
datatable(head(bltabl), options = list(scrollX = TRUE, scrollY = "400px", scrollCollapse = TRUE, paging = FALSE))
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-789ac957d4247e05615e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-789ac957d4247e05615e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["Gma_1G0000010.1","Gma_1G0000020.1","Gma_1G0000030.1","Gma_1G0000040.1","Gma_1G0000070.1","Gma_1G0000090.1"],["sp","sp","sp","sp","sp","sp"],["P22735","Q9JI35","A9UMG5","Q9BXW6","Q6AXK4","P43136"],["TGM1_HUMAN","HRH3_CAVPO","IMPTB_XENTR","OSBL1_HUMAN","BABA1_DANRE","NR2F6_MOUSE"],[50.659,54.684,66.556,63.028,71.011,70.07299999999999],[683,395,302,971,376,411],[318,160,94,333,92,101],[7,4,3,8,4,3],[328,136,22,7,1,31],[2334,1266,924,2907,1122,1260],[109,50,6,2,1,1],[786,443,301,950,361,390],[0,1.98e-140,3.47e-132,0,1.51e-174,4.15e-160],[688,411,382,1193,493,459]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>V1<\/th>\n      <th>V2<\/th>\n      <th>V3<\/th>\n      <th>V4<\/th>\n      <th>V5<\/th>\n      <th>V6<\/th>\n      <th>V7<\/th>\n      <th>V8<\/th>\n      <th>V9<\/th>\n      <th>V10<\/th>\n      <th>V11<\/th>\n      <th>V12<\/th>\n      <th>V13<\/th>\n      <th>V14<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"scrollY":"400px","scrollCollapse":true,"paging":false,"columnDefs":[{"className":"dt-right","targets":[5,6,7,8,9,10,11,12,13,14]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"V1","targets":1},{"name":"V2","targets":2},{"name":"V3","targets":3},{"name":"V4","targets":4},{"name":"V5","targets":5},{"name":"V6","targets":6},{"name":"V7","targets":7},{"name":"V8","targets":8},{"name":"V9","targets":9},{"name":"V10","targets":10},{"name":"V11","targets":11},{"name":"V12","targets":12},{"name":"V13","targets":13},{"name":"V14","targets":14}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
datatable(head(spgo), options = list(scrollX = TRUE, scrollY = "400px", scrollCollapse = TRUE, paging = FALSE))
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-b69fe2c2333b55aa7069" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-b69fe2c2333b55aa7069">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6"],["A0A023I7E1","A0A024B7W1","A0A024SC78","A0A024SH76","A0A026W182","A0A044RE18"],["reviewed","reviewed","reviewed","reviewed","reviewed","reviewed"],["ENG1_RHIMI","POLG_ZIKVF","CUTI1_HYPJR","GUX2_HYPJR","ORCO_OOCBI","BLI_ONCVO"],["Glucan endo-1,3-beta-D-glucosidase 1 (Endo-1,3-beta-glucanase 1) (EC 3.2.1.39) (Laminarinase) (RmLam81A)","Genome polyprotein [Cleaved into: Capsid protein C (Capsid protein) (Core protein); Protein prM (Precursor membrane protein); Peptide pr (Peptide precursor); Small envelope protein M (Matrix protein); Envelope protein E; Non-structural protein 1 (NS1); Non-structural protein 2A (NS2A); Serine protease subunit NS2B (Flavivirin protease NS2B regulatory subunit) (Non-structural protein 2B); Serine protease NS3 (EC 3.4.21.91) (EC 3.6.1.15) (EC 3.6.4.13) (Flavivirin protease NS3 catalytic subunit) (Non-structural protein 3); Non-structural protein 4A (NS4A); Peptide 2k; Non-structural protein 4B (NS4B); RNA-directed RNA polymerase NS5 (EC 2.1.1.56) (EC 2.1.1.57) (EC 2.7.7.48) (NS5)]","Cutinase (EC 3.1.1.74)","Exoglucanase 2 (EC 3.2.1.91) (1,4-beta-cellobiohydrolase) (Cellobiohydrolase 6A) (Cel6A) (Exocellobiohydrolase II) (CBHII) (Exoglucanase II)","Odorant receptor coreceptor","Endoprotease bli (EC 3.4.21.75) (Blisterase)"],["ENG1 LAM81A","","M419DRAFT_76732","cbh2 M419DRAFT_122470","Orco X777_12371","Bli"],["Rhizomucor miehei","Zika virus (isolate ZIKV/Human/French Polynesia/10087PF/2013) (ZIKV)","Hypocrea jecorina (strain ATCC 56765 / BCRC 32924 / NRRL 11460 / Rut C-30) (Trichoderma reesei)","Hypocrea jecorina (strain ATCC 56765 / BCRC 32924 / NRRL 11460 / Rut C-30) (Trichoderma reesei)","Ooceraea biroi (Clonal raider ant) (Cerapachys biroi)","Onchocerca volvulus"],[796,3423,248,471,478,693],["glucan endo-1,3-beta-D-glucosidase activity [GO:0042973]; glucan endo-1,3-beta-glucanase activity, C-3 substituted reducing group [GO:0052861]; glucan endo-1,4-beta-glucanase activity, C-3 substituted reducing group [GO:0052862]","4 iron, 4 sulfur cluster binding [GO:0051539]; ATP binding [GO:0005524]; ATP hydrolysis activity [GO:0016887]; double-stranded RNA binding [GO:0003725]; exogenous protein binding [GO:0140272]; GTP binding [GO:0005525]; metal ion binding [GO:0046872]; mRNA (guanine-N7-)-methyltransferase activity [GO:0004482]; mRNA (nucleoside-2'-O-)-methyltransferase activity [GO:0004483]; protein dimerization activity [GO:0046983]; RNA helicase activity [GO:0003724]; RNA-dependent RNA polymerase activity [GO:0003968]; serine-type endopeptidase activity [GO:0004252]; structural molecule activity [GO:0005198]","cutinase activity [GO:0050525]","cellulose 1,4-beta-cellobiosidase activity [GO:0016162]; cellulose binding [GO:0030248]","odorant binding [GO:0005549]; olfactory receptor activity [GO:0004984]","metal ion binding [GO:0046872]; serine-type endopeptidase activity [GO:0004252]"],["extracellular region [GO:0005576]; glucan endo-1,3-beta-D-glucosidase activity [GO:0042973]; glucan endo-1,3-beta-glucanase activity, C-3 substituted reducing group [GO:0052861]; glucan endo-1,4-beta-glucanase activity, C-3 substituted reducing group [GO:0052862]; cell wall organization [GO:0071555]; polysaccharide catabolic process [GO:0000272]","extracellular region [GO:0005576]; host cell endoplasmic reticulum membrane [GO:0044167]; host cell nucleus [GO:0042025]; host cell perinuclear region of cytoplasm [GO:0044220]; membrane [GO:0016020]; viral capsid [GO:0019028]; viral envelope [GO:0019031]; virion membrane [GO:0055036]; 4 iron, 4 sulfur cluster binding [GO:0051539]; ATP binding [GO:0005524]; ATP hydrolysis activity [GO:0016887]; double-stranded RNA binding [GO:0003725]; exogenous protein binding [GO:0140272]; GTP binding [GO:0005525]; metal ion binding [GO:0046872]; mRNA (guanine-N7-)-methyltransferase activity [GO:0004482]; mRNA (nucleoside-2'-O-)-methyltransferase activity [GO:0004483]; protein dimerization activity [GO:0046983]; RNA helicase activity [GO:0003724]; RNA-dependent RNA polymerase activity [GO:0003968]; serine-type endopeptidase activity [GO:0004252]; structural molecule activity [GO:0005198]; clathrin-dependent endocytosis of virus by host cell [GO:0075512]; fusion of virus membrane with host endosome membrane [GO:0039654]; induction by virus of host autophagy [GO:0039520]; proteolysis [GO:0006508]; suppression by virus of host JAK-STAT cascade via inhibition of host TYK2 activity [GO:0039574]; suppression by virus of host JAK-STAT cascade via inhibition of STAT1 activity [GO:0039563]; suppression by virus of host JAK-STAT cascade via inhibition of STAT2 activity [GO:0039564]; suppression by virus of host transcription [GO:0039653]; suppression by virus of host type I interferon-mediated signaling pathway [GO:0039502]; viral RNA genome replication [GO:0039694]; virion attachment to host cell [GO:0019062]","extracellular region [GO:0005576]; cutinase activity [GO:0050525]","extracellular region [GO:0005576]; cellulose 1,4-beta-cellobiosidase activity [GO:0016162]; cellulose binding [GO:0030248]; cellulose catabolic process [GO:0030245]","plasma membrane [GO:0005886]; odorant binding [GO:0005549]; olfactory receptor activity [GO:0004984]; antennal development [GO:0007469]; detection of chemical stimulus involved in sensory perception of smell [GO:0050911]; detection of pheromone [GO:0043695]; olfactory behavior [GO:0042048]; response to pheromone [GO:0019236]; signal transduction [GO:0007165]; social behavior [GO:0035176]","extracellular region [GO:0005576]; metal ion binding [GO:0046872]; serine-type endopeptidase activity [GO:0004252]; dibasic protein processing [GO:0090472]; zymogen activation [GO:0031638]"],["cell wall organization [GO:0071555]; polysaccharide catabolic process [GO:0000272]","clathrin-dependent endocytosis of virus by host cell [GO:0075512]; fusion of virus membrane with host endosome membrane [GO:0039654]; induction by virus of host autophagy [GO:0039520]; proteolysis [GO:0006508]; suppression by virus of host JAK-STAT cascade via inhibition of host TYK2 activity [GO:0039574]; suppression by virus of host JAK-STAT cascade via inhibition of STAT1 activity [GO:0039563]; suppression by virus of host JAK-STAT cascade via inhibition of STAT2 activity [GO:0039564]; suppression by virus of host transcription [GO:0039653]; suppression by virus of host type I interferon-mediated signaling pathway [GO:0039502]; viral RNA genome replication [GO:0039694]; virion attachment to host cell [GO:0019062]","","cellulose catabolic process [GO:0030245]","antennal development [GO:0007469]; detection of chemical stimulus involved in sensory perception of smell [GO:0050911]; detection of pheromone [GO:0043695]; olfactory behavior [GO:0042048]; response to pheromone [GO:0019236]; signal transduction [GO:0007165]; social behavior [GO:0035176]","dibasic protein processing [GO:0090472]; zymogen activation [GO:0031638]"],["extracellular region [GO:0005576]","extracellular region [GO:0005576]; host cell endoplasmic reticulum membrane [GO:0044167]; host cell nucleus [GO:0042025]; host cell perinuclear region of cytoplasm [GO:0044220]; membrane [GO:0016020]; viral capsid [GO:0019028]; viral envelope [GO:0019031]; virion membrane [GO:0055036]","extracellular region [GO:0005576]","extracellular region [GO:0005576]","plasma membrane [GO:0005886]","extracellular region [GO:0005576]"],["GO:0000272; GO:0005576; GO:0042973; GO:0052861; GO:0052862; GO:0071555","GO:0003724; GO:0003725; GO:0003968; GO:0004252; GO:0004482; GO:0004483; GO:0005198; GO:0005524; GO:0005525; GO:0005576; GO:0006508; GO:0016020; GO:0016887; GO:0019028; GO:0019031; GO:0019062; GO:0039502; GO:0039520; GO:0039563; GO:0039564; GO:0039574; GO:0039653; GO:0039654; GO:0039694; GO:0042025; GO:0044167; GO:0044220; GO:0046872; GO:0046983; GO:0051539; GO:0055036; GO:0075512; GO:0140272","GO:0005576; GO:0050525","GO:0005576; GO:0016162; GO:0030245; GO:0030248","GO:0004984; GO:0005549; GO:0005886; GO:0007165; GO:0007469; GO:0019236; GO:0035176; GO:0042048; GO:0043695; GO:0050911","GO:0004252; GO:0005576; GO:0031638; GO:0046872; GO:0090472"],["","","","","",""],["3.2.1.39","2.1.1.56; 2.1.1.57; 2.7.7.48; 3.4.21.91; 3.6.1.15; 3.6.4.13","3.1.1.74","3.2.1.91","","3.4.21.75"],["","","","","",""],["","","","","",""],["IPR005200;IPR040720;IPR040451;","IPR011492;IPR043502;IPR000069;IPR038302;IPR013755;IPR001122;IPR037172;IPR027287;IPR026470;IPR038345;IPR001157;IPR000752;IPR000487;IPR000404;IPR001528;IPR046811;IPR002535;IPR038688;IPR000336;IPR001850;IPR014412;IPR011998;IPR036253;IPR038055;IPR013756;IPR014001;IPR001650;IPR014756;IPR026490;IPR027417;IPR009003;IPR000208;IPR007094;IPR002877;IPR029063;","IPR029058;IPR000675;IPR043580;IPR043579;IPR011150;","IPR016288;IPR036434;IPR035971;IPR000254;IPR001524;","IPR004117;","IPR008979;IPR034182;IPR002884;IPR000209;IPR036852;IPR023827;IPR022398;IPR023828;IPR015500;IPR032815;IPR038466;"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Entry<\/th>\n      <th>Reviewed<\/th>\n      <th>Entry.Name<\/th>\n      <th>Protein.names<\/th>\n      <th>Gene.Names<\/th>\n      <th>Organism<\/th>\n      <th>Length<\/th>\n      <th>Gene.Ontology..molecular.function.<\/th>\n      <th>Gene.Ontology..GO.<\/th>\n      <th>Gene.Ontology..biological.process.<\/th>\n      <th>Gene.Ontology..cellular.component.<\/th>\n      <th>Gene.Ontology.IDs<\/th>\n      <th>Interacts.with<\/th>\n      <th>EC.number<\/th>\n      <th>Reactome<\/th>\n      <th>UniPathway<\/th>\n      <th>InterPro<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"scrollY":"400px","scrollCollapse":true,"paging":false,"columnDefs":[{"className":"dt-right","targets":7},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Entry","targets":1},{"name":"Reviewed","targets":2},{"name":"Entry.Name","targets":3},{"name":"Protein.names","targets":4},{"name":"Gene.Names","targets":5},{"name":"Organism","targets":6},{"name":"Length","targets":7},{"name":"Gene.Ontology..molecular.function.","targets":8},{"name":"Gene.Ontology..GO.","targets":9},{"name":"Gene.Ontology..biological.process.","targets":10},{"name":"Gene.Ontology..cellular.component.","targets":11},{"name":"Gene.Ontology.IDs","targets":12},{"name":"Interacts.with","targets":13},{"name":"EC.number","targets":14},{"name":"Reactome","targets":15},{"name":"UniPathway","targets":16},{"name":"InterPro","targets":17}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
datatable(
  left_join(bltabl, spgo,  by = c("V3" = "Entry")) %>%
  select(V1, V3, V13, Protein.names, Organism, Gene.Ontology..biological.process., Gene.Ontology.IDs) 
 # %>% mutate(V1 = str_replace_all(V1,pattern = "solid0078_20110412_FRAG_BC_WHITE_WHITE_F3_QV_SE_trimmed", replacement = "Ab"))
)
```

<div class="datatables html-widget html-fill-item" id="htmlwidget-3caca9a71bb05ed325bc" style="width:100%;height:auto;"></div>

``` r
annot_tab <-
  left_join(bltabl, spgo,  by = c("V3" = "Entry")) %>%
  select(V1, V3, V13, Protein.names, Organism, Gene.Ontology..biological.process., Gene.Ontology.IDs)

write.table(annot_tab, file = "../output/03-transcriptome-annotation/G_macrocephalus_IDmapping_2024_04_17.tab", sep = "\t",
            row.names = TRUE, col.names = NA)
```

``` bash
head -n 3 ../output/03-transcriptome-annotation/G_macrocephalus_IDmapping_2024_04_17.tab
```

``` r
# Read dataset
#dataset <- read.csv("../output/blast_annot_go.tab", sep = '\t')  # Replace with the path to your dataset

# Select the column of interest
column_name <- "Organism"  # Replace with the name of the column of interest
column_data <- annot_tab[[column_name]]

# Count the occurrences of the strings in the column
string_counts <- table(column_data)

# Convert to a data frame, sort by count, and select the top 10
string_counts_df <- as.data.frame(string_counts)
colnames(string_counts_df) <- c("String", "Count")
string_counts_df <- string_counts_df[order(string_counts_df$Count, decreasing = TRUE), ]
top_10_strings <- head(string_counts_df, n = 10)

# Plot the top 10 most common strings using ggplot2
ggplot(top_10_strings, aes(x = reorder(String, -Count), y = Count, fill = String)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Top 10 Species hits",
       x = column_name,
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip()
```

<img src="03-transcriptome-annotation_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

``` r
#data <- read.csv("../output/blast_annot_go.tab", sep = '\t')

# Rename the `Gene.Ontology..biological.process.` column to `Biological_Process`
colnames(annot_tab)[colnames(annot_tab) == "Gene.Ontology..biological.process."] <- "Biological_Process"

# Separate the `Biological_Process` column into individual biological processes
data_separated <- unlist(strsplit(annot_tab$Biological_Process, split = ";"))

# Trim whitespace from the biological processes
data_separated <- gsub("^\\s+|\\s+$", "", data_separated)

# Count the occurrences of each biological process
process_counts <- table(data_separated)
process_counts <- data.frame(Biological_Process = names(process_counts), Count = as.integer(process_counts))
process_counts <- process_counts[order(-process_counts$Count), ]

# Select the 20 most predominant biological processes
top_20_processes <- process_counts[1:20, ]

# Create a color palette for the bars
bar_colors <- rainbow(nrow(top_20_processes))

# Create a staggered vertical bar plot with different colors for each bar
barplot(top_20_processes$Count, names.arg = rep("", nrow(top_20_processes)), col = bar_colors,
        ylim = c(0, max(top_20_processes$Count) * 1.25),
        main = "Occurrences of the 20 Most Predominant Biological Processes", xlab = "Biological Process", ylab = "Count")
```

<img src="03-transcriptome-annotation_files/figure-gfm/go-1.png" style="display: block; margin: auto;" />

``` r
# Create a separate plot for the legend
png("../output/GOlegend.png", width = 800, height = 600)
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = top_20_processes$Biological_Process, fill = bar_colors, cex = 1, title = "Biological Processes")
dev.off()
```

    ## png 
    ##   2

``` r
knitr::include_graphics("../output/GOlegend.png")
```

<img src="../output/GOlegend.png" width="800" style="display: block; margin: auto;" />