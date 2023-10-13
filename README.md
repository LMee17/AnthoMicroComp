# Comparing Micobial Community Compositions Across Anthophila

## Step 1: Transcriptomic Sample Processing

All scripts in this section are found in `Scripts/Bash/`. 

### 1.1 Non-_Apis mellifera_ Samples
CZID.org (Kalantar _et al_ 2020) has different "Host" genome options that the user selects to map their nucleotide sequences against, with the end goal of mining the unmapped reads for microbial identification. CZID has a "Bee" option on its webpage at the time of writing (February 2023) which aligns given reads against the _Apis mellifera_ genome. For the many non-_A.mellifera_ samples in this study there were a number of preprocessing steps that had to take place. The script `Process_PreMapper_2.sh` runs these processes which essentially consists of:

1. Downloading the given SRA sample using `SRA Toolkit` (Katz _et al_ 2022)
2. Unpacking the sample
3. Mapping it against the nearest host genome using `STAR` (Dobin _et al_ 2013)
4. Retaining unmapped reads
5. Putting reads that mapped at < 50% to the side for further processing

This script requires a list of SRA sample accessions (one per line) and a directory of genomes as input, as well as the above packages installed, to work.

The genomes used in this analysis were: 

- Andrena spp. : GCA_929108735.1
- Andrena camellia : GCA_910592295.1
- Andrena cineraria : GCA_944738655.1
- Andrena fulva : GCA_910592295.1
- Andrena haemorrhoa : GCA_910592295.1
- Andrena vaga : GCA_944738655.1
- Anthophora plumipes : GCA_001263275.1
- Apis cerana : GCA_001442555.1
- Bombus spp. : GCA_910591885.2
- Bombus breviceps : GCA_014825925.1
- Bombus confusus : GCA_014737475.1
- Bombus consobrinus : GCA_014737475.1
- Bombus difficillimus : GCA_014737525.1
- Bombus haemorrhoidalis : GCA_014825975.1
- Bombus ignitus : GCA_014825875.1
- Bombus lucorum : GCA_910591885.2
- Bombus opulentus : GCA_014737405.1
- Bombus pascuorum : GCA_905332965.1
- Bombus picipes : GCA_014737485.1
- Bombus pyrosoma : GCA_014825855.1
- Bombus rupestris : GCA_014737355.1
- Bombus sibiricus : GCA_014737505.1
- Bombus soroeensis : GCA_014737365.1
- Bombus superbus : GCA_014737385.1
- Bombus terrestris : GCA_910591885.2
- Bombus terricola : GCA_910591885.2
- Bombus turneri : GCA_014825825.1
- Bombus waltoni : GCA_014737395.1
- Ceratina australensis : GCA_004307685.1
- Dufourea novaeangliae : GCA_001272555.1
- Epeolus variegatus : GCA_907165295.1
- Eufriesea mexicana : GCA_001483705.1
- Euglossa dilemma : GCA_001483705.1
- Euglossa viridissima : GCA_001483705.1
- Exoneura spp. : GCA_019453415.1
- Habropoda laboriosa : GCA_001263275.1
- Halictus sexcinctus : GCA_916610255.1
- Lasioglossum spp. : GCA_916610255.1
- Megalopta genalis : GCA_011865705.1
- Nomada lathburiana : GCA_907165295.1
- Osmia bicornis : GCA_907164935.1
- Osmia cornuta : GCA_907164935.1
- Tetragonisca angustula : GCA_011634685.1
- Tetragonula carbonaria : GCA_011634685.1

Any samples that had less than 50% reads mapped were re-mapped using more relaxed parameters and `ReProcess.sh`. All unmapped reads were uploaded to CZID using `CZidUpload.sh`, which requires the installation of the CZID command-line interface (`CLI`).

### 1.2 _A. mellifera_ Samples
As there is already an _A. mellifera_ genome available on CZID.org, then these samples were downloaded from the SRA, unpacked and uploaded to CZID using the single script `ProcessSRA_v3.sh`. This requires similar input and the same dependencies as the above scripts.

### 1.3 CZID Pipeline
The CZID pipeline runs automatically from the moment samples are uploaded (see the website for details / instructions). Once all are complete, taxon reports can be downloaded in bulk for the next stage of the anlaysis.

## Step Two: Community Composition Analyses

All of these scripts are ran in `R` and can be found in `Scripts/R`. Required `R` packages are indicated in the "Libraries" sections of the scripts at the outset.

Metadata is found in `input/Metadata/`, taxonomic classifications in `input/Phylo_Misc/`, CZID statistics in `input/CZ_Overview/`, CZID taxon reports in `input/TaxonReports/` and counts are saved in `input/Counts/` once produced.

### 2.1 Taxonomy

All taxonomic classifications in this analysis were taken from a `.dmp` file from NCBI Taxonomy accessed 18th October 2022 (Schoch _et al_ 2020). The taxonomy used throughout this analysis was made using `Scripts/Bash/SplitTax.sh` and `MakeTaxa.R` (instructions are within the annotations of the latter). 

### 2.2 Concatenating and Filtering Data

`ReadnFilt_Dec22.R` takes a folder of CZID taxon report files as input and produces count tables. This is tailored specifically for these data and requires the prefiltered sample metadata and a text file on the discrepancies between the taxon report taxa and the current taxonomy (according to NCBI, October 2022). As output it produces count tables with samples as columns and microbial hit as rows. This is done for 1) Prokaryotes 2) Eukaryotes and 3) Viruses. 

### 2.3 Community Composition Anlaysis

`Analysis_Prokaryote3.R` uses `vegan` (Dixon _et al_ 2003) packages to compute dissimilarity matrices, produce NMDS and use PERMANOVA (`adonis2`) to assess the impact of phylogeny, location and social lifestyle on prokaryote (bacterial) community composition. `Analysis_Eukaryote3.R` and `Analysis_Viral3.R` do the same for eukaryotic and viral counts, respectively.

### 2.4 Assessing the Presence of Well Characterised Bee Microbiota and Tribe-Microbe Associations

To see if there are candidates for microbes associated with particular tribes, `DescribeTribe.R` identifies microbes found in over 50% of tribal members at an average relative abundance above 0.01%. 

## Program Versions

`CZID CLI` v4.1.2 \
`CZID Pipeline` v7.1 \
`R` v4.2.2 \
`SRA Toolkit` v3.0.0 \
`vegan` v2.6-4  

## References

Dixon, P., 2003. VEGAN, a package of R functions for community ecology. _Journal of vegetation science_, 14(6), pp.927-930.\
Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M. and Gingeras, T.R., 2013. STAR: ultrafast universal RNA-seq aligner. _Bioinformatics_, 29(1), pp.15-21.\
Kalantar, K.L., Carvalho, T., de Bourcy, C.F., Dimitrov, B., Dingle, G., Egger, R., Han, J., Holmes, O.B., Juan, Y.F., King, R. and Kislyuk, A., 2020. IDseq—An open source cloud-based pipeline and analysis service for metagenomic pathogen detection and monitoring. _Gigascience_, 9(10), p.giaa111. \
Katz, K., Shutov, O., Lapoint, R., Kimelman, M., Brister, J.R. and O’Sullivan, C., 2022. The sequence read archive: a decade more of explosive growth. _Nucleic acids research_, 50(D1), pp.D387-D390.\
Schoch, C.L., Ciufo, S., Domrachev, M., Hotton, C.L., Kannan, S., Khovanskaya, R., Leipe, D., Mcveigh, R., O’Neill, K., Robbertse, B. and Sharma, S., 2020. NCBI Taxonomy: a comprehensive update on curation, resources and tools. _Database_, 2020.
