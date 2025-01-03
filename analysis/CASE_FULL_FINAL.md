CASE Full FINAL
================

# Setup

``` r
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(qvalue)
library(stringr)
library(pcadapt)
library(poolSeq)
```

    ## Loading required package: data.table

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## Loading required package: foreach

    ## Loading required package: stringi

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## The following object is masked from 'package:plyr':
    ## 
    ##     count

    ## Loading required package: Rcpp

``` r
library(ACER)
library(ggman)
```

    ## Loading required package: ggrepel

``` r
library(VennDiagram)
```

    ## Loading required package: grid

    ## Loading required package: futile.logger

``` r
library(scales)
library(data.table)
library(patchwork)


# Custom theme for black background figures
theme_black = function(base_size = 16, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "white", size = 1.5),
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "transparent"),  
      legend.key = element_rect(color = "white",  fill = "transparent"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "transparent", color  =  NA),  
      panel.border = element_blank(),  
      panel.grid.major = element_line(color = "transparent"),  
      panel.grid.minor = element_line(color = "transparent"),  
      #panel.margin = unit(0.5, "lines"),   
      panel.spacing= unit(0.5, "lines"),
      # Specify facetting options
      strip.background = element_rect(fill = "transparent", color = "transparent"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "transparent", fill = "transparent"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

# Function to convert p-values to q-values
Qvalue_convert <- function(table) {


table <- spread(table,GROUP,PVAL)  
  
table$CHR <- table$CHROM

table %>% 
  mutate(CHR = str_replace(CHR, "NC_035780.1", "1")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035781.1", "2")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035782.1", "3")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035783.1", "4")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035784.1", "5")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035785.1", "6")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035786.1", "7")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035787.1", "8")) %>%
  mutate(CHR = str_replace(CHR, "NC_035788.1", "9")) %>% 
  mutate(CHR = str_replace(CHR, "NC_035789.1", "10"))  -> pp1

pp1 <- unite(pp1, "SNP", c("CHROM","BP"),remove = FALSE)

pp1$CHR <- as.numeric(pp1$CHR)  

# Set pi0=1 to equal Benjamini Hochberg
pp1$QCA <- qvalue(pp1$PCA, pi0 = 1)$qvalues
pp1$QCON <- qvalue(pp1$PCON, pi0 = 1)$qvalues
pp1$QSE <- qvalue(pp1$PSE, pi0 = 1)$qvalues
pp1$QCASE <- qvalue(pp1$PCASE, pi0 = 1)$qvalues

#pp1$QCA <- qvalue(pp1$PCA)$qvalues
#pp1$QCON <- qvalue(pp1$PCON)$qvalues
#pp1$QSE <- qvalue(pp1$PSE)$qvalues
#pp1$QCASE <- qvalue(pp1$PCASE)$qvalues

return(pp1)
}

# Function to determine significant loci in various outpus
Significat_subset <- function(pv, alpha, alpha2) {
  # Add new logical columns based on conditions
  pv$Sig.CASE <- pv$QCASE < alpha
  pv$Sig.CA <- pv$QCA < alpha
  pv$Sig.SE <- pv$QSE < alpha
  
  # Subset the dataframe
  ppsig <- subset(pv, QCA < alpha | QCASE < alpha | QSE < alpha)
  ppsig <- subset(ppsig, QCON > alpha2)
  
  # Print diagnostics
  print(nrow(pv))
  print(nrow(ppsig))
  print(nrow(ppsig) / nrow(pv))
  
  # Return the modified dataframe
  return(ppsig)
}

group_and_average <- function(df) {
  df %>%
    dplyr::group_by(SNP) %>%
    dplyr::summarize(
      # Use a conditional inside across to handle all-NA groups
      across(
        c(PCA, PCASE, PCON, PSE, CHR, QCA, QCON, QSE, QCASE), 
        ~ if (all(is.na(.))) NA else min(., na.rm = TRUE)),
      # Combine logical columns, retain TRUE if any TRUE exists in the group
      Sig.CASE = any(Sig.CASE, na.rm = TRUE),
      Sig.CA = any(Sig.CA, na.rm = TRUE),
      Sig.SE = any(Sig.SE, na.rm = TRUE),
      CHROM = first(CHROM),
      BP = first(BP),
      .groups = "drop"  # Ungroup the result
    )
}


cbPaletteSmall <- c("#E69F00", "#009E73", "#999999","sky blue","#0072B2")
cbPaletteSmall4 <- c("#999999","#E69F00", "#0072B2", "#009E73")
cbPaletteSmall3 <- c("#E69F00", "#0072B2","#009E73")
```

``` r
sessionInfo()
```

    ## R version 3.6.0 (2019-04-26)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Scientific Linux 7.4 (Nitrogen)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/local/lib64/R/lib/libRblas.so
    ## LAPACK: /usr/local/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.1.1     scales_1.3.0        VennDiagram_1.7.0  
    ##  [4] futile.logger_1.4.3 ggman_0.99.0        ggrepel_0.9.5.9999 
    ##  [7] ACER_1.0            poolSeq_0.3.5       Rcpp_1.0.12        
    ## [10] matrixStats_0.57.0  stringi_1.7.8       foreach_1.5.0      
    ## [13] data.table_1.14.2   pcadapt_4.3.3       stringr_1.4.0      
    ## [16] qvalue_2.18.0       dplyr_1.0.7         plyr_1.8.6         
    ## [19] tidyr_1.1.4         ggplot2_3.5.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtools_3.8.2         tidyselect_1.1.1     xfun_0.31           
    ##  [4] purrr_1.0.2          reshape2_1.4.4       splines_3.6.0       
    ##  [7] colorspace_2.1-0     vctrs_0.6.5          generics_0.1.1      
    ## [10] htmltools_0.5.2      yaml_2.2.1           utf8_1.2.4          
    ## [13] blob_1.2.1           rlang_1.1.3          pillar_1.9.0        
    ## [16] glue_1.7.0           withr_3.0.0          DBI_1.1.0           
    ## [19] lambda.r_1.2.4       lifecycle_1.0.4      munsell_0.5.0       
    ## [22] gtable_0.3.4         codetools_0.2-16     evaluate_0.15       
    ## [25] knitr_1.39           fastmap_1.1.0        fansi_1.0.6         
    ## [28] formatR_1.11         digest_0.6.29        cli_3.6.2           
    ## [31] tools_3.6.0          magrittr_2.0.3       tibble_3.2.1        
    ## [34] futile.options_1.0.1 pkgconfig_2.0.3      assertthat_0.2.1    
    ## [37] rmarkdown_2.12       rstudioapi_0.13      iterators_1.0.12    
    ## [40] R6_2.5.1             compiler_3.6.0

## Download Popoolation2 scripts

``` bash
cd ../scripts
git clone https://github.com/ToBoDev/assessPool.git
```

## Create conda (mamba) environment

``` bash
mamba env create --file ../CASE_environment.yaml
mamba env create --file ../random_draw_environment.yaml
```

# Filtering

``` bash
source activate CASE
cd ../raw.vcf
bash ../scripts/CASE_PVCF_filter2.sh CASE 64 20
```

``` bash
source activate CASE
cd ../raw.vcf 
bcftools index CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz
bcftools query -l CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz > samples
```

## Break up into Spawns (Blocks)

Each sequencing run had some extra samples that need to be removed
before starting the analysis Each block was filtered for missing data
less than 10& and MAF of \> 0.015 based on read counts

``` bash
source activate CASE

cd ../raw.vcf
bcftools view -S <(grep B11 samples | grep -v 1.G) CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz | bcftools view -i 'F_MISSING<0.1'  | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' -O z -o B11.CASE.FIL.vcf.gz

bcftools index B11.CASE.FIL.vcf.gz

bcftools query -f '%CHROM\t%POS\n' B11.CASE.FIL.vcf.gz > B11.pos

bcftools view -S <(grep B10 samples | grep -v J17 | grep -v J05) CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz | bcftools view -i 'F_MISSING<0.1' | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' -O z -o B10.CASE.FIL.vcf.gz
bcftools index B10.CASE.FIL.vcf.gz
bcftools query -f '%CHROM\t%POS\n' B10.CASE.FIL.vcf.gz > B10.pos

bcftools view -S <(grep B12 samples) CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz | bcftools view -i 'F_MISSING<0.1' | bcftools view -M 4 -m 2 | bcftools +fill-tags -- -t 'AAF:1=sum(FORMAT/AO)/sum(FORMAT/DP)' | bcftools +fill-tags -- -t  FORMAT/VAF | bcftools view --threads 40 -i 'AAF > 0.015 && AAF < 0.985' -O z -o B12.CASE.FIL.vcf.gz
bcftools index B12.CASE.FIL.vcf.gz
bcftools query -f '%CHROM\t%POS\n' B12.CASE.FIL.vcf.gz > B12.pos

cat B*.pos | sort | uniq > fil.pos

bcftools view -R fil.pos --threads 40 -m2 -M 4 CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz -O z -o SNP.CASE.TRSdp.20.B90.2a.perp.vcf.gz

bcftools view --threads 40 SNP.CASE.TRSdp.20.B90.2a.perp.vcf.gz | mawk '!/#/' | cut -f 1,2 | mawk '{print $1"\t"$2-1"\t"$2}' > total.snp.bed
#bedtools intersect -wb -a total.snp.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > CASE.study.background.LOC
```

### Create Sync files

``` bash
source activate CASE

bcftools view --threads 40 ../raw.vcf/B10.CASE.FIL.vcf.gz | mawk '!/\.:\.:\./' > temp.vcf
python2 ../scripts/VCFtoPopPool.py temp.vcf CASE.Block10.sync 
rm temp.vcf

bcftools view --threads 40 ../raw.vcf/B11.CASE.FIL.vcf.gz | mawk '!/\.:\.:\./' > temp.vcf
python2 ../scripts/VCFtoPopPool.py temp.vcf CASE.Block11.sync 
rm temp.vcf

bcftools view --threads 40 ../raw.vcf/B12.CASE.FIL.vcf.gz  | mawk '!/\.:\.:\./' > temp.vcf
python2 ../scripts/VCFtoPopPool.py temp.vcf CASE.Block12.sync 
rm temp.vcf

bcftools view --threads 40 ../raw.vcf/SNP.CASE.TRSdp.20.B90.2a.perp.vcf.gz  | mawk '!/\.:\.:\./' > temp.vcf
python2 ../scripts/VCFtoPopPool.py temp.vcf CASE.All.Blocks.sync 
rm temp.vcf
```

## Add coverage stats to sync file and filter by minimum coverage

``` bash
source activate CASE
mawk -f ../scripts/add_cov_sync CASE.Block10.sync | mawk '$20 > 9 && $22 > 24'> CASE.dp20.Block10.cov.sync &
mawk -f ../scripts/add_cov_sync CASE.Block11.sync | mawk '$24 > 9 && $26 > 24'> CASE.dp20.Block11.cov.sync &
mawk -f ../scripts/add_cov_sync CASE.Block12.sync | mawk '$24 > 9 && $26 > 24'> CASE.dp20.Block12.cov.sync &
mawk -f ../scripts/add_cov_sync CASE.All.Blocks.sync | mawk '$66 > 9 && $68 > 24' > CASE.All.Blocks.cov.sync
```

# Visualize across 10,000 random loci

## B10

``` bash
source activate CASE
timeout 1s mawk '!/CHR/' CASE.dp20.Block10.cov.sync | cut -f 1-19 | shuf | shuf | head -q -n 11000 > input.b10.rand.sync 
../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input input.b10.rand.sync --output-prefix B10.dp20.rand --max-coverage 50000 

mawk -f ../scripts/polarize_freqs B10.dp20.rand_rc | grep -Ei '(A{16}|C{16}|G{16}|T{16})'| head -q -n 10000| cut -f10-25 | mawk '!/maa/'   > B10.dp20.rand.pool 
```

``` r
# Read data from a file named "B10.dp20.rand.pool" into a data frame.
pool.data.b10 <- read.table("B10.dp20.rand.pool")

# Apply a function to evaluate each cell of the data frame as an R expression.
# This converts the text content of the cells into their evaluated form.
df.pool.b10 <- apply(pool.data.b10, c(1, 2), function(x) eval(parse(text = x)))

# Filter rows from the data frame where all elements in a row are the same.
# The resulting data frame only keeps rows with at least one differing value.
df_filtered <- df.pool.b10[!apply(df.pool.b10, 1, function(row) all(row == row[1])), ]

# Transpose the filtered data to prepare it for input into the PCAdapt package.
pool.data2.b10 <- t(df_filtered)

# Read the transposed data and prepare it for PCAdapt analysis, specifying that the data type is "pool".
filename.b10 <- read.pcadapt(pool.data2.b10, type = "pool")

# Perform a PCAdapt analysis to detect population structure using 4 principal components (K=5)
# and a minimum minor allele frequency (MAF) threshold of 0.01.
res.b10.rand <- pcadapt(filename.b10, K = 4, min.maf = 0.01)

# Define a vector of population labels, where each label corresponds to samples in the dataset.
poplist.names <- c(rep("CA", 3), rep("CASE", 3), rep("CON", 3), rep("IS", 4), rep("SE", 3))

# Extract the PC scores from the PCAdapt result and store them in a data frame.
p1.b10.rand.df <- data.frame(res.b10.rand$scores)

# Rename the columns of the data frame to indicate the principal components (PC1, PC2, etc.).
colnames(p1.b10.rand.df) <- c("PC1", "PC2", "PC3", "PC4")

# Add a column to the data frame for population labels.
p1.b10.rand.df$POP <- poplist.names

# Create a scatter plot of the first two principal components (PC1 vs. PC2).
# Points are filled based on population (POP) and have a transparent appearance.
ggplot(p1.b10.rand.df, aes(x = PC1, y = PC2, fill = POP)) + 
  geom_point(aes(alpha = 0.2), size = 10, shape = 21, col = "black") + 
  scale_fill_manual(values = cbPaletteSmall, name = "Treatment") + 
  guides(alpha = FALSE) + 
  theme(axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24))
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
png(filename="PC_B10_rand.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.b10.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

## B11

``` bash
source activate CASE

mawk '!/CHR/' CASE.dp20.Block11.cov.sync | cut -f1-23 | shuf | shuf -n 11000 > input.b11.rand.sync
../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input input.b11.rand.sync --output-prefix B11.dp20.rand --max-coverage 50000 

mawk -f ../scripts/polarize_freqs B11.dp20.rand_rc | grep -Ei '(A{20}|C{20}|G{20}|T{20})' | mawk '!/maa/' | cut -f10-29 | head -q -n 10000 > B11.dp20.rand.pool
```

``` r
pool.data.b11 <- read.table("B11.dp20.rand.pool")
df.pool.b11 <- apply(pool.data.b11, c(1, 2), function(x) eval(parse(text = x)))
df_filtered <- df.pool.b11[!apply(df.pool.b11, 1, function(row) all(row == row[1])), ]

pool.data2.b11 <- t(df_filtered)
filename.b11 <- read.pcadapt(pool.data2.b11, type = "pool")

res.b11.rand <- pcadapt(filename.b11, K =5,min.maf = 0.01)

poplist.names <- c(rep("CA", 4),rep("CASE", 4),rep("CON", 4),rep("IS", 4),rep("SE", 4))
p1.b11.rand.df <- data.frame(res.b11.rand$scores)
colnames(p1.b11.rand.df) <- c("PC1","PC2", "PC3", "PC4")
p1.b11.rand.df$POP <- poplist.names
ggplot(p1.b11.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="black")+  scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
png(filename="PC_B11_rand.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.b11.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

## B12

``` bash
source activate CASE

mawk '!/CHR/' CASE.dp20.Block12.cov.sync | cut -f1-23 > temp1
shuf temp1 | shuf -n 11000  > input.b12.rand.sync
rm temp1
../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input input.b12.rand.sync --output-prefix B12.dp20.rand --max-coverage 50000 

mawk -f ../scripts/polarize_freqs B12.dp20.rand_rc |  grep -Ei '(A{20}|C{20}|G{20}|T{20})' -m 10000 | cut -f10-29  > B12.dp20.rand.pool
```

``` r
pool.data.b12 <- read.table("B12.dp20.rand.pool")
df.pool.b12 <- apply(pool.data.b12, c(1, 2), function(x) eval(parse(text = x)))
df_filtered <- df.pool.b12[!apply(df.pool.b12, 1, function(row) all(row == row[1])), ]

pool.data2.b12 <- t(df_filtered)
filename.b12 <- read.pcadapt(pool.data2.b12, type = "pool")

res.b12.rand <- pcadapt(filename.b12, K =5,min.maf = 0.01)

poplist.names <- c(rep("CA", 4),rep("CASE", 4),rep("CON", 4),rep("IS", 4),rep("SE", 4))
p1.b12.rand.df <- data.frame(res.b12.rand$scores)
colnames(p1.b12.rand.df) <- c("PC1","PC2", "PC3", "PC4")
p1.b12.rand.df$POP <- poplist.names
ggplot(p1.b12.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="black")+  scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
png(filename="PC_B12_rand.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.b12.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

## All

``` bash
source activate CASE

mawk '!/CHR/' CASE.All.Blocks.cov.sync | mawk '$66 > 9' | cut --complement -f15,4,17,29,40,65- | shuf -n 50000 | shuf -n 11000  > input.all.rand.sync
../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input input.all.rand.sync --output-prefix all.dp20.rand --max-coverage 50000 

mawk -f ../scripts/polarize_freqs all.dp20.rand_rc | grep -Ei '(A{20}|C{20}|G{20}|T{20})' -m 10000 | cut -f10-65  > all.dp20.rand.pool
```

``` r
pool.data.all <- read.table("all.dp20.rand.pool")
df.pool.all <- apply(pool.data.all, c(1, 2), function(x) eval(parse(text = x)))
df_filtered.all <- df.pool.all[!apply(df.pool.all, 1, function(row) all(row == row[1])), ]

pool.data2.all <- t(df_filtered.all)
filename.all <- read.pcadapt(pool.data2.all, type = "pool")

res.all.rand <- pcadapt(filename.all, K =5,min.maf = 0.0001)

poplist.names <- c(rep("CA", 11),rep("CASE", 11),rep("CON", 11),rep("IS", 12),rep("SE", 11))
poplist.names <- c("B12","B11","B12","B10","B11","B10","B11","B12","B10","B12","B11","B10","B11","B12","B12","B10","B12","B11","B10","B11","B11","B12","B10","B11","B12","B11","B10","B12","B10","B11","B12","B11","B12","B10","B11","B10","B11","B10","B11","B10","B11","B12","B12","B12","B12","B11","B12","B10","B10","B11","B11","B10","B12","B12","B11","B12")
p1.all.rand.df <- data.frame(res.all.rand$scores)
colnames(p1.all.rand.df) <- c("PC1","PC2", "PC3", "PC4")
p1.all.rand.df$POP <- poplist.names
ggplot(p1.all.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="black")+  scale_fill_manual(values=cbPaletteSmall,name="Spawn")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
png(filename="PC_all_rand.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.all.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Spawn")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

## All Initial Samples

``` bash
source activate CASE

mawk -f ../scripts/polarize_freqs all.dp20.rand_rc | grep -Ei '(A{20}|C{20}|G{20}|T{20})' -m 10000 | cut -f43-54  > all.is.dp20.rand.pool
```

``` r
pool.data.all.is <- read.table("all.is.dp20.rand.pool")
df.pool.all.is <- apply(pool.data.all.is, c(1, 2), function(x) eval(parse(text = x)))
df_filtered.all.is <- df.pool.all.is[!apply(df.pool.all.is, 1, function(row) all(row == row[1])), ]

pool.data2.all.is <- t(df_filtered.all.is)
filename.all.is <- read.pcadapt(pool.data2.all.is, type = "pool")

res.all.is.rand <- pcadapt(filename.all.is, K =5,min.maf = 0.0001)

poplist.names <- c("B10","B11","B10","B11","B10","B11","B10","B11","B12","B12","B12","B12")
p1.all.rand.df <- data.frame(res.all.is.rand$scores)
colnames(p1.all.rand.df) <- c("PC1","PC2", "PC3", "PC4")
p1.all.rand.df$POP <- poplist.names
ggplot(p1.all.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="black")+  scale_fill_manual(values=cbPaletteSmall,name="Spawn")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
png(filename="PC_all_IS_rand.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.all.rand.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Spawn")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

### Venn Diagram for called loci

``` bash
cut -f1,2 CASE.dp20.Block10.cov.sync > Block10.pos
cut -f1,2 CASE.dp20.Block11.cov.sync > Block11.pos
cut -f1,2 CASE.dp20.Block12.cov.sync > Block12.pos

cat <(echo "SNP") <(cat Block1*.pos | cut -f1,2 | sort | uniq -c | mawk '$1 <2' | mawk '{print $2"_"$3}') > singleton.loci
```

``` r
# Read and preprocess data
B10.loci <- read.table("Block10.pos", header = TRUE)
colnames(B10.loci) <- c("CHROM", "BP")
B10.loci <- B10.loci %>% mutate(BLOCK_10 = TRUE)

B11.loci <- read.table("Block11.pos", header = TRUE)
colnames(B11.loci) <- c("CHROM", "BP")
B11.loci <- B11.loci %>% mutate(BLOCK_11 = TRUE)

B12.loci <- read.table("Block12.pos", header = TRUE)
colnames(B12.loci) <- c("CHROM", "BP")
B12.loci <- B12.loci %>% mutate(BLOCK_12 = TRUE)

# Full join all data frames to preserve all CHROM and BP combinations
pos_df <- B10.loci %>%
  full_join(B11.loci %>% select(CHROM, BP, BLOCK_11), by = c("CHROM", "BP")) %>%
  full_join(B12.loci %>% select(CHROM, BP, BLOCK_12), by = c("CHROM", "BP"))

# Replace NA values in the BLOCK_* columns with FALSE
pos_df <- pos_df %>%
  mutate(
    BLOCK_10 = ifelse(is.na(BLOCK_10), FALSE, BLOCK_10),
    BLOCK_11 = ifelse(is.na(BLOCK_11), FALSE, BLOCK_11),
    BLOCK_12 = ifelse(is.na(BLOCK_12), FALSE, BLOCK_12)
  )


loci.b10 <- nrow(subset(pos_df,BLOCK_10 == TRUE))
loci.b11 <- nrow(subset(pos_df,BLOCK_11 == TRUE))
loci.b12 <- nrow(subset(pos_df,BLOCK_12 == TRUE))

loci.b10.b11 <- nrow(subset(pos_df,BLOCK_10 == TRUE & BLOCK_11 ==TRUE))
loci.b10.b12 <- nrow(subset(pos_df,BLOCK_10 == TRUE & BLOCK_12 ==TRUE))
loci.b11.b12 <- nrow(subset(pos_df,BLOCK_11 == TRUE & BLOCK_12 ==TRUE))
loci.b10.b11.b12 <- nrow(subset(pos_df,BLOCK_10 == TRUE & BLOCK_11 ==TRUE & BLOCK_12 ==TRUE))

png(filename="VenSNPs.png", type="cairo",units="px", width=4000, height=4000, res=500, bg="transparent")

venn.plot <- draw.triple.venn(
  area1=loci.b10, area2=loci.b11, area3=loci.b12, 
  n12= loci.b10.b11, n13=loci.b10.b12, 
  n23=loci.b11.b12,scaled=TRUE,
  n123=loci.b10.b11.b12, cex=rep(2, 7),
  category = c("Block 10", "Block 11", "Block 12"),
  fill = cbPaletteSmall3,cat.cex = rep(2, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

# Individual Blocks

#### B10

``` bash
source activate CASE

bash ../scripts/sum.sh <(cut -f1,2,3,13-16 CASE.dp20.Block10.cov.sync) > CASE.B10.IS.sync
paste CASE.B10.IS.sync <(cut -f4 CASE.B10.IS.sync) <(cut -f4 CASE.B10.IS.sync)  > CASE.B10.ISsum3.sync

AV_COV=$(mawk -f ../scripts/add_cov_sync <(cut -f13-16 --complement CASE.dp20.Block10.cov.sync) | mawk '{sum=sum+$21} END {print sum/NR}')

echo $AV_COV

mawk -f ../scripts/add_cov_sync CASE.B10.ISsum3.sync | mawk -v x=$AV_COV '$7 >= x' | mawk '!/CHR/' | cut -f1-6 > CASE.B10.ISsum3nh.sync

mawk -f ../scripts/add_cov_sync_IS CASE.B10.ISsum3.sync | mawk -v x=$AV_COV '$7 >= x'     > CASE.B10.ISsum3nh.cov.sync


paste CASE.dp20.Block10.cov.sync <(mawk -f ../scripts/add_cov_sync CASE.B10.ISsum3.sync ) | mawk -v x=$AV_COV '$29 >= x' | cut -f1-22 > temp.sync

mv temp.sync CASE.dp20.Block10.COV.sync
#cp CASE.dp20.Block10.cov.sync CASE.dp20.Block10.COV.sync
mawk -f ../scripts/add_cov_sync <(cut -f13-16,20- --complement CASE.dp20.Block10.COV.sync) > temp.cov


paste temp.cov <(cut -f7 CASE.B10.ISsum3nh.cov.sync) | mawk -v OFS='\t' '{if ($18 > $19) {$18=$19; print $0} else {print $0}}' | cut -f1-18> CASE.dp20.Block10.COV
```

``` bash
source activate random_draw

#python ../scripts/sub_sample.py CASE.dp20.Block10.COV CASE.B10.ISsum3nh.sync CASE.B10.ISs.sync
python ../scripts/sub_sample.py CASE.B10.ISsum3nh.cov.sync CASE.B10.ISsum3nh.sync CASE.B10.ISs.sync

cat <(echo -e "CHROM\tPOS\tREF\tIS_RS1\tIS_RS2\tIS_RS3") CASE.B10.ISs.sync > CASE.B10.ISsum.sync

paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync) <(cut -f4 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f5 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f6 CASE.dp20.Block10.COV.sync)> CA.B10.sync

head -1 CA.B10.sync

mawk '!/CHR/' CA.B10.sync > CA.B10.input

paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync)  <(cut -f7 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f8 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f9 CASE.dp20.Block10.COV.sync)> CASE.B10.sync

head -1 CASE.B10.sync
tail -n +2 CASE.B10.sync > CASE.B10.input


paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync) <(cut -f17 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f19 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f18 CASE.dp20.Block10.COV.sync) > SE.B10.sync

head -1 SE.B10.sync
tail -n +2 SE.B10.sync > SE.B10.input


paste <(cut -f1,2,3,4 CASE.B10.ISsum.sync) <(cut -f10 CASE.dp20.Block10.COV.sync) <(cut -f5 CASE.B10.ISsum.sync) <(cut -f11 CASE.dp20.Block10.COV.sync) <(cut -f6 CASE.B10.ISsum.sync) <(cut -f12 CASE.dp20.Block10.COV.sync)> CON.B10.sync

head -1 CON.B10.sync
tail -n +2 CON.B10.sync | mawk '!/0:0:0:0:0:0/' > CON.B10.input
```

#### B11

``` bash
AV_COV=$(mawk -f ../scripts/add_cov_sync <(cut -f16-19 --complement CASE.dp20.Block11.cov.sync) | mawk '{sum=sum+$22} END {print sum/NR}')

echo $AV_COV

bash ../scripts/sum.sh <(cut -f1,2,3,16-19 CASE.dp20.Block11.cov.sync) > CASE.B11.IS.sync
paste CASE.B11.IS.sync <(cut -f4 CASE.B11.IS.sync) <(cut -f4 CASE.B11.IS.sync) <(cut -f4 CASE.B11.IS.sync) > CASE.B11.ISsum4.sync

mawk -f ../scripts/add_cov_sync CASE.B11.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x' | mawk '!/CHR/' | cut -f1-7 > CASE.B11.ISsum4nh.sync

paste CASE.dp20.Block11.cov.sync <(mawk -f ../scripts/add_cov_sync CASE.B11.ISsum4.sync ) | mawk -v x=$AV_COV '$34 >= x' | cut -f1-23 > temp.sync

mv temp.sync CASE.dp20.Block11.COV.sync

mawk -f ../scripts/add_cov_sync_IS CASE.B11.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x'  > CASE.B11.ISsum4.cov.sync
mawk -f ../scripts/add_cov_sync <(cut -f16-19 --complement CASE.dp20.Block11.COV.sync) > CASE.dp20.Block11.COV

source activate random_draw

#python ../scripts/sub_sample.py CASE.dp20.Block11.COV CASE.B11.ISsum4nh.sync CASE.B11.ISs.sync
python ../scripts/sub_sample.py CASE.B11.ISsum4.cov.sync CASE.B11.ISsum4nh.sync CASE.B11.ISs.sync

cat <(echo -e "CHROM\tPOS\tREF\tIS_RS1\tIS_RS2\tIS_RS3\tIS_R4") CASE.B11.ISs.sync > CASE.B11.ISsum.sync

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f4 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 5 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f6 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f7 CASE.dp20.Block11.COV.sync) > CA.B11.sync

head -1 CA.B11.sync
tail -n +2 CA.B11.sync > CA.B11.input

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f8 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 9 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f10 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f11 CASE.dp20.Block11.COV.sync) > CASE.B11.sync

head -1 CASE.B11.sync
tail -n +2 CASE.B11.sync > CASE.B11.input

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f20 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 21 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f22 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f23 CASE.dp20.Block11.COV.sync) > SE.B11.sync

head -1 SE.B11.sync
tail -n +2 SE.B11.sync > SE.B11.input

paste <(cut -f1,2,3,4 CASE.B11.ISsum.sync) <(cut -f12 CASE.dp20.Block11.COV.sync) <(cut -f5 CASE.B11.ISsum.sync) <(cut -f 13 CASE.dp20.Block11.COV.sync) <(cut -f6 CASE.B11.ISsum.sync) <(cut -f14 CASE.dp20.Block11.COV.sync) <(cut -f7 CASE.B11.ISsum.sync) <(cut -f15 CASE.dp20.Block11.COV.sync) > CON.B11.sync

head -1 CON.B11.sync
tail -n +2 CON.B11.sync > CON.B11.input
```

#### B12

``` bash
source activate CASE

bash ../scripts/sum.sh <(cut -f1,2,3,16-19 CASE.dp20.Block12.cov.sync) > CASE.B12.IS.sync
paste CASE.B12.IS.sync <(cut -f4 CASE.B12.IS.sync) <(cut -f4 CASE.B12.IS.sync) <(cut -f4 CASE.B12.IS.sync) > CASE.B12.ISsum4.sync

AV_COV=$(mawk -f ../scripts/add_cov_sync <(cut -f16-19 --complement CASE.dp20.Block12.cov.sync) | mawk '{sum=sum+$22} END {print sum/NR}')

echo $AV_COV

mawk -f ../scripts/add_cov_sync CASE.B12.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x' | mawk '!/CHR/' | cut -f1-7 > CASE.B12.ISsum4nh.sync

paste CASE.dp20.Block12.cov.sync <(mawk -f ../scripts/add_cov_sync CASE.B12.ISsum4.sync ) | mawk -v x=$AV_COV '$34 >= x' | cut -f1-26 > temp.sync

mv -f temp.sync CASE.dp20.Block12.COV.sync

mawk -f ../scripts/add_cov_sync_IS CASE.B12.ISsum4.sync | mawk -v x=$AV_COV '$8 >= x'  > CASE.B12.ISsum4nh.cov.sync

mawk -f ../scripts/add_cov_sync <(cut -f16-19,24- --complement CASE.dp20.Block12.COV.sync) > CASE.dp20.Block12.COV


paste CASE.dp20.Block12.COV <(cut -f8 CASE.B12.ISsum4nh.cov.sync) | mawk -v OFS='\t' '{if ($22 > $23) {$22=$23; print $0} else {print $0}}' > temp.cov.sync

cut -f1-22 temp.cov.sync > CASE.dp20.Block12.COV
rm temp.cov.sync
```

``` bash
source activate random_draw

#python ../scripts/sub_sample.py CASE.dp20.Block12.COV CASE.B12.ISsum4nh.sync CASE.B12.ISs.sync
python ../scripts/sub_sample.py CASE.B12.ISsum4nh.cov.sync CASE.B12.ISsum4nh.sync CASE.B12.ISs.sync

cat <(echo -e "CHROM\tPOS\tREF\tIS_RS1\tIS_RS2\tIS_RS3\tIS_R4") CASE.B12.ISs.sync > CASE.B12.ISsum.sync

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f4 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 5 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f6 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f7 CASE.dp20.Block12.COV.sync) > CA.B12.sync

head -1 CA.B12.sync
tail -n +2 CA.B12.sync > CA.B12.input

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f8 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 9 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f10 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f11 CASE.dp20.Block12.COV.sync) > CASE.B12.sync

head -1 CASE.B12.sync
tail -n +2 CASE.B12.sync > CASE.B12.input

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f20 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 21 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f22 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f23 CASE.dp20.Block12.COV.sync) > SE.B12.sync

head -1 SE.B12.sync
tail -n +2 SE.B12.sync > SE.B12.input

paste <(cut -f1,2,3,4 CASE.B12.ISsum.sync) <(cut -f12 CASE.dp20.Block12.COV.sync) <(cut -f5 CASE.B12.ISsum.sync) <(cut -f 13 CASE.dp20.Block12.COV.sync) <(cut -f6 CASE.B12.ISsum.sync) <(cut -f14 CASE.dp20.Block12.COV.sync) <(cut -f7 CASE.B12.ISsum.sync) <(cut -f15 CASE.dp20.Block12.COV.sync) > CON.B12.sync

head -1 CON.B12.sync
tail -n +2 CON.B12.sync > CON.B12.input
```

### Calculate P-values

#### Modified CMH

``` r
ca.b10.sync <- read.sync(file="CA.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.b10.cov <- coverage(ca.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
ca.b10.af <- af(ca.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))

case.b10.sync <- read.sync(file="CASE.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.b10.cov <- coverage(case.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
case.b10.af <- af(case.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))

se.b10.sync <- read.sync(file="SE.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.b10.cov <- coverage(se.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
se.b10.af <- af(se.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))

con.b10.sync <- read.sync(file="CON.B10.input", gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.b10.cov <- coverage(con.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
con.b10.af <- af(con.b10.sync,gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3))
```

``` r
ca.b10.pvals <- adapted.cmh.test(freq=ca.b10.af, coverage=ca.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(ca.b10.af)/2), MeanStart = TRUE, mincov =20)

pcadf <- as.data.table(cbind(row.names(ca.b10.af),ca.b10.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"


case.b10.pvals <- adapted.cmh.test(freq=case.b10.af, coverage=case.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(case.b10.af)/2), MeanStart = TRUE, mincov =20)

pcasedf <- as.data.table(cbind(row.names(case.b10.af),case.b10.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"


se.b10.pvals <- adapted.cmh.test(freq=se.b10.af, coverage=se.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(se.b10.af)/2), MeanStart = TRUE, mincov =20)

psedf <- as.data.table(cbind(row.names(se.b10.af),se.b10.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"


con.b10.pvals <- adapted.cmh.test(freq=con.b10.af, coverage=con.b10.cov, Ne=rep(250, 3), gen=c(0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3,3), poolSize=rep(c(40000,10000), ncol(con.b10.af)/2), MeanStart = TRUE, mincov =20)

pcondf <- as.data.table(cbind(row.names(con.b10.af),con.b10.pvals))
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"

B10.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
B10.pval.table <- B10.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
B10.pval.table$PVAL <- as.numeric(B10.pval.table$PVAL)
```

``` r
ca.b11.sync <- read.sync(file="CA.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.b11.cov <- coverage(ca.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
ca.b11.af <- af(ca.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

case.b11.sync <- read.sync(file="CASE.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.b11.cov <- coverage(case.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
case.b11.af <- af(case.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

se.b11.sync <- read.sync(file="SE.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.b11.cov <- coverage(se.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
se.b11.af <- af(se.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

con.b11.sync <- read.sync(file="CON.B11.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.b11.cov <- coverage(con.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
con.b11.af <- af(con.b11.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

``` r
ca.b11.pvals <- adapted.cmh.test(freq=ca.b11.af, coverage=ca.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(ca.b11.af)/2), MeanStart = TRUE, mincov =14)

pcadf <- as.data.table(cbind(row.names(ca.b11.af),ca.b11.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"

case.b11.pvals <- adapted.cmh.test(freq=case.b11.af, coverage=case.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(case.b11.af)/2), MeanStart = TRUE, mincov =14)

pcasedf <- as.data.table(cbind(row.names(case.b11.af),case.b11.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"

se.b11.pvals <- adapted.cmh.test(freq=se.b11.af, coverage=se.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(se.b11.af)/2), MeanStart = TRUE, mincov =14)

psedf <- as.data.table(cbind(row.names(se.b11.af),se.b11.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"

con.b11.pvals <- adapted.cmh.test(freq=con.b11.af, coverage=con.b11.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(con.b11.af)/2), MeanStart = TRUE, mincov =14)

pcondf <- as.data.table(cbind(row.names(con.b11.af),con.b11.pvals))
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"


B11.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
B11.pval.table <- B11.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
B11.pval.table$PVAL <- as.numeric(B11.pval.table$PVAL)
```

``` r
ca.b12.sync <- read.sync(file="CA.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.b12.cov <- coverage(ca.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
ca.b12.af <- af(ca.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

case.b12.sync <- read.sync(file="CASE.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.b12.cov <- coverage(case.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
case.b12.af <- af(case.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

se.b12.sync <- read.sync(file="SE.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.b12.cov <- coverage(se.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
se.b12.af <- af(se.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))

con.b12.sync <- read.sync(file="CON.B12.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.b12.cov <- coverage(con.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
con.b12.af <- af(con.b12.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4 ))
```

``` r
ca.b12.pvals <- adapted.cmh.test(freq=ca.b12.af, coverage=ca.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(ca.b12.af)/2), MeanStart = TRUE, mincov =20)

pcadf <- as.data.table(cbind(row.names(ca.b12.af),ca.b12.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"

case.b12.pvals <- adapted.cmh.test(freq=case.b12.af, coverage=case.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(case.b12.af)/2), MeanStart = TRUE, mincov =20)

pcasedf <- as.data.table(cbind(row.names(case.b12.af),case.b12.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"

se.b12.pvals <- adapted.cmh.test(freq=se.b12.af, coverage=se.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(se.b12.af)/2), MeanStart = TRUE, mincov =20)

psedf <- as.data.table(cbind(row.names(se.b12.af),se.b12.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"



con.b12.pvals <- adapted.cmh.test(freq=con.b12.af, coverage=con.b12.cov, Ne=rep(1000, 4), gen=c(0, 1, 0, 1, 0, 1, 0, 1), repl=c(1, 1, 2, 2, 3, 3, 4, 4), poolSize=rep(c(40000,10000), ncol(con.b12.af)/2), MeanStart = TRUE, mincov =20)

pcondf <- as.data.table(cbind(row.names(con.b12.af),con.b12.pvals))
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"

B12.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
B12.pval.table <- B12.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
B12.pval.table$PVAL <- as.numeric(B12.pval.table$PVAL)
```

``` r
Qvalue_convert(B10.pval.table) ->B10.pv
Qvalue_convert(B11.pval.table) ->B11.pv
Qvalue_convert(B12.pval.table) ->B12.pv
```

``` r
pp1 <- rbind(B10.pv,B11.pv,B12.pv)

alpha = 0.1
B11.alpha = alpha * 1
alpha2 =0.1

ppsig.B10.1 <- Significat_subset(B10.pv,alpha,alpha2)
```

    ## [1] 495239
    ## [1] 85765
    ## [1] 0.173179

``` r
ppsig.B11.1 <- Significat_subset(B11.pv,B11.alpha,alpha2)
```

    ## [1] 175041
    ## [1] 14925
    ## [1] 0.08526574

``` r
ppsig.B12.1 <- Significat_subset(B12.pv,alpha,alpha2)
```

    ## [1] 580904
    ## [1] 77241
    ## [1] 0.1329669

``` r
ppsig.B10.1$BLOCK <- 10
ppsig.B11.1$BLOCK <- 11
ppsig.B12.1$BLOCK <- 12

ppsig3.FDR10 <- bind_rows(ppsig.B10.1,ppsig.B11.1,ppsig.B12.1)

all_sig.CA<- subset(ppsig3.FDR10, Sig.CA == "TRUE") %>% group_by(SNP) %>% filter(n()>2)
all_sig.CA$Sig.CASE<- FALSE
all_sig.CA$Sig.SE<- FALSE

all_sig.CASE<- subset(ppsig3.FDR10, Sig.CASE == "TRUE") %>% group_by(SNP) %>% filter(n()>2)
all_sig.CASE$Sig.CA<- FALSE
all_sig.CASE$Sig.SE<- FALSE

all_sig.SE<- subset(ppsig3.FDR10, Sig.SE == "TRUE") %>% group_by(SNP) %>% filter(n()>2)
all_sig.SE$Sig.CASE<- FALSE
all_sig.SE$Sig.CA<- FALSE

#all_sig <- ppsig3.FDR10 %>% group_by(SNP) %>% filter(n()>2)
all_sig <- bind_rows(all_sig.CA,all_sig.CASE,all_sig.SE)

write.table(all_sig, "Sig.Loci.3", sep="\t", row.names = FALSE, quote = FALSE)

alpha = 0.05
B11.alpha = alpha * 2
alpha2 =0.1

ppsig.B10.05 <- Significat_subset(B10.pv,alpha,alpha2)
```

    ## [1] 495239
    ## [1] 59545
    ## [1] 0.1202349

``` r
ppsig.B11.05 <- Significat_subset(B11.pv,B11.alpha,alpha2)
```

    ## [1] 175041
    ## [1] 14925
    ## [1] 0.08526574

``` r
ppsig.B12.05 <- Significat_subset(B12.pv,alpha,alpha2)
```

    ## [1] 580904
    ## [1] 42661
    ## [1] 0.07343898

``` r
ppsig.B10.05$BLOCK <- 10
ppsig.B11.05$BLOCK <- 11
ppsig.B12.05$BLOCK <- 12

ppsig3.FDR05 <- rbind(ppsig.B10.05,ppsig.B11.05,ppsig.B12.05)

multi_sig.CA<- subset(ppsig3.FDR05, Sig.CA == "TRUE") %>% group_by(SNP) %>% filter(n()>1)
multi_sig.CA$Sig.CASE <- FALSE
multi_sig.CA$Sig.SE <- FALSE

multi_sig.CASE<- subset(ppsig3.FDR05, Sig.CASE == "TRUE") %>% group_by(SNP) %>% filter(n()>1)
multi_sig.CASE$Sig.CA <- FALSE
multi_sig.CASE$Sig.SE <- FALSE

multi_sig.SE<- subset(ppsig3.FDR05, Sig.SE == "TRUE") %>% group_by(SNP) %>% filter(n()>1)
multi_sig.SE$Sig.CA <- FALSE
multi_sig.SE$Sig.CASE <- FALSE

multi_sig <- rbind(multi_sig.CA,multi_sig.CASE,multi_sig.SE)

#multi_sig.CA<- subset(ppsig3.FDR05, QCA < B11.alpha) %>% group_by(SNP) %>% filter(n()>1)
#multi_sig.CASE<- subset(ppsig3.FDR05, QCASE < B11.alpha) %>% group_by(SNP) %>% filter(n()>1)
#multi_sig.SE<- subset(ppsig3.FDR05, QSE < B11.alpha) %>% group_by(SNP) %>% filter(n()>1)

#multi_sig <- rbind(ppsig3.FDR05 %>% group_by(SNP) %>% filter(n()>1))
#multi_sig <- rbind(multi_sig.CA,multi_sig.CASE,multi_sig.SE)

write.table(multi_sig, "Sig.Loci.2", sep="\t", row.names = FALSE, quote = FALSE)

alpha = 0.01
alpha2 =0.1
B11.alpha = alpha *2

ppsig.B10.01 <- Significat_subset(B10.pv,alpha,alpha2)
```

    ## [1] 495239
    ## [1] 30021
    ## [1] 0.06061922

``` r
ppsig.B11.01 <- Significat_subset(B11.pv,B11.alpha,alpha2)
```

    ## [1] 175041
    ## [1] 3137
    ## [1] 0.01792152

``` r
ppsig.B12.01 <- Significat_subset(B12.pv,alpha,alpha2)
```

    ## [1] 580904
    ## [1] 12544
    ## [1] 0.02159393

``` r
conpp1 <- subset(pp1, QCON < alpha)

ppsig.B10.01$BLOCK <- 10
ppsig.B11.01$BLOCK <- 11
ppsig.B12.01$BLOCK <- 12

ppsig3.FDR01 <- rbind(ppsig.B10.01,ppsig.B11.01,ppsig.B12.01)

write.table(ppsig3.FDR01, "Sig.Loci.FDR01.1", sep="\t", row.names = FALSE, quote = FALSE)
write.table(ppsig.B10.01, "Sig.Loci.FDR01.Block10", sep="\t", row.names = FALSE, quote = FALSE)
write.table(ppsig.B11.01, "Sig.Loci.FDR01.Block11", sep="\t", row.names = FALSE, quote = FALSE)
write.table(ppsig.B12.01, "Sig.Loci.FDR01.Block12", sep="\t", row.names = FALSE, quote = FALSE)

multi_sig.CA.01<- subset(ppsig3.FDR01, Sig.CA == "TRUE") %>% group_by(SNP) %>% filter(n()>1)
multi_sig.CA.01$Sig.CASE <- FALSE
multi_sig.CA.01$Sig.SE <- FALSE

multi_sig.CASE.01<- subset(ppsig3.FDR01, Sig.CASE == "TRUE") %>% group_by(SNP) %>% filter(n()>1)
multi_sig.CASE.01$Sig.CA <- FALSE
multi_sig.CASE.01$Sig.SE <- FALSE

multi_sig.SE.01<- subset(ppsig3.FDR01, Sig.SE == "TRUE") %>% group_by(SNP) %>% filter(n()>1)
multi_sig.SE.01$Sig.CA <- FALSE
multi_sig.SE.01$Sig.CASE <- FALSE

multi_sig.01 <- rbind(multi_sig.CA.01,multi_sig.CASE.01,multi_sig.SE.01)

write.table(multi_sig.01, "Sig.Loci.FDR01.2", sep="\t", row.names = FALSE, quote = FALSE)

alpha = 0.01
alpha2 =0.1
B11.alpha = alpha * 2


ppsig.B10 <- Significat_subset(B10.pv,alpha,alpha2)
```

    ## [1] 495239
    ## [1] 30021
    ## [1] 0.06061922

``` r
ppsig.B10 <- subset(ppsig.B10, Sig.CASE == TRUE & QCASE < quantile(ppsig.B10$QCASE, na.rm = TRUE, probs = 0.01)
                    | Sig.CA == TRUE & QCA < quantile(ppsig.B10$QCA, na.rm = TRUE, probs = 0.01)
                    | Sig.SE == TRUE & QSE < quantile(ppsig.B10$QSE, na.rm = TRUE, probs = 0.01)
                    )
ppsig.B11 <- Significat_subset(B11.pv,B11.alpha,alpha2)
```

    ## [1] 175041
    ## [1] 3137
    ## [1] 0.01792152

``` r
ppsig.B11 <- subset(ppsig.B11, Sig.CASE == TRUE & QCASE < quantile(ppsig.B11$QCASE, na.rm = TRUE, probs = 0.01)
                    | Sig.CA == TRUE & QCA < quantile(ppsig.B11$QCA, na.rm = TRUE, probs = 0.01)
                    | Sig.SE == TRUE & QSE < quantile(ppsig.B11$QSE, na.rm = TRUE, probs = 0.01)
                    )
ppsig.B12 <- Significat_subset(B12.pv,alpha,alpha2)
```

    ## [1] 580904
    ## [1] 12544
    ## [1] 0.02159393

``` r
ppsig.B12 <- subset(ppsig.B12, Sig.CASE == TRUE & QCASE < quantile(ppsig.B12$QCASE, na.rm = TRUE, probs = 0.01)
                    | Sig.CA == TRUE & QCA < quantile(ppsig.B12$QCA, na.rm = TRUE, probs = 0.01)
                    | Sig.SE == TRUE & QSE < quantile(ppsig.B12$QSE, na.rm = TRUE, probs = 0.01)
                    )

ppsig.B10$BLOCK <- 10
ppsig.B11$BLOCK <- 11
ppsig.B12$BLOCK <- 12

singleton <- read.table("singleton.loci",header=TRUE)

ppsig1.s <- semi_join(rbind(ppsig.B10,ppsig.B11,ppsig.B12),singleton, by = "SNP")
```

# Across all blocks

#### Create sync file

``` bash
source activate CASE

bash ../scripts/sum.sh <(cut -f1,2,3,42,44,46,48 CASE.All.Blocks.cov.sync) > CASE.B10.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,43,45,47,49 CASE.All.Blocks.cov.sync) > CASE.B11.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,50-53 CASE.All.Blocks.cov.sync) > CASE.B12.AB.ISsum.sync


paste CASE.All.Blocks.cov.sync <(mawk -f ../scripts/add_cov_sync CASE.B10.AB.ISsum.sync | cut -f5 ) <(mawk -f ../scripts/add_cov_sync CASE.B11.AB.ISsum.sync | cut -f 5) <(mawk -f ../scripts/add_cov_sync CASE.B12.AB.ISsum.sync | cut -f5) | mawk '$69 >72 && $70 > 72 && $71 > 72 && $66 > 1' | cut -f1-68 > CASE.All.Blocks.cov.fil.sync

bash ../scripts/sum.sh <(cut -f1,2,3,42,44,46 CASE.All.Blocks.cov.fil.sync) > CASE.B10.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,43,45,47,49 CASE.All.Blocks.cov.fil.sync) > CASE.B11.AB.ISsum.sync
bash ../scripts/sum.sh <(cut -f1,2,3,50-53 CASE.All.Blocks.cov.fil.sync) > CASE.B12.AB.ISsum.sync

mawk -f ../scripts/add_cov_sync <(cut -f1-3,4,8,10,18,22,25,30,34,36,56,57,60  CASE.All.Blocks.cov.fil.sync ) > B10.AB.cov.sync
mawk -f ../scripts/add_cov_sync <(cut -f1-3,6,11,16,19,26,27,31,33,39,54,59,63 CASE.All.Blocks.cov.fil.sync) > B11.AB.cov.sync
mawk -f ../scripts/add_cov_sync <(cut -f1-3,5,7,14,20,21,23,32,35,38,61,62,64 CASE.All.Blocks.cov.fil.sync) > B12.AB.cov.sync

#mawk -f ../scripts/add_cov_sync <(cut -f1-3,4,8,10,18,22,25,30,34,36,42,44,46,48,56,57,60  CASE.All.Blocks.cov.fil.sync ) > B10.AB.cov.sync
#mawk -f ../scripts/add_cov_sync <(cut -f1-3,6,11,16,19,26,27,31,33,39,43,45,47,49,54,59,63 CASE.All.Blocks.cov.fil.sync) > B11.AB.cov.sync
#mawk -f ../scripts/add_cov_sync <(cut -f1-3,5,7,14,20,21,23,32,35,38,50,51,52,53,61,62,64 CASE.All.Blocks.cov.fil.sync) > B12.AB.cov.sync

paste CASE.B10.AB.ISsum.sync <(cut -f4 CASE.B10.AB.ISsum.sync) <(cut -f4 CASE.B10.AB.ISsum.sync)  > CASE.B10.AB.ISsum3.sync

paste CASE.B11.AB.ISsum.sync <(cut -f4 CASE.B11.AB.ISsum.sync) <(cut -f4 CASE.B11.AB.ISsum.sync)  > CASE.B11.AB.ISsum4.sync

paste CASE.B12.AB.ISsum.sync <(cut -f4 CASE.B12.AB.ISsum.sync) <(cut -f4 CASE.B12.AB.ISsum.sync)  > CASE.B12.AB.ISsum4.sync

mawk -f ../scripts/add_cov_sync_IS CASE.B10.AB.ISsum3.sync  > CASE.B10.AB.ISsum3.cov.sync
mawk -f ../scripts/add_cov_sync_IS CASE.B11.AB.ISsum4.sync  > CASE.B11.AB.ISsum4.cov.sync
mawk -f ../scripts/add_cov_sync_IS CASE.B12.AB.ISsum4.sync  > CASE.B12.AB.ISsum4.cov.sync

paste B10.AB.cov.sync <(cut -f7 CASE.B10.AB.ISsum3.cov.sync) | mawk -v OFS='\t' '{if (NR > 1 && $18 > $19) {$18=$19; print $0} else {print $0}}' > temp.cov.sync
cut -f1-18 temp.cov.sync > B10.AB.cov.sync

paste B11.AB.cov.sync <(cut -f7 CASE.B11.AB.ISsum4.cov.sync) | mawk -v OFS='\t' '{if (NR > 1 && $18 > $19) {$18=$19; print $0} else {print $0}}' > temp.cov.sync
cut -f1-18 temp.cov.sync > B11.AB.cov.sync

paste B12.AB.cov.sync <(cut -f7 CASE.B12.AB.ISsum4.cov.sync) | mawk -v OFS='\t' '{if (NR > 1 && $18 > $19) {$18=$19; print $0} else {print $0}}' > temp.cov.sync
cut -f1-18 temp.cov.sync > B12.AB.cov.sync
rm temp.cov.sync
```

``` bash
source activate random_draw
python ../scripts/sub_sample.py CASE.B10.AB.ISsum3.cov.sync <( mawk '!/CHR/' CASE.B10.AB.ISsum3.sync) CASE.B10.ISs.AB.sync
python ../scripts/sub_sample.py CASE.B11.AB.ISsum4.cov.sync <( mawk '!/CHR/' CASE.B11.AB.ISsum4.sync) CASE.B11.ISs.AB.sync
python ../scripts/sub_sample.py CASE.B12.AB.ISsum4.cov.sync <(mawk '!/CHR/' CASE.B12.AB.ISsum4.sync) CASE.B12.ISs.AB.sync


#python ../scripts/sub_sample.py B10.AB.cov.sync <( mawk '!/CHR/' CASE.B10.AB.ISsum3.sync) CASE.B10.ISs.AB.sync
#python ../scripts/sub_sample.py B11.AB.cov.sync <( mawk '!/CHR/' CASE.B11.AB.ISsum4.sync) CASE.B11.ISs.AB.sync
#python ../scripts/sub_sample.py B12.AB.cov.sync <(mawk '!/CHR/' CASE.B12.AB.ISsum4.sync) CASE.B12.ISs.AB.sync


wait 

cat <(echo -e "CHROM\tPOS\tREF\tISB11_RS1\tISB11_RS2\tISB11_RS3") CASE.B11.ISs.AB.sync > CASE.B11.ISsum.AB.sync
cat <(echo -e "CHROM\tPOS\tREF\tISB10_RS1\tISB10_RS2\tISB10_RS3") CASE.B10.ISs.AB.sync > CASE.B10.ISsum.AB.sync
cat <(echo -e "CHROM\tPOS\tREF\tISB12_RS1\tISB12_RS2\tISB12_RS3") CASE.B12.ISs.AB.sync > CASE.B12.ISsum.AB.sync

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f18 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f22 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f25 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f19 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f26 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f27 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f20 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f21 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f23 CASE.All.Blocks.cov.fil.sync) > CASE.AB.sync

tail -n +2 CASE.AB.sync | mawk '$2 != 74609903 && $2 != 93711267 && $2 != 36030885' > CASE.AB.input 

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f4 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f8 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f10 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f6 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f11 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f16 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f5 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f7 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f14 CASE.All.Blocks.cov.fil.sync) > CA.AB.sync

tail -n +2 CA.AB.sync | mawk '$2 != 93711267 && $2 != 59024150 && $2 != 2532550'> CA.AB.input 

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f30 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f34 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f36 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f31 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f33 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f39 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f32 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f35 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f38 CASE.All.Blocks.cov.fil.sync) > CON.AB.sync

tail -n +2 CON.AB.sync | mawk '$2 != 34804390 && $2 != 54895214 && $2 !=  49000140 && $2 != 15140923' > CON.AB.input 

paste <(cut -f1-4 CASE.B10.ISsum.AB.sync) <(cut -f56 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B10.ISsum.AB.sync) <(cut -f57 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B10.ISsum.AB.sync) <(cut -f60 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B11.ISsum.AB.sync) <(cut -f54 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B11.ISsum.AB.sync) <(cut -f59 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B11.ISsum.AB.sync) <(cut -f63 CASE.All.Blocks.cov.fil.sync) <(cut -f4 CASE.B12.ISsum.AB.sync) <(cut -f61 CASE.All.Blocks.cov.fil.sync) <(cut -f5 CASE.B12.ISsum.AB.sync) <(cut -f62 CASE.All.Blocks.cov.fil.sync) <(cut -f6 CASE.B12.ISsum.AB.sync) <(cut -f64 CASE.All.Blocks.cov.fil.sync) > SE.AB.sync

tail -n +2 SE.AB.sync | mawk '$2 != 67902213 && $2 != 19828063 && $2 != 93711267 ' > SE.AB.input 
```

#### Calculate p-values

``` r
ca.ab.sync <- read.sync(file="CA.AB.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
ca.ab.cov <- coverage(ca.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
ca.ab.af <- af(ca.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))

case.ab.sync <- read.sync(file="CASE.AB.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
case.ab.cov <- coverage(case.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
case.ab.af <- af(case.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))

se.ab.sync <- read.sync(file="SE.AB.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
se.ab.cov <- coverage(se.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
se.ab.af <- af(se.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))

con.ab.sync <- read.sync(file="CON.AB.input", gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
```

    ## Reading sync file ...
    ## Extracting biallelic counts ...
    ## Creating result object ...

``` r
con.ab.cov <- coverage(con.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
con.ab.af <- af(con.ab.sync,gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ))
```

``` r
ca.ab.pvals <- adapted.cmh.test(freq=ca.ab.af, coverage=ca.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000), gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ), poolSize=rep(c(40000,10000), ncol(ca.ab.af)/2), MeanStart = TRUE, mincov =14)

pcadf <- as.data.table(cbind(row.names(ca.ab.af),ca.ab.pvals))
colnames(pcadf) <- c("Locus","PVAL")
pcadf$GROUP <- "PCA"

case.ab.pvals <- adapted.cmh.test(freq=case.ab.af, coverage=case.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000), gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ), poolSize=rep(c(40000,10000), ncol(case.ab.af)/2), MeanStart = TRUE, mincov =14)

pcasedf <- as.data.table(cbind(row.names(case.ab.af),case.ab.pvals))
colnames(pcasedf) <- c("Locus","PVAL")
pcasedf$GROUP <- "PCASE"

se.ab.pvals <- adapted.cmh.test(freq=se.ab.af, coverage=se.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000), gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ), poolSize=rep(c(40000,10000), ncol(se.ab.af)/2), MeanStart = TRUE, mincov =14)

psedf <- as.data.table(cbind(row.names(se.ab.af),se.ab.pvals))
colnames(psedf) <- c("Locus","PVAL")
psedf$GROUP <- "PSE"

con.ab.pvals <- adapted.cmh.test(freq=con.ab.af, coverage=con.ab.cov, Ne=c(250,250,250,1000,1000,1000,1000,1000,1000), gen=c(0, 1, 0, 1, 0, 1, 0, 1,0,1,0,1,0,1,0,1,0,1), repl=c(1, 1, 2, 2, 3, 3, 4 , 4,5,5,6,6,7,7,8,8,9,9 ), poolSize=rep(c(40000,10000), ncol(con.ab.af)/2), MeanStart = TRUE, mincov =14)

pcondf <- as.data.table(cbind(row.names(con.ab.af),con.ab.pvals)) 
colnames(pcondf) <- c("Locus","PVAL")
pcondf$GROUP <- "PCON"


AB.pval.table <- bind_rows(pcadf,pcasedf,pcondf,psedf)
AB.pval.table <- AB.pval.table %>%
  separate(Locus, into = c("CHROM", "BP"), sep = "\\.(?=[^.]+$)")
AB.pval.table$PVAL <- as.numeric(AB.pval.table$PVAL)
```

``` r
Qvalue_convert(AB.pval.table) -> AB.pv
```

``` r
alpha = 0.01
alpha2 =0.1

ppsig.AB <- Significat_subset(AB.pv,alpha,alpha2)
```

    ## [1] 393734
    ## [1] 15320
    ## [1] 0.03890952

``` r
all_sig.AB.CA<- subset(ppsig.AB, Sig.CA == "TRUE") 
all_sig.AB.CA$Sig.CASE <- FALSE
all_sig.AB.CA$Sig.SE <- FALSE

all_sig.AB.CASE<- subset(ppsig.AB, Sig.CASE == "TRUE") 
all_sig.AB.CASE$Sig.SE <- FALSE
all_sig.AB.CASE$Sig.CA <- FALSE

all_sig.AB.SE<- subset(ppsig.AB, Sig.SE == "TRUE") 
all_sig.AB.SE$Sig.CASE <- FALSE
all_sig.AB.SE$Sig.CA <- FALSE

all_sig.AB.AT <- rbind(all_sig.AB.CA,all_sig.AB.CASE,all_sig.AB.SE)

write.table(all_sig.AB.AT, "Sig.Loci.FDR01.AB", sep="\t", row.names = FALSE, quote = FALSE)


alphaAB = 0.0001
alpha2 =0.1


ppsig1.ab <-Significat_subset(AB.pv,alphaAB,alpha2)
```

    ## [1] 393734
    ## [1] 1813
    ## [1] 0.004604632

``` r
write.table(ppsig1.ab, "Sig.Loci.FDR01.AB", sep="\t", row.names = FALSE, quote = FALSE)

conpp.AB <- subset(AB.pv, QCON < 0.001)
all_con <- bind_rows(conpp1,conpp.AB) %>% group_by(SNP) %>% filter(n()>1)
#all_con <- bind_rows(conpp1,conpp.AB) 
all_con <- subset(all_con,QCON / QCASE < 10 & QCON / QCA < 10 & QCON/ QSE < 10)

write.table(all_con, "ALL.CON.SNPS", sep="\t", row.names = FALSE, quote = FALSE)

ppsig.AB$BLOCK <- 33
ppsig1.ab$BLOCK <-33

#multiple.sig.o<- bind_rows(ppsig1 %>% distinct(SNP, .keep_all = TRUE),ppsig1.ab) %>% group_by(SNP) %>% filter(n()>1)

multiple.sig.SE <- bind_rows(subset(ppsig3.FDR01, Sig.SE == TRUE), subset(ppsig1.ab, Sig.SE == TRUE)) %>% group_by(SNP) %>% filter(n()>1)
multiple.sig.SE$Sig.CA <- FALSE
multiple.sig.SE$Sig.CASE <- FALSE

multiple.sig.CASE <- bind_rows(subset(ppsig3.FDR01, Sig.CASE == TRUE), subset(ppsig1.ab, Sig.CASE == TRUE)) %>% group_by(SNP) %>% filter(n()>1)
multiple.sig.CASE$Sig.CA <- FALSE
multiple.sig.CASE$Sig.SE <- FALSE

multiple.sig.CA<- bind_rows(subset(ppsig3.FDR01, Sig.CA == TRUE), subset(ppsig1.ab, Sig.CA == TRUE)) %>% group_by(SNP) %>% filter(n()>1)
multiple.sig.CA$Sig.SE <- FALSE
multiple.sig.CA$Sig.CASE <- FALSE

multiple.sig.o <- bind_rows(multiple.sig.SE, multiple.sig.CASE, multiple.sig.CA)

total.sig <- anti_join(bind_rows(multiple.sig.o,all_sig,multi_sig, multi_sig.01, ppsig1.s),all_con, by = "SNP")

write.table(total.sig, "Total.Significant.uf.Loci", sep="\t", row.names = FALSE, quote = FALSE)
```

``` bash
source activate CASE

cut -f2,3 Total.Significant.uf.Loci | tail -n +2 | sort -k1,2 | uniq | mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > CASE.Total.Significant.uf.loci.bed

bedtools intersect -wb -a  <(mawk '{print $2"\t"$3-1"\t"$3}' ALL.CON.SNPS ) -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > CON.LOCI

bedtools intersect -wb -a  <(mawk '{print $2"\t"$3-1"\t"$3}' ALL.CON.SNPS ) -b ~/CASE/analysis/sorted.ref3.0.gene.bed | cut -f4- | bedtools sort | bedtools merge > con.loci.filter.bed

cat <(echo -e "SNP\tCHROM\tBP") <(bedtools intersect -a CASE.Total.Significant.uf.loci.bed -b con.loci.filter.bed | mawk '{print $1"_"$3"\t"$1"\t"$3}') > ToFilter.Con.Loci.txt
```

``` r
con_filter <- read.table("ToFilter.Con.Loci.txt",header=TRUE)

total.sig <- unique(anti_join(total.sig,con_filter, by = "SNP"))

write.table(total.sig, "Total.Significant.Loci", sep="\t", row.names = FALSE, quote = FALSE)


grouped_total.sig <- group_and_average(total.sig)

write.table(grouped_total.sig, "Total.Significant.Grouped.Loci", sep="\t", row.names = FALSE, quote = FALSE)

totalsig.CA <- subset(total.sig, Sig.CA == TRUE & Sig.SE == FALSE )
totalsig.SE <- subset(total.sig, Sig.CA == FALSE & Sig.SE == TRUE )
totalsig.both <- subset(total.sig, Sig.CA == TRUE & Sig.SE == TRUE )
totalsig.CASE <- subset(total.sig, Sig.CASE == TRUE )

write.table(totalsig.CA, "sig.ca", sep="\t", row.names = FALSE, quote = FALSE)
write.table(totalsig.SE, "sig.se", sep="\t", row.names = FALSE, quote = FALSE)
write.table(totalsig.both, "sig.both", sep="\t", row.names = FALSE, quote = FALSE)
write.table(totalsig.CASE, "sig.case", sep="\t", row.names = FALSE, quote = FALSE)
```

``` bash
cat <(echo -e "SNP\tGroup") <(mawk '{if (NR > 1) print $1"\tSE"}' sig.se ) <(mawk '{if (NR > 1) print $1"\tCA"}' sig.ca ) <(mawk '{if (NR > 1) print $1"\tBoth"}' sig.both ) > sig.table
```

``` bash
source activate CASE

cut -f2,3 Total.Significant.Loci | tail -n +2 | sort -k1,2 | uniq | mawk '{print $1 "\t" $2-1 "\t" $2}' | bedtools sort -i - > CASE.Total.Significant.loci.bed
bedtools slop -b 100 -i CASE.Total.Significant.loci.bed  -g ../reference.fasta.fai  | bedtools merge -i - >  CASE.Total.Significant.loci.intervals.bed


bedtools intersect -wb -a CASE.Total.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.TOTAL.CASE.LOC
```

``` bash
source activate CASE

mawk '$14 == "TRUE" ' Total.Significant.Loci | sort -k1,2 | mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CA.Significant.loci.bed
bedtools intersect -wb -a CA.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.1.CA.LOC

mawk '$15 =="TRUE"' Total.Significant.Loci | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}'| uniq  > SE.Significant.loci.bed
bedtools intersect -wb -a SE.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.1.SE.LOC

mawk '$13 == "TRUE"' Total.Significant.Loci | sort -k1,2 |   mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.loci.bed
bedtools intersect -wb -a CASE.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.1.CASE.LOC

mawk '$11 == "TRUE" && $13 == "TRUE"  && $12 == "FALSE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.CASE.Significant.loci.bed

bedtools intersect -wb -a SE.CASE.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.SE.CASE.LOC

mawk '$11 == "TRUE" && $13 == "TRUE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.CASE.CA.Significant.loci.bed

bedtools intersect -wb -a SE.CASE.CA.Significant.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.SE.CASE.CA.LOC

mawk '$11 == "FALSE" && $13 == "TRUE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.CA.Significant.loci.bed

mawk '$11 == "TRUE" && $13 == "FALSE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > CA.CASE.Significant.loci.bed

mawk '$11 == "FALSE" && $13 == "TRUE"  && $12 == "FALSE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > SE.ONLY.Significant.loci.bed

mawk '$11 == "TRUE" && $13 == "FALSE"  && $12 == "FALSE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > CASE.ONLY.Significant.loci.bed

mawk '$11 == "FALSE" && $13 == "FALSE"  && $12 == "TRUE"' Total.Significant.Grouped.Loci | mawk '{print $14 "\t" $15-1 "\t" $15}'| uniq  > CA.ONLY.Significant.loci.bed

cat CA.ONLY.Significant.loci.bed CA.CASE.Significant.loci.bed| sort|  uniq  > CA.factor.Significant.loci.bed

cat SE.ONLY.Significant.loci.bed SE.CASE.Significant.loci.bed| sort | uniq  > SE.factor.Significant.loci.bed
```

``` bash
source activate CASE


#B10
mawk '$14 == "TRUE"' Sig.Loci.FDR01.Block10 | sort -k1,2  | mawk '{print $2 "\t" $3-1 "\t" $3}'| uniq  > CA.Significant.B10.loci.bed
bedtools intersect -wb -a CA.Significant.B10.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B10.CA.LOC

mawk '$15 == "TRUE" ' Sig.Loci.FDR01.Block10 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > SE.Significant.B10.loci.bed
bedtools intersect -wb -a SE.Significant.B10.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B10.SE.LOC

mawk '$13 == "TRUE" ' Sig.Loci.FDR01.Block10 | sort -k1,2 |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.B10.loci.bed
bedtools intersect -wb -a CASE.Significant.B10.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B10.CASE.LOC

#B11
mawk '$14 == "TRUE"' Sig.Loci.FDR01.Block11 | sort -k1,2  | mawk '{print $2 "\t" $3-1 "\t" $3}'  | uniq > CA.Significant.B11.loci.bed
bedtools intersect -wb -a CA.Significant.B11.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B11.CA.LOC

mawk '$15 == "TRUE" ' Sig.Loci.FDR01.Block11 | sort -k1,2 |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq > SE.Significant.B11.loci.bed
bedtools intersect -wb -a SE.Significant.B11.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B11.SE.LOC

mawk '$13 == "TRUE" ' Sig.Loci.FDR01.Block11 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.B11.loci.bed
bedtools intersect -wb -a CASE.Significant.B11.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B11.CASE.LOC

#B12
mawk '$14 == "TRUE" ' Sig.Loci.FDR01.Block12 | sort -k1,2 | mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CA.Significant.B12.loci.bed
bedtools intersect -wb -a CA.Significant.B12.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B12.CA.LOC

mawk '$15 == "TRUE" ' Sig.Loci.FDR01.Block12 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > SE.Significant.B12.loci.bed
bedtools intersect -wb -a SE.Significant.B12.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B12.SE.LOC

mawk '$13 == "TRUE" ' Sig.Loci.FDR01.Block12 | sort -k1,2  |  mawk '{print $2 "\t" $3-1 "\t" $3}' | uniq  > CASE.Significant.B12.loci.bed
bedtools intersect -wb -a CASE.Significant.B12.loci.bed -b ~/CASE/analysis/sorted.ref3.0.gene.bed | grep -oh "gene=LOC.*;g" | sed 's/gene=//g' | sed 's/;g//g' | sort | uniq > Sig.loci.B12.CASE.LOC
```

``` r
source("../gprofiler2W.R")
set_base_url("http://biit.cs.ut.ee/gprofiler_beta")

gostres <- gost(query = list("CASE" = read.table("Sig.loci.1.CASE.LOC",header=FALSE), "SE" = read.table("Sig.loci.1.SE.LOC",header=FALSE), "CA" = read.table("Sig.loci.1.CA.LOC",header=FALSE)), organism = "cvgca002022765v4", multi_query = FALSE,user_threshold = 0.05,domain_scope = "known")

p <-gostplot(gostres, capped = TRUE, interactive = FALSE)

p <- p+ theme_black() + theme(legend.position = "none")

pp <- publish_gostplot(p, width = NA, height = NA, filename = NULL )

png(filename="TestGO.png", type="cairo",units="px", width=4000, height=2500, res=300, bg="transparent")

pp 
dev.off()
```

    ## png 
    ##   2

# Manhattan Plots

#### CASE

``` r
b10.h <- subset(total.sig,  BLOCK == "10" & Sig.CASE == TRUE)
b10.c <- subset(B10.pv, QCON < 0.1)


B10.pp <- B10.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b10.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B10.pp, pvalue="QCASE", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmo <- semi_join(dfm, b10.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})
b10.CASE <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[4], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[4])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size = FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  scale_size_continuous(range = c(0.25,3)) +
  ggtitle("Block 10") +
  labs(x = "Chromosome") +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 

b11.h <- subset(total.sig,  BLOCK == "11" & Sig.CASE == TRUE)
b11.c <- subset(B11.pv, QCON < 0.1)

B11.pp <- B11.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b11.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B11.pp, pvalue="QCASE", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))


dfmo <- semi_join(dfm, b11.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b11.CASE <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+  
  #geom_point(aes(x=index, y=log10(QCON), colour = as.factor(chrom_alt), alpha=0.5, shape = Group))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[4], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[4])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  scale_size_continuous(range = c(0.25,3)) +
  labs(x = "Chromosome") +
  ggtitle("Block 11") +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 


b12.h <- subset(total.sig,  BLOCK == "12" & Sig.CASE == TRUE)
b12.c <- subset(B12.pv, QCON < 0.1)

B12.pp <- B12.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b12.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B12.pp, pvalue="QCASE", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmo <- semi_join(dfm, b12.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b12.CASE <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[4], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[4])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  #scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  scale_y_continuous(expand = c(0,0),limits=c(0,20))+ #Removes two outliers near end of CHR5
  labs(x = "Chromosome") +
  ggtitle("Block 12") +
  scale_size_continuous(range = c(0.25,3), limits=c(0,20)) +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 

png(filename="CASEmpHighlighted.png", type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")

(b10.CASE / b11.CASE /b12.CASE )  + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') 
#(b10.CASE / b11.CASE /b12.CASE )  + plot_annotation(tag_levels = 'A') & theme_black()
#man+scale_fill_manual(values = c("#009E73","#E69F00", "#0072B2"),name = "Legend")
null <- dev.off()
```

#### SE

``` r
b10.h <- subset(total.sig,  BLOCK == "10" & Sig.SE == TRUE)
b10.c <- subset(B10.pv, QCON < 0.1)

B10.pp <- B10.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b10.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B10.pp, pvalue="QSE", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmo <- semi_join(dfm, b10.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b10.SE <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[3], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[3])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  ggtitle("Block 10") +
  labs(x = "Chromosome") +
  scale_size_continuous(range = c(0.25,3)) +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 

b11.h <- subset(total.sig,  BLOCK == "11" & Sig.SE == TRUE)
b11.c <- subset(B11.pv, QCON < 0.1)

B11.pp <- B11.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b11.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B11.pp, pvalue="QSE", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))


dfmo <- semi_join(dfm, b11.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b11.SE <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+  
  #geom_point(aes(x=index, y=log10(QCON), colour = as.factor(chrom_alt), alpha=0.5, shape = Group))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[3], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[3])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  labs(x = "Chromosome") +
  ggtitle("Block 11") +
  theme_classic()+
  scale_size_continuous(range = c(0.25,3)) +
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 


b12.h <- subset(total.sig,  BLOCK == "12" & Sig.SE == TRUE)
b12.c <- subset(B12.pv, QCON < 0.1)

B12.pp <- B12.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b12.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B12.pp, pvalue="QSE", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmo <- semi_join(dfm, b12.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b12.SE <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[3], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[3])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  #scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  scale_y_continuous(expand = c(0,0),limits=c(0,15))+ #Removes two singleton loci at end of CHR 5
  labs(x = "Chromosome") +
  scale_size_continuous(range = c(0.25,3), limits= c(0,15)) +
  ggtitle("Block 12") +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 

png(filename="SEmpHighlighted.png", type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")

(b10.SE / b11.SE /b12.SE )  + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') 
#(b10.SE / b11.SE /b12.SE )  + plot_annotation(tag_levels = 'A') & theme_black()
#man+scale_fill_manual(values = c("#009E73","#E69F00", "#0072B2"),name = "Legend")
null <- dev.off()
```

#### CA

``` r
b10.h <- subset(total.sig,  BLOCK == "10" & Sig.CA == TRUE)
b10.c <- subset(B10.pv, QCON < 0.1)

B10.pp <- B10.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b10.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B10.pp, pvalue="QCA", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmo <- semi_join(dfm, b10.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b10.CA <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[2], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[2])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  #scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  scale_y_continuous(expand = c(0,0),limits=c(0,32))+ #Limit manually set, removes 3 control filtered loci in CHR 2
  ggtitle("Block 10") +
  labs(x = "Chromosome") +
  scale_size_continuous(range = c(0.25,3), limits=c(0,32)) +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 

b11.h <- subset(total.sig,  BLOCK == "11" & Sig.CA == TRUE)
b11.c <- subset(B11.pv, QCON < 0.1)

B11.pp <- B11.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b11.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B11.pp, pvalue="QCA", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))


dfmo <- semi_join(dfm, b11.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b11.CA <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+
  #geom_point(aes(x=index, y=log10(QCON), colour = as.factor(chrom_alt), alpha=0.5, shape = Group))+  
  geom_point(data=dfmo, fill=cbPaletteSmall4[2], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[2])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  labs(x = "Chromosome") +
  ggtitle("Block 11") +
  scale_size_continuous(range = c(0.25,3)) +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 


b12.h <- subset(total.sig,  BLOCK == "12" & Sig.CA == TRUE)
b12.c <- subset(B12.pv, QCON < 0.1)

B12.pp <- B12.pv %>%
  mutate(
    Group = case_when(
      SNP %in% singleton$SNP ~ "Singleton",
      SNP %in% b12.c$SNP ~ "Control Treatment Filtered",
      SNP %in% all_con$SNP ~ "Control Treatment Filtered",
      SNP %in% con_filter$SNP ~ "Control Treatment Filtered",
      TRUE ~ "Shared SNP"
    )
  )

man <-ggman(B12.pp, pvalue="QCA", chr="CHR", snp="SNP",bp="BP", relative.positions = TRUE ,sigLine = NA, pointSize=1.5,title = "", alpha = 0.5)

dfm <- man[[1]]

dfm$chrom_alt <- factor(dfm$chrom_alt, levels=c(0,1))

dfmo <- semi_join(dfm, b12.h, by = "SNP")
dfmo$Group <- "Outlier"

dfmsplit <- split(dfm, dfm$chrom)
xbreaks <- sapply(dfmsplit,function(x) {
    midpoint <- length(x$index)/2
    if(midpoint <1) midpoint <- 1
    return(x$index[midpoint])
})

b12.CA <-ggplot(dfm, aes(x= index, y=marker, colour = as.factor(chrom_alt)))+
  geom_point(aes(alpha=0.5, shape=Group, size = marker))+ 
  geom_point(data=dfmo, fill=cbPaletteSmall4[2], aes(x=index, y=marker, shape=Group, size=marker),  color = cbPaletteSmall4[2])+
  scale_x_continuous(breaks = xbreaks, labels = names(xbreaks),expand = c(0,0),limits = c(0,max(dfm$index)+10)) +
  guides(colour = FALSE,alpha=FALSE, fill = FALSE, size=FALSE) +
  #scale_y_continuous(expand = c(0,0),limits=c(0,max(dfm$marker+1,na.rm=TRUE)))+
  scale_y_continuous(expand = c(0,0),limits=c(0,20))+ #removes two singleton loci near end of CHR 5
  labs(x = "Chromosome") +
  scale_size_continuous(range = c(0.25,3), limits=c(0,20)) +
  ggtitle("Block 12") +
  theme_classic()+
  labs(y=expression(-Log10(q-value))) +
  scale_shape_manual(values=c(15,21,16,17))+
  scale_color_manual(values=c("grey", "dark grey"), guide =FALSE) 

png(filename="CAmpHighlighted.png", type="cairo",units="px", width=5600, height=3000, res=300, bg="transparent")

(b10.CA / b11.CA /b12.CA )  + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') 
#(b10.CA / b11.CA /b12.CA )  + plot_annotation(tag_levels = 'A') & theme_black()
#man+scale_fill_manual(values = c("#009E73","#E69F00", "#0072B2"),name = "Legend")
null <- dev.off()
```

# Venn Diagrams

``` r
library(VennDiagram)

sig.CA <- nrow(unique(subset(grouped_total.sig, Sig.CA == "TRUE",c(SNP))))
sig.SE <- nrow(unique(subset(grouped_total.sig, Sig.SE == "TRUE",c(SNP))))
sig.CASE <- nrow(unique(subset(grouped_total.sig, Sig.CASE == "TRUE",c(SNP))))

sig.CA.SE <- nrow(unique(subset(grouped_total.sig, Sig.CA == "TRUE" & Sig.SE == "TRUE",c(SNP) )))
sig.CA.CASE <- nrow(unique(subset(grouped_total.sig,Sig.CA == "TRUE" & Sig.CASE == "TRUE",c(SNP) )))

sig.SE.CASE <- nrow(unique(subset(grouped_total.sig, Sig.CASE == "TRUE" & Sig.SE == "TRUE" ,c(SNP))))

sig.all <- nrow(unique(subset(grouped_total.sig, Sig.CASE == "TRUE" & Sig.CA == "TRUE" & Sig.SE == "TRUE",c(SNP))))
#sig.all <- nrow(subset(pp1, Sig.CA == TRUE & Sig.SE == TRUE  & Sig.CASE == TRUE & QPCADAPT < alpha))

png(filename="Ven.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")

venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

``` r
png(filename="VenBO.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("black", 7))
dev.off()
```

    ## png 
    ##   2

### B10 Only

``` r
sig.CA <- nrow(unique(subset(ppsig.B10.01, Sig.CA == TRUE,c(SNP))))
sig.SE <- nrow(unique(subset(ppsig.B10.01, Sig.SE == TRUE ,c(SNP))))
sig.CASE <- nrow(unique(subset(ppsig.B10.01, Sig.CASE == TRUE,c(SNP))))

sig.CA.SE <- nrow(unique(subset(ppsig.B10.01,Sig.CA == TRUE & Sig.SE == TRUE,c(SNP) )))
sig.CA.CASE <- nrow(unique(subset(ppsig.B10.01,Sig.CA == TRUE & Sig.CASE == TRUE ,c(SNP) )))

sig.SE.CASE <- nrow(unique(subset(ppsig.B10.01, Sig.SE == TRUE & Sig.CASE == TRUE ,c(SNP))))

sig.all <- nrow(unique(subset(ppsig.B10.01, Sig.CA == TRUE & Sig.SE == TRUE & Sig.CASE == TRUE,c(SNP))))
#sig.all <- nrow(subset(pp1, Sig.CA == TRUE & Sig.SE == TRUE  & Sig.CASE == TRUE & QPCADAPT < alpha))

png(filename="VenB10.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")

venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

``` r
png(filename="VenBOB10.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("black", 7))
dev.off()
```

    ## png 
    ##   2

### B11 Only

``` r
sig.CA <- nrow(unique(subset(ppsig.B11.01, Sig.CA == TRUE,c(SNP))))
sig.SE <- nrow(unique(subset(ppsig.B11.01, Sig.SE == TRUE,c(SNP))))
sig.CASE <- nrow(unique(subset(ppsig.B11.01, Sig.CASE == TRUE,c(SNP))))

sig.CA.SE <- nrow(unique(subset(ppsig.B11.01, Sig.CA == TRUE & Sig.SE == TRUE,c(SNP) )))
sig.CA.CASE <- nrow(unique(subset(ppsig.B11.01,Sig.CA == TRUE & Sig.CASE == TRUE ,c(SNP) )))

sig.SE.CASE <- nrow(unique(subset(ppsig.B11.01, Sig.SE == TRUE & Sig.CASE == TRUE ,c(SNP))))

sig.all <- nrow(unique(subset(ppsig.B11.01, Sig.CA == TRUE & Sig.SE == TRUE & Sig.CASE == TRUE,c(SNP))))
#sig.all <- nrow(subset(pp1, Sig.CA == TRUE & Sig.SE == TRUE  & Sig.CASE == TRUE & QPCADAPT < alpha))

png(filename="VenB11.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")

venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

``` r
png(filename="VenBOB11.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("black", 7))
dev.off()
```

    ## png 
    ##   2

### B12 Only

``` r
library(VennDiagram)
alpha = 0.01
alpha2 = 0.1

#pp2 <- subset(pp1, QCON > alpha2)
sig.CA <- nrow(unique(subset(ppsig.B11.01, Sig.CA == TRUE,c(SNP))))
sig.SE <- nrow(unique(subset(ppsig.B11.01, Sig.SE == TRUE,c(SNP))))
sig.CASE <- nrow(unique(subset(ppsig.B11.01, Sig.CASE == TRUE,c(SNP))))

sig.CA.SE <- nrow(unique(subset(ppsig.B11.01, Sig.CA == TRUE & Sig.SE == TRUE,c(SNP) )))
sig.CA.CASE <- nrow(unique(subset(ppsig.B11.01,Sig.CA == TRUE & Sig.CASE == TRUE ,c(SNP) )))

sig.SE.CASE <- nrow(unique(subset(ppsig.B11.01, Sig.SE == TRUE & Sig.CASE == TRUE ,c(SNP))))

sig.all <- nrow(unique(subset(ppsig.B11.01, Sig.CA == TRUE & Sig.SE == TRUE & Sig.CASE == TRUE,c(SNP))))

png(filename="VenB12.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")

venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

``` r
png(filename="VenBOB12.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("black", 7))
dev.off()
```

    ## png 
    ##   2

### LOC Only

``` bash
source activate CASE

cat Sig.loci.1.CA.LOC Sig.loci.1.SE.LOC | sort | uniq -c | mawk '$1 > 1'  | mawk '{print $2}' > Sig.loci.1.CA.SE.LOC
cat Sig.loci.1.CA.LOC Sig.loci.1.CASE.LOC | sort | uniq -c | mawk '$1 > 1'| mawk '{print $2}' > Sig.loci.1.CA.CASE.LOC
cat Sig.loci.1.SE.LOC Sig.loci.1.CASE.LOC | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}'> Sig.loci.1.SE.CASE.LOC
cat Sig.loci.1.SE.LOC Sig.loci.1.CASE.LOC Sig.loci.1.CA.LOC | sort | uniq -c | mawk '$1 > 2'| mawk '{print $2}' > Sig.loci.1.SE.CASE.CA.LOC

cat Sig.loci.1.CASE.LOC <(cat Sig.loci.1.SE.LOC Sig.loci.1.CASE.LOC Sig.loci.1.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.CASE.ONLY.LOC
cat Sig.loci.1.SE.LOC <(cat Sig.loci.1.SE.LOC Sig.loci.1.CASE.LOC Sig.loci.1.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.SE.ONLY.LOC
cat Sig.loci.1.CA.LOC <(cat Sig.loci.1.SE.LOC Sig.loci.1.CASE.LOC Sig.loci.1.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.CA.ONLY.LOC

grep -v -f Sig.loci.1.CA.LOC Sig.loci.1.SE.CASE.LOC > Sig.loci.CASE.SE.ONLY.LOC



cat Sig.loci.B10.CA.LOC Sig.loci.B10.SE.LOC | sort | uniq -c | mawk '$1 > 1'  | mawk '{print $2}' > Sig.loci.B10.CA.SE.LOC
cat Sig.loci.B10.CA.LOC Sig.loci.B10.CASE.LOC | sort | uniq -c | mawk '$1 > 1'| mawk '{print $2}' > Sig.loci.B10.CA.CASE.LOC
cat Sig.loci.B10.SE.LOC Sig.loci.B10.CASE.LOC | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}'> Sig.loci.B10.SE.CASE.LOC
cat Sig.loci.B10.SE.LOC Sig.loci.B10.CASE.LOC Sig.loci.B10.CA.LOC | sort | uniq -c | mawk '$1 > 2'| mawk '{print $2}' > Sig.loci.B10.SE.CASE.CA.LOC

cat Sig.loci.B10.CASE.LOC <(cat Sig.loci.B10.SE.LOC Sig.loci.B10.CASE.LOC Sig.loci.B10.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B10.CASE.ONLY.LOC
cat Sig.loci.B10.SE.LOC <(cat Sig.loci.B10.SE.LOC Sig.loci.B10.CASE.LOC Sig.loci.B10.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B10.SE.ONLY.LOC
cat Sig.loci.B10.CA.LOC <(cat Sig.loci.B10.SE.LOC Sig.loci.B10.CASE.LOC Sig.loci.B10.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B10.CA.ONLY.LOC

grep -v -f Sig.loci.B10.CA.LOC Sig.loci.B10.SE.CASE.LOC > Sig.loci.B10.CASE.SE.ONLY.LOC

cat Sig.loci.B11.CA.LOC Sig.loci.B11.SE.LOC | sort | uniq -c | mawk '$1 > 1'  | mawk '{print $2}' > Sig.loci.B11.CA.SE.LOC
cat Sig.loci.B11.CA.LOC Sig.loci.B11.CASE.LOC | sort | uniq -c | mawk '$1 > 1'| mawk '{print $2}' > Sig.loci.B11.CA.CASE.LOC
cat Sig.loci.B11.SE.LOC Sig.loci.B11.CASE.LOC | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}'> Sig.loci.B11.SE.CASE.LOC
cat Sig.loci.B11.SE.LOC Sig.loci.B11.CASE.LOC Sig.loci.B11.CA.LOC | sort | uniq -c | mawk '$1 > 2'| mawk '{print $2}' > Sig.loci.B11.SE.CASE.CA.LOC

cat Sig.loci.B11.CASE.LOC <(cat Sig.loci.B11.SE.LOC Sig.loci.B11.CASE.LOC Sig.loci.B11.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B11.CASE.ONLY.LOC
cat Sig.loci.B11.SE.LOC <(cat Sig.loci.B11.SE.LOC Sig.loci.B11.CASE.LOC Sig.loci.B11.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B11.SE.ONLY.LOC
cat Sig.loci.B11.CA.LOC <(cat Sig.loci.B11.SE.LOC Sig.loci.B11.CASE.LOC Sig.loci.B11.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B11.CA.ONLY.LOC

grep -v -f Sig.loci.B11.CA.LOC Sig.loci.B11.SE.CASE.LOC > Sig.loci.B11.CASE.SE.ONLY.LOC

cat Sig.loci.B12.CA.LOC Sig.loci.B12.SE.LOC | sort | uniq -c | mawk '$1 > 1'  | mawk '{print $2}' > Sig.loci.B12.CA.SE.LOC
cat Sig.loci.B12.CA.LOC Sig.loci.B12.CASE.LOC | sort | uniq -c | mawk '$1 > 1'| mawk '{print $2}' > Sig.loci.B12.CA.CASE.LOC
cat Sig.loci.B12.SE.LOC Sig.loci.B12.CASE.LOC | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}'> Sig.loci.B12.SE.CASE.LOC
cat Sig.loci.B12.SE.LOC Sig.loci.B12.CASE.LOC Sig.loci.B12.CA.LOC | sort | uniq -c | mawk '$1 > 2'| mawk '{print $2}' > Sig.loci.B12.SE.CASE.CA.LOC

cat Sig.loci.B12.CASE.LOC <(cat Sig.loci.B12.SE.LOC Sig.loci.B12.CASE.LOC Sig.loci.B12.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B12.CASE.ONLY.LOC
cat Sig.loci.B12.SE.LOC <(cat Sig.loci.B12.SE.LOC Sig.loci.B12.CASE.LOC Sig.loci.B12.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B12.SE.ONLY.LOC
cat Sig.loci.B12.CA.LOC <(cat Sig.loci.B12.SE.LOC Sig.loci.B12.CASE.LOC Sig.loci.B12.CA.LOC | sort | uniq -c | mawk '$1 < 2' | mawk '{print $2}') | sort | uniq -c | mawk '$1 > 1' | mawk '{print $2}' > Sig.loci.B12.CA.ONLY.LOC

grep -v -f Sig.loci.B12.CA.LOC Sig.loci.B12.SE.CASE.LOC > Sig.loci.B12.CASE.SE.ONLY.LOC
```

``` r
library(VennDiagram)
alpha = 0.05
alpha2 = 0.1

Sig.loci.1.CA.LOC <- read.table("Sig.loci.1.CA.LOC", header =FALSE)
Sig.loci.1.SE.LOC <- read.table("Sig.loci.1.SE.LOC", header =FALSE)
Sig.loci.1.CASE.LOC <- read.table("Sig.loci.1.CASE.LOC", header =FALSE)

sig.CA <- nrow(unique(Sig.loci.1.CA.LOC))
sig.SE <- nrow(unique(Sig.loci.1.SE.LOC))
sig.CASE <- nrow(unique(Sig.loci.1.CASE.LOC))

sig.CA.SE <- nrow(unique(read.table("Sig.loci.1.CA.SE.LOC", header =FALSE)))
sig.CA.CASE <- nrow(unique(read.table("Sig.loci.1.CA.CASE.LOC", header =FALSE)))
sig.SE.CASE <- nrow(unique(read.table("Sig.loci.1.SE.CASE.LOC", header =FALSE)))

sig.all <- nrow(unique(read.table("Sig.loci.1.SE.CASE.CA.LOC", header =FALSE)))


png(filename="VenLOC.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

### B10

``` r
library(VennDiagram)
alpha = 0.05
alpha2 = 0.1

Sig.loci.B10.CA.LOC <- read.table("Sig.loci.B10.CA.LOC", header =FALSE)
Sig.loci.B10.SE.LOC <- read.table("Sig.loci.B10.SE.LOC", header =FALSE)
Sig.loci.B10.CASE.LOC <- read.table("Sig.loci.B10.CASE.LOC", header =FALSE)

sig.CA <- nrow(unique(Sig.loci.B10.CA.LOC))
sig.SE <- nrow(unique(Sig.loci.B10.SE.LOC))
sig.CASE <- nrow(unique(Sig.loci.B10.CASE.LOC))

sig.CA.SE <- nrow(unique(read.table("Sig.loci.B10.CA.SE.LOC", header =FALSE)))
sig.CA.CASE <- nrow(unique(read.table("Sig.loci.B10.CA.CASE.LOC", header =FALSE)))
sig.SE.CASE <- nrow(unique(read.table("Sig.loci.B10.SE.CASE.LOC", header =FALSE)))

sig.all <- nrow(unique(read.table("Sig.loci.B10.SE.CASE.CA.LOC", header =FALSE)))


png(filename="VenB10LOC.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

### B11

``` r
library(VennDiagram)
alpha = 0.05
alpha2 = 0.1

Sig.loci.B11.CA.LOC <- read.table("Sig.loci.B11.CA.LOC", header =FALSE)
Sig.loci.B11.SE.LOC <- read.table("Sig.loci.B11.SE.LOC", header =FALSE)
Sig.loci.B11.CASE.LOC <- read.table("Sig.loci.B11.CASE.LOC", header =FALSE)

sig.CA <- nrow(unique(Sig.loci.B11.CA.LOC))
sig.SE <- nrow(unique(Sig.loci.B11.SE.LOC))
sig.CASE <- nrow(unique(Sig.loci.B11.CASE.LOC))

sig.CA.SE <- nrow(unique(read.table("Sig.loci.B11.CA.SE.LOC", header =FALSE)))
sig.CA.CASE <- nrow(unique(read.table("Sig.loci.B11.CA.CASE.LOC", header =FALSE)))
sig.SE.CASE <- nrow(unique(read.table("Sig.loci.B11.SE.CASE.LOC", header =FALSE)))

sig.all <- nrow(unique(read.table("Sig.loci.B11.SE.CASE.CA.LOC", header =FALSE)))


png(filename="VenB11LOC.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

### B12

``` r
library(VennDiagram)
alpha = 0.05
alpha2 = 0.1

Sig.loci.B12.CA.LOC <- read.table("Sig.loci.B12.CA.LOC", header =FALSE)
Sig.loci.B12.SE.LOC <- read.table("Sig.loci.B12.SE.LOC", header =FALSE)
Sig.loci.B12.CASE.LOC <- read.table("Sig.loci.B12.CASE.LOC", header =FALSE)

sig.CA <- nrow(unique(Sig.loci.B12.CA.LOC))
sig.SE <- nrow(unique(Sig.loci.B12.SE.LOC))
sig.CASE <- nrow(unique(Sig.loci.B12.CASE.LOC))

sig.CA.SE <- nrow(unique(read.table("Sig.loci.B12.CA.SE.LOC", header =FALSE)))
sig.CA.CASE <- nrow(unique(read.table("Sig.loci.B12.CA.CASE.LOC", header =FALSE)))
sig.SE.CASE <- nrow(unique(read.table("Sig.loci.B12.SE.CASE.LOC", header =FALSE)))

sig.all <- nrow(unique(read.table("Sig.loci.B12.SE.CASE.CA.LOC", header =FALSE)))


png(filename="VenB12LOC.png", type="cairo",units="px", width=3000, height=3000, res=200, bg="transparent")
venn.plot <- draw.triple.venn(
  area1=sig.CA, area2=sig.SE, area3=sig.CASE, 
  n12= sig.CA.SE, n13=sig.CA.CASE, 
  n23=sig.SE.CASE,
  n123=sig.all, cex=5,
  category = c("CA", "SE", "CASE"),
  fill = cbPaletteSmall3,cat.cex = rep(5, 3),
  cat.col = cbPaletteSmall3, label.col = rep("white", 7))
dev.off()
```

    ## png 
    ##   2

# PCA

``` bash
source activate CASE



bcftools view --threads 40 ../raw.vcf/B10.CASE.FIL.vcf.gz -R <(cut -f2,3 Total.Significant.Loci | tail -n +2 | sort | uniq )| mawk '!/\.:\.:\./' > B10.total.outlier.CASE.dp20.vcf &

bcftools view --threads 40 ../raw.vcf/B11.CASE.FIL.vcf.gz -R <(cut -f2,3 Total.Significant.Loci | tail -n +2 | sort | uniq )| mawk '!/\.:\.:\./' > B11.total.outlier.CASE.dp20.vcf &

bcftools view --threads 40 ../raw.vcf/B11.CASE.FIL.vcf.gz -R <(mawk '$16 == 11 || $16 ==12' Total.Significant.Loci | cut -f2,3 | sort | uniq )| mawk '!/\.:\.:\./' > B11.B11.outlier.CASE.dp20.vcf

bcftools view --threads 40 ../raw.vcf/B12.CASE.FIL.vcf.gz -R <(cut -f2,3 Total.Significant.Loci | tail -n +2 | sort | uniq )| mawk '!/\.:\.:\./' > B12.total.outlier.CASE.dp20.vcf &

bcftools view --threads 40 -S <(grep -v 1.G ../raw.vcf/samples | grep -v J17B10 | grep -v J05B10) ../raw.vcf/CASE.TRSdp.20.g5.nDNA.FIL.vcf.gz -R <(cut -f2,3 Total.Significant.Loci | tail -n +2 | sort | uniq )| mawk '!/\.:\.:\./' > total.outlier.CASE.dp20.vcf 

wait

python2 ~/CASE/VCFtoPopPool.py B10.total.outlier.CASE.dp20.vcf B10.total.outlier.sync
python2 ~/CASE/VCFtoPopPool.py B11.total.outlier.CASE.dp20.vcf B11.total.outlier.sync
python2 ~/CASE/VCFtoPopPool.py B12.total.outlier.CASE.dp20.vcf B12.total.outlier.sync 

python2 ~/CASE/VCFtoPopPool.py B11.B11.outlier.CASE.dp20.vcf B11.B11.outlier.sync


python2 ~/CASE/VCFtoPopPool.py total.outlier.CASE.dp20.vcf AB.total.outlier.sync


mawk '!/CHR/' AB.total.outlier.sync > AB.total.outlier.input
mawk '!/CHR/' B10.total.outlier.sync > B10.total.outlier.input
mawk '!/CHR/' B11.total.outlier.sync > B11.total.outlier.input
mawk '!/CHR/' B11.B11.outlier.sync > B11.B11.outlier.input
mawk '!/CHR/' B12.total.outlier.sync > B12.total.outlier.input

../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input  B10.total.outlier.input --min-count 1 --min-coverage 1 --output-prefix B10.total.outlier --max-coverage 50000 &

../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input  B11.total.outlier.input --min-count 1 --min-coverage 1 --output-prefix B11.total.outlier --max-coverage 50000 &

../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input  B11.B11.outlier.input --min-count 1 --min-coverage 1 --output-prefix B11.B11.outlier --max-coverage 50000 &

../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input  B12.total.outlier.input --min-count 1 --min-coverage 1 --output-prefix B12.total.outlier --max-coverage 50000 &

../scripts/assessPool/scripts/p2/snp-frequency-diff.pl --input  AB.total.outlier.input --min-count 1 --min-coverage 1 --output-prefix AB.total.outlier --max-coverage 50000 

mawk -f ../scripts/polarize_freqs AB.total.outlier_rc | cut -f10-65 | mawk '!/maa/' > AB.total.outlier.CASE.dp20.pool

mawk -f ../scripts/polarize_freqs B10.total.outlier_rc | cut -f10-25 | mawk '!/maa/' > B10.total.outlier.CASE.dp20.pool

mawk -f ../scripts/polarize_freqs B11.total.outlier_rc | cut -f10-29 | mawk '!/maa/' > B11.total.outlier.CASE.dp20.pool

mawk -f ../scripts/polarize_freqs B11.B11.outlier_rc | cut -f10-29 | mawk '!/maa/' > B11.B11.outlier.CASE.dp20.pool

mawk -f ../scripts/polarize_freqs B12.total.outlier_rc | cut -f10-29 | mawk '!/maa/' > B12.total.outlier.CASE.dp20.pool
```

``` r
pool.data <- read.table("AB.total.outlier.CASE.dp20.pool")
pool.data.b10 <- read.table("B10.total.outlier.CASE.dp20.pool")
pool.data.b11 <- read.table("B11.total.outlier.CASE.dp20.pool")
pool.data.b11o <- read.table("B11.B11.outlier.CASE.dp20.pool")
pool.data.b12 <- read.table("B12.total.outlier.CASE.dp20.pool")

df.pool <- apply(pool.data, c(1, 2), function(x) eval(parse(text = x)))
df.pool.b10 <- apply(pool.data.b10, c(1, 2), function(x) eval(parse(text = x)))
df.pool.b11 <- apply(pool.data.b11, c(1, 2), function(x) eval(parse(text = x)))
df.pool.b11o <- apply(pool.data.b11o, c(1, 2), function(x) eval(parse(text = x)))
df.pool.b12 <- apply(pool.data.b12, c(1, 2), function(x) eval(parse(text = x)))

pool.data2 <- t(df.pool)
pool.data2.b10 <- t(df.pool.b10)
pool.data2.b11 <- t(df.pool.b11)
pool.data2.b11o <- t(df.pool.b11o)
pool.data2.b12 <- t(df.pool.b12)


filename <- read.pcadapt(pool.data2, type = "pool")
res <- pcadapt(filename, min.maf = 0.1)

filename.b10 <- read.pcadapt(pool.data2.b10, type = "pool")

filename.b11 <- read.pcadapt(pool.data2.b11, type = "pool")
filename.b11o <- read.pcadapt(pool.data2.b11o, type = "pool")

filename.b12 <- read.pcadapt(pool.data2.b12, type = "pool")



par(mfrow = c(2, 2))
for (i in 1:4)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

``` r
plot(res,option="screeplot")
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

``` r
res <- pcadapt(filename, K =4,min.maf = 0.001)
poplist.names <- c(rep("CA", 11),rep("CASE", 11),rep("CON", 11),rep("IS", 12),rep("SE", 11))
p1 <- plot(res, option = "scores", i = 1, j = 2, pop = poplist.names)
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

``` r
p1.df <- data.frame(res$scores)
colnames(p1.df) <- c("PC1","PC2", "PC3", "PC4")
p1.df$POP <- poplist.names
png(filename="PC1.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatmen")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

``` r
plot(res, option = "scores", i = 1, j = 2, pop = poplist.names)
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-68-2.png)<!-- -->

``` r
poplist.names <- c("B12","B11","B12","B10","B11","B10","B11","B12","B10","B12","B11","B10","B11","B12","B12","B10","B12","B11","B10","B11","B11","B12","B10","B11","B12","B11","B10","B12","B10","B11","B12","B11","B12","B10","B11","B10","B11","B10","B11","B10","B11","B12","B12","B12","B12","B11","B12","B10","B10","B11","B11","B10","B12","B12","B11","B12")
p1 <- plot(res, option = "scores", i = 1, j = 2, pop = poplist.names)
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-68-3.png)<!-- -->

``` r
p1.df <- data.frame(res$scores)
colnames(p1.df) <- c("PC1","PC2", "PC3", "PC4")
p1.df$POP <- poplist.names
png(filename="PC1b.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =10,shape=21, col="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Spawn")+ 
guides(alpha=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

``` r
res.b10 <- pcadapt(filename.b10, K = 3, min.maf=0.00)
poplist.names <- c(rep("CA", 3),rep("CASE", 3),rep("CON", 3),rep("IS", 4),rep("SE", 3))
p1 <- plot(res.b10, option = "scores", i = 1, j = 2, pop = poplist.names)
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

``` r
p1.df <- data.frame(res.b10$scores)
colnames(p1.df) <- c("PC1","PC2", "PC3")
p1.df$POP <- poplist.names
png(filename="PCb10.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent") 
ggplot(p1.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =12, shape =21, color="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ guides(alpha=FALSE, fill=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

``` r
res.b11 <- pcadapt(filename.b11, K = 3, min.maf=0.01)
poplist.names <- c(rep("CA", 4),rep("CASE", 4),rep("CON", 4),rep("IS", 4),rep("SE", 4))
p1 <- plot(res.b11, option = "scores", i = 1, j = 2, pop = poplist.names)
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-69-2.png)<!-- -->

``` r
p1.df <- data.frame(res.b11$scores)
p1.df <- p1.df[-2,]
colnames(p1.df) <- c("PC1","PC2", "PC3")
p1.df$POP <- poplist.names[-2]
png(filename="PCb11.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.df, aes(x=PC1, y= PC2, fill = POP)) + geom_point(aes(alpha=0.2), shape=21, size =12, color="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ guides(alpha=FALSE,fill=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

``` r
res.b11o <- pcadapt(filename.b11o, K = 3, min.maf=0.01)
poplist.names <- c(rep("CA", 4),rep("CASE", 4),rep("CON", 4),rep("IS", 4),rep("SE", 4))
p1 <- plot(res.b11o, option = "scores", i = 1, j = 2, pop = poplist.names)
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-69-3.png)<!-- -->

``` r
p1.df <- data.frame(res.b11o$scores)
p1.df <- p1.df[-2,]
colnames(p1.df) <- c("PC1","PC2", "PC3")
p1.df$POP <- poplist.names[-2]
png(filename="PCb11o.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.df, aes(x=PC1, y= PC2, fill = POP)) + geom_point(aes(alpha=0.2), shape=21, size =12, color="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ guides(alpha=FALSE,fill=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2

``` r
res.b12 <- pcadapt(filename.b12, K = 3, min.maf=0.0)
poplist.names <- c(rep("CA", 4),rep("CASE", 4),rep("CON", 4),rep("IS", 4),rep("SE", 4))
p1 <- plot(res.b12, option = "scores", i = 1, j = 2, pop = poplist.names)
```

![](CASE_FULL_FINAL_files/figure-gfm/unnamed-chunk-69-4.png)<!-- -->

``` r
p1.df <- data.frame(res.b12$scores)
colnames(p1.df) <- c("PC1","PC2", "PC3")
p1.df$POP <- poplist.names
png(filename="PCb12.png", type="cairo",units="px", width=5400, height=3000, res=300, bg="transparent")
ggplot(p1.df, aes(x=PC1, y= PC2, fill= POP)) + geom_point(aes(alpha=0.2), size =12,shape=21, color="white")+ theme_black() + scale_fill_manual(values=cbPaletteSmall,name="Treatment")+ guides(alpha=FALSE, fill=FALSE) +theme(axis.title.x = element_text(size = 24),axis.title.y = element_text(size = 24))
dev.off()
```

    ## png 
    ##   2
