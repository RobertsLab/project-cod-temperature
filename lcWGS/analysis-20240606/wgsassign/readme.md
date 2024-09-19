WGSassign is a program that uses a pre-determined set of markers/sites/SNPs that can predict population structure. Pipeline steps include: 

- Create training & test sets of bam lists from reference fish of known populations  
- Downsample populations with very high number of reference fish (as that can sway assignments - more likely to assign fish to those populations)    
- get site allele frequences (SAFs), and for each pairwise population combination calculate 2-dimentional site frequency spectrum (2dSFS) & fixation index (Fst)     
- Create multiple subsets of SNPs from the top N snps from each pairwise population combinaton based on (highest) Fst (e.g. top 50, 100, 1000 from each)  
- Assign test reference individuals and assess accuracy to identify how many top n SNPs to use for experimental assignment    
- Assign experimental fish to marine regions  

I do these steps twice, 1) using all spawning locations at the finest resolution, and 2) grouping spawning locations into 5 marine reigions  
