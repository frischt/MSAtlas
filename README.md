# MS Atlas
The MS Atlas is a database providing transcriptome information of multiple sclerosis (MS) white matter lesion. We include NAWM, active, chronic active, inactive and remyelinating lesions.

The study is based on post mortem brain tissue from 10 MS and 5 control (non-neurological disease) patients. Altogether, 100 samples were classified, RNA was extracted and next-generation sequencing was applied. The resulting sequencing data were analyzed with edgeR.

Based on the normalized read count, a generalized linear model (accounting for age, sex and lesion distribution) was trained for every lesion type. Finally, the calculated p-values were normalized with FDR-correction (Benjamini-Hochberg).

The MS Atlas is able to visualize differentially expressed genes and extract mechanistic markers based on de novo network enrichment (KeyPathwayMiner). 

# docker build

The webservice can be build and run with docker, based on the following commands:

1) ```docker build -t msatlas .```
2) ```docker run -p 3844:3838 msatlas``` 

Afterwards, the website can be opened in a browser ```http://127.0.0.1:3844/```. Please replace the local port 3844 with any free port on your machine. 
