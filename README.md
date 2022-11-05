
# MetDNA2Vis

<!-- badges: start -->
<!-- badges: end -->

[**Knowledge-guided multi-layer network (KGMN)**](http://metdna.zhulab.cn/) is a new approach leveraging knowledge-guided multi-layer networks to annotate known and unknown metabolites in untargeted metabolomics data. The goal of **MetDNA2Vis** package is to visualize, reproduce and investigate putatively annotated known and unknown metabolites from KGMN. <br>

**Detailed tutorial** can be found [**here**](https://github.com/ZhuMetLab/MetDNA2_Web/blob/main/Tutorials/Tutorial_visualization.pdf)


## Installation

You can install the development version of MetDNA2InSilicoTool like so:

``` r
if(!require(devtools)){
install.packages("devtools")
}

if(!require(BiocManager)){
install.packages("BiocManager")
}

# Install CRAN/Bioconductor packages
required_pkgs <- c("dplyr","tidyr","readr","CHNOSZ","igraph",
  "magrittr","ggplot2","ggraph","tidygraph")
list_installed <- installed.packages()

new_pkgs <- required_pkgs[!(required_pkgs %in% list_installed[,'Package'])]
if (length(new_pkgs) > 0) {
  BiocManager::install(new_pkgs)
} else {
  cat('Required CRAN/Bioconductor packages installed\n')
}


# Install ZhuLab packages
devtools::install_github("ZhuMetLab/SpectraTools")
devtools::install_github("ZhuMetLab/MetDNA2Vis")
```

## Example

Here is a example which contains codes to help to reproduce above analysis quickly.

``` r
library(MetDNA2Vis)
library(SpectraTools)

# set working directory
setwd('D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/Demo_MetDNA2_NIST_urine_pos/06_visualization/')

# Export global networks 
# construct network 1
reconstructNetwork1(is_unknown_annotation = TRUE)

# construct network 2
annotation_table <- reformatTable1()
reconstructNetwork2(annotation_table = annotation_table)

# construct network 3
reconstructNetwork3()

# Export subnetworks -----------------------------------------------------------
# network 1 of unknown peak subnetwork
# Note: the folder_output should keep same among different layer subnetworks
retrieveSubNetwork1(centric_met = c('C00082', 'KeggExd000923'), 
  is_unknown_annotation = TRUE, 
  folder_output = c('M182T541_M262T526'))


# network 2 of unknown peak subnetwork
retrieveSubNetwork2(from_peak = 'M182T541', 
  end_peak = 'M262T526', 
  folder_output = c('M182T541_M262T526'))

# network 3 of unknown peak subnetwork
retrieveSubNetwork3(base_peaks = c('M182T541', 'M262T526'),
  base_adducts = c('[M+H]+', '[M+H]+'),
  folder_output = c('M182T541_M262T526'))


# merge subnetwork
mergeSubnetwork(from_peak = 'M182T541', 
  end_peak = 'M262T526', 
  folder_output = 'M182T541_M262T526')

```

## Citation
This free open-source software implements academic research by the authors and co-workers. If you use it, please support the project by citing the appropriate journal articles. 

Zhiwei Zhou†, Mingdu Luo†, Haosong Zhang, Yandong Yin, Yuping Cai, and Zheng-Jiang Zhu*, Metabolite annotation from knowns to unknowns through knowledge-guided multi-layer metabolic networking, **Nature Communications**, **2022**, 13: 6656 [**Link**](https://www.nature.com/articles/s41467-022-34537-6)

## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)

