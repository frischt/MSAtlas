<div style="max-width: 65%; margin-left: auto; margin-right: auto">
  
  <h1>Step-by-step guide</h1>

  <div class="well toc">
    <p><b>Contents</b></p>
    <ul>
      <li><a href="#paramSelection">Parameter Selection</a></li>
      <ul>
        <li><a href="#lesionType">Selection of Lesion Type</a></li>
        <li><a href="#filtering">Filtering by FDR and Log<sub>2</sub>FC</a></li>
        <li><a href="#network_enrichment">Network Enrichment</a></li>
        <li><a href="#geneSymbols">Gene Symbols & Export</a></li>
      </ul>
      <li><a href="#res">Results</a></li>
      <ul>
        <li><a href="#heatmap">Heatmap</a></li>
        <li><a href="#network">Network</a></li>
        <li><a href="#volcano">Volcano</a></li>
      </ul>
    </ul>
</div>
  <h2 id="paramSelection">Parameter Selection</h2> 

  The database contains all genes as well as the corresponding Log<sub>2</sub>FC and p-Values (FDR corrected) from the statistical analysis. Depending on the lesion type of interest, and the thresholds you consider to be significant, you can select a subset of all genes. In the following section we describe the parameter you can select and elaborate on the consequences. In the first step of the standard workflow you will start by selecting the set of genes you are interested in. Afterwards, you can select the parameters for <i>de novo</i> network enrichment and visualize the resulting IID subnetwork.


  <h3 id="lesionType", style="margin-top: 100px; margin-bottom: 25px", align="center">Selection of Lesion type</h3>
  <div class="row">
  	<div class="col-lg-6">
  		In the first step, you can select the lesion type of interest. It is possible to filter for up-, down-, or overall de-regulated genes. For example, you might want to look at genes that are up-regulated in active lesion and, at the same time, down-regulated in chronic active lesion.
   </div>
    <div class="col-lg-6", align="center"><img src=images/sidePane_1.png></div>
  </div>


  <h3 id="filtering", style="margin-top: 50px; margin-bottom: 25px", align="center"">Filtering by FDR and Log<sub>2</sub>FC</h3>
  <div class="row">
    <div class="col-lg-6", align="center", style="margin=auto"><img src=images/sidePane_2.png></div>
  	<div class="col-lg-6", style="margin=auto">
      In the second filtering step you can select the p-value (FDR corrected) and Log<sub>2</sub>FC. By default, you filter for a FDR < 0.05 allowing 5 percent of the remaining genes to be false positive. Additionally, it is possible to adjust the Log<sub>2</sub>FC if you are interested only in genes showing a high fold-change between control and lesion type.
      Finally, you have to decide if the selected parameters should be applied to only one or all the selected lesion types. The option "all", requires the genes to fulfill the criteria in all selected lesion types, which is much more strict.
  </div>
  </div>

  <h3 id="geneSymbols", style="margin-top: 50px; margin-bottom: 25px", align="center"">Gene Symbols & Export</h3>
  <div class="row">
    <div class="col-lg-6">
      You have the possibility to select genes of interest based on a list of gene symbols. You can either insert them comma separated into the text field or upload a CSV file (<a href="geneList.csv" download>example file</a>). You should be aware, that this is the last filter applied to the gene list. If none of the genes you enter here are significant regarding the filtering steps described above, you might end up without any significant genes. By pressing the "Start" button the selection is confirmed, all parameters are applied, and the results are visualized. If required, the gene sets that fulfill your criteria can be exported into a CSV file, including the gene names and the Log<sub>2</sub>FC for the corresponding lesion.
    </div>
    <div class="col-lg-6", align="center"><img src=images/sidePane_3.png></div>
  </div>

  <h3 id="heatmap", style="margin-top: 50px; margin-bottom: 25px", align="center">Visualization of Selected Genes</h3>
  By pressing the start button, all selected parameters are applied, and the result is consequently visualized in a heatmap (see below). The rows and columns are ordered based on a Hierarchical Clustering using Euclidean distance. Additionally, on the columns a dendrogram shows which lesion types cluster closer together. The color coding is based on the Log<sub>2</sub>FC, were the color red indicates genes that are upregulated in the lesion, while blue genes are downregulated. When hovering over the heatmap with the mouse, the gene name and the corresponding Log<sub>2</sub>FC is shown.
  <div align="center">
    <img src="images/screen_heatmap.png" alt="heatmap", style="max-width: 100%; margin-top: 25px">
  </div>

 <h2 id="network_enrichment", style="margin-top: 75">Network Enrichment</h2> 
  Next, you can push on the bottom "Network", and visualize the selected genes projected onto a human brain specific protein-protein interaction network. Here, the integrated <i>de novo</i> network enrichment method <a href="https://keypathwayminer.compbio.sdu.dk/keypathwayminer/">KeyPathwayMiner</a>, which extracts sub-network that distinguish on a mechanistic level between MS lesion types and, thus, provides first hints on how lesion evolution is driven and controlled on a system biological level. You can select the number of exception genes (k), which do not necessarily have to be significantly differentially expressed between lesion types (i.e. outliers), but still play a central role in the interaction network. A mouse klick on a node/gene in the network reveals additional information. The key networks can be exported in SIF format for downstream analyses in Cytoscape or as PNG image file.


  <h3 id="network", style="margin-top: 50px; margin-bottom: 25px", align="center">Network</h3>
   <div class="row", style="margin-bottom: 25px">
    <div class="col-lg-6">
     In order to extract useful information out of the sparsely connected network you can either visualize the biggest connected component or run <a href="https://keypathwayminer.compbio.sdu.dk/keypathwayminer/">KeyPathwayMiner</a>. When you are interested in the connected components you'll have to option to select the first, second,... connected component. For <a href="https://keypathwayminer.compbio.sdu.dk/keypathwayminer/">KeyPathwayMiner</a> you have to select the parameter k indicating how many outliers you want to allow.
   </div>
    <div class="col-lg-6", align="center"><img src=images/sidePane_2_1.png></div>
  </div>

  The resulting network is shown in the figure below. Every node represents one gene and the corresponding name will be shown when you hover your mouse above the node. The color coding shows betweenness centrality where red indicates nodes that are more central and hence more important for the network.
  <img src="images/screen_network.png" alt="network", style="max-width: 100%; margin-top: 25px">

  <h3 id="volcano", style="margin-top: 50px; margin-bottom: 25px", align="center">Volcano</h3>
  Finally, you can visualize volcano plots for all genes of a selected lesion. After selecting the lesion of interest, you will get an on-the-fly visualization of the transcriptional landscape of genes detected between that lesion type and control. Deregulated genes below FDR 0.05 are indicated in bright red and roange, where orange also indicated absolute Log<sub>2</sub>FC greater 1.

  <img src="images/screen_vulcano.png" alt="heatmap", style="max-width: 100%; margin-top: 25px">

</div>

