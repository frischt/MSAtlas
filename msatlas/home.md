<div style="max-width: 65% ;margin-left: auto; margin-right: auto">
<div class="row">
<div class="col-sm-8">
	<h1>MS Atlas</h1>

<div id="bereich">
<h3>Database Update</h3>
<p>
On July 10th 2019 we updated our database, thereby introducing substantial changes to the underlying data. Any analysis performed before this date is invalid and has to be repeated.
</p>
</div>

<p> The MS Atlas is a database providing transcriptome information of multiple sclerosis (MS) white matter lesion. We include NAWM, active, chronic active, inactive and remyelinating lesions.</p>

<p> The study is based on <i>post mortem</i> brain tissue from 10 MS and 5 control (non-neurological disease) patients. Altogether, 100 samples were classified, RNA was extracted and next-generation sequencing was applied. The resulting sequencing data were analyzed with edgeR.</p>

<p> Based on the normalized read count, a generalized linear model (accounting for age, sex and lesion distribution) was trained for every lesion type. Finally, the calculated p-values were normalized with FDR-correction (Benjamini-Hochberg). </p>

<p> The MS Atlas is able to visualize differentially expressed genes and extract mechanistic markers based on <i>de novo</i> network enrichment (KeyPathwayMiner). </p>

<p> For more information how to use the MS Atlas we recommend to watch the screencast you can find on the rigth side. Under "Guide" you can find a similar tutorial in text form. </p>

</div>
<div class="col-sm-4">
<h4>Screencast</h4>
<div class="embed-responsive embed-responsive-16by9">
<iframe class="embed-responsive-item" width="560" height="315" src="https://www.youtube.com/embed/HUfMFrnCIJc" frameborder="0" allowfullscreen>
</iframe>
</div>
</div>
</div>


<div class="row", style="margin-bottom: 100px;">
<div class="col-sm-4 col-xs-4"><a href="images/screen_heatmap.png" target="_blank" class="thumbnail"><img src="images/screen_heatmap_thumb.png" class="img-responsive"></a></div>
<div class="col-sm-4 col-xs-4"><a href="images/screen_network.png" target="_blank" class="thumbnail"><img src="images/screen_network_thumb.png" class="img-responsive"></a></div>
<div class="col-sm-4 col-xs-4"><a href="images/screen_vulcano.png" target="_blank" class="thumbnail"><img src="images/screen_vulcano_thumb.png" class="img-responsive"></a></div>
</div>


Until publication of the MS Atlas database paper, please kindly cite the following paper when using MS Atlas (data) for your research:
> Tobias Frisch, Maria L. Elkjaer, Richard Reynolds, Tanja Maria Michel, Tim Kacprowski, Mark Burton, Torben A. Kruse, Mads Thomassen, Jan Baumbach, Zsolt Illes<br>
  "MS Atlas - A molecular map of brain lesion stages in progressive multiple sclerosis"<br>
   bioRxiv 584920; doi: https://doi.org/10.1101/584920
  

</div>
