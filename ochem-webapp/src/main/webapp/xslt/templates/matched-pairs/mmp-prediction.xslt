<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="pair-zoom.xslt" />
	<xsl:include href="mmp-helpers.xslt" />
	
	<xsl:template name="content">
	
		<xsl:call-template name="transformation-stats-scripts"/>
	
		<title>MMP analysis of a predicted dataset</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><img src="img/icons/mmp-48.png"/><span class="setcompare">MatchedPairs</span>: MMP analysis of a predicted dataset</h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					The analysis below is based on predictions of <b><xsl:value-of select="//modelMapping/property/@name"/></b> for set <b><xsl:value-of select="//param[@key='dataset']" /></b>
					<xsl:if test="//others/mmpStatus/unindexedMolecules != 0">
						(excluding <xsl:value-of select="//others/mmpStatus/unindexedMolecules"/> molecules not indexed for MMP analysis)
					</xsl:if>
					<div>
						<br/>
						Minimal pair similarity:
						<select name="similarity" filter="1" onchange="transformationStats.changeSimilarity(); return false;">
							<option value="0">Any</option>
							<option value="25">25</option>
							<option value="50" selected="1">50</option>
							<option value="60">60</option>
							<option value="70">70</option>
							<option value="75">75</option>
						</select> 
						<a class="infolink" title="Minimal Tanimoto similarity between two molecules in a pair. Used to hide matched pairs that with too dissimilar molecules."></a>
						<br/>
						<div class="transformations-browser">
							<input type="checkbox" filter="1" name="bootstrapStatistics"/> Take the estimated prediction accuracy into account <a class="infolink" help="ad-help"></a>
							<div id="ad-help" class="invisible">
								As you might know, each OCHEM prediction is complemented with the individual accuracy estimate.<br/><br/>
								In the prediction-based MMP analysis, each prediction is randomly "perturbated" based on the estimated prediction accuracy. The statistics is bootstrapped over multiple replicas of the perturbated predictions.
								This allows to identify only the transformations that affect the property in a consistent and reliable manner.
							</div>
						</div>
					</div>
					<br/>
					
					<table class="plot-container" width="97%">
						<tr>
							<td height="300" class="invisible">
								<div class="yaxis">Predicted transformation effect (&#916;<sub>pair</sub>)</div>
								<span style="position: absolute" class="roller"><img src="img/roller-big.gif"/><br/>Loading...<br/></span>
								<div id="mmp-plot">
								</div>
							</td>
							<td rowspan="3" width="100%" align="left" style="text-align: left; padding-left: 15px;">
								<xsl:call-template name="transformation-stats-layout"/>		
							</td>
						</tr>
						<tr>
							<td>Prediction value</td>
						</tr>
						<tr>
							<td height="100%">&#160;</td>
						</tr>
					</table>
					
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			transformationStats.chartURL = "getValueDeltaChart";
			transformationStats.getPointArray = function(p) {
				return [p[3][0], 1*p[3][1] - 1*p[3][0], p[0], p[1], p[2]];
			}
		</script>
		
		<xsl:call-template name="pair-zoom"/>
		
	</xsl:template>
</xsl:stylesheet>
