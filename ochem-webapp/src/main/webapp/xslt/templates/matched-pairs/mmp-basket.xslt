<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="pair-zoom.xslt" />
	<xsl:include href="mmp-helpers.xslt" />
	
	<xsl:template name="content">
	
		<xsl:call-template name="transformation-stats-scripts"/>
	
		<title>MMP analysis of a basket</title>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<h1><img src="img/icons/mmp-48.png"/><span class="setcompare">MatchedPairs</span>: MMP analysis of a dataset (experimental) <a class="infolink" target="_blank" href="https://docs.ochem.eu/display/MAN/Molecular+Matched+Pairs" title="Click to read more about molecular matched pairs (aka MMPs) and their use in OCHEM"></a></h1>
				</td></tr>
			<tr>
				<td class="itunes-right">
					The info is based on dataset <b><xsl:value-of select="/model/basket/@name"/></b> for property <b><xsl:value-of select="/model/others/property/@name"/></b>
					<xsl:if test="//others/mmpStatus/unindexedMolecules != 0">
						(excluding <xsl:value-of select="//others/mmpStatus/unindexedMolecules"/> molecules not indexed for MMP analysis)
					</xsl:if>
					<div>
						Minimal pair similarity:
						<select name="similarity" filter="1" onchange="transformationStats.changeSimilarity(); return false;">
							<option value="0">Any</option>
							<option value="25">25</option>
							<option value="50" selected="1">50</option>
							<option value="75">75</option>
						</select> 
						<a class="infolink" title="Minimal Tanimoto similarity between two molecules in a pair. Used to hide matched pairs that with too dissimilar molecules."></a>
					</div>
					<br/><br/>
					<div class="yaxis">Transformation effect (&#916;<sub>pair</sub>)</div>
					<table class="plot-container" width="97%">
						<tr>
							<td height="300">
								<div id="mmp-plot">
								</div>
							</td>
							<td rowspan="3" width="100%" align="left" style="text-align: left; padding-left: 15px;">
								<xsl:call-template name="transformation-stats-layout"/>		
							</td>
						</tr>
						<tr>
							<td>Property value</td>
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
