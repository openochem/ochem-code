<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - J48 Decision Tree</title>
		<h1>Configure J48 Decision Tree method</h1>
		<table class="configuration">
			<tr>
				<td>Use unpruned tree</td>
				<td><input type="checkbox" name="use-unpruned"/></td>
			</tr>
			<tr>
				<td>Confidence threshold for pruning [0.001 - 0.999]:</td>
				<td><input type="text" class="small" name="confidence" value="0.25"/></td>
			</tr>
			<tr>
				<td>Minimum number of instances per leaf:</td>
				<td><input type="text" class="small" name="instances-per-leaf" value="2"/></td>
			</tr>			
			<tr>
				<td>Use reduced error pruning:</td>
				<td><input type="checkbox" name="reduced-pruning"/></td>
			</tr>
			<tr>
				<td>Number of folds for reduced error pruning:</td>
				<td><input type="text" class="small" name="num-reduced-pruning-folds" value="3"/></td>
			</tr>
			<tr>
				<td>Use binary splits for nominal attributes:</td>
				<td><input type="checkbox" name="use-binary-splits"/></td>
			</tr>
			<tr>
				<td>Don't perform subtree raising:</td>
				<td><input type="checkbox" name="dont-perform-raising"/></td>
			</tr>
			<tr>
				<td>Don't clean up after the tree has been built:</td>
				<td><input type="checkbox" name="no-cleanup"/></td>
			</tr>
			<tr>
				<td>Use Laplace smoothing for predicted probabilites:</td>
				<td><input type="checkbox" name="use-laplase-smoothing"/></td>
			</tr>			
			<tr>
				<td>Random seed:</td>
				<td><input type="text" class="small" name="seed" value="1005"/></td>
			</tr>			
		</table>
		
		<script language="javascript">
			<xsl:if test="//method = 'WEKA-J48'">
				setCheckbox("use-unpruned", '<xsl:value-of select="//attachment/configuration/modelConfiguration/useUnpruned"/>');
				setValue("confidence", '<xsl:value-of select="//attachment/configuration/modelConfiguration/confidence"/>');
				setValue("instances-per-leaf", '<xsl:value-of select="//attachment/configuration/modelConfiguration/instancesPerLeaf"/>');
				setCheckbox("reduced-pruning", '<xsl:value-of select="//attachment/configuration/modelConfiguration/reducedPruning"/>');
				setValue("num-reduced-pruning-folds", '<xsl:value-of select="//attachment/configuration/modelConfiguration/numReducedPruningFolds"/>');
				setCheckbox("use-binary-splits", '<xsl:value-of select="//attachment/configuration/modelConfiguration/useBinarySplits"/>');
				setCheckbox("dont-perform-raising", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dontPerformRaising"/>');
				setCheckbox("no-cleanup", '<xsl:value-of select="//attachment/configuration/modelConfiguration/noCleanup"/>');
				setCheckbox("use-laplase-smoothing", '<xsl:value-of select="//attachment/configuration/modelConfiguration/useLaplaseSmoothing"/>');
				setValue("seed", '<xsl:value-of select="//attachment/configuration/modelConfiguration/seed"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>