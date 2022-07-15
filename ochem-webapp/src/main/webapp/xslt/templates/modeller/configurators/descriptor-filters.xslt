<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - descriptors prefiltering</title>
		<style type="text/css">
			.float {float: left; margin-top: 0px !important;}
			.params {margin-left: 20px; background-color: #FFC; padding: 10px;}
		</style>		
		<h1 class="float">Select filters of descriptors</h1><br/><br/>

		<input type="checkbox" name="unique-values" checked="true"/>
			<label>Eliminate descriptors with less than <input type="text" name="unique-values-threshold" class="small" value="2"/> unique values</label><br/><br/>

		<input type="checkbox" name="maximum" checked="true"/>
			<label>Delete descriptors that have absolute values larger than</label>
			<input type="text" name="maximum-threshold" class="small" value="999999"/><br/><br/>

		<input type="checkbox" name="std" checked="true"/>
			<label>Delete descriptors that have variance smaller than</label>
			<input type="text" name="std-threshold" class="small" value="0.01"/><br/><br/>

		<input type="checkbox" name="decorrelation" checked="true"/>
			<label>Group descriptors, that have pair-wise correlations Pearson's correlation coefficient <i>R</i> larger than</label>
			<input type="text" name="decorrelation-threshold" class="small" value="0.95"/><br/><br/>
			
		<input type="checkbox" name="ufs"/> 
			<label>Use Unsupervised Forward Selection to delete variables using the above value of multiple correlation coefficient <i>R</i></label><br/><br/>
		
						
		<input type="checkbox" name="manual"/><label>After filtering, I want to select necessary descriptors myself  <font color="red">(advanced)</font></label>
		<div style="margin-left: 30px; margin-top: 5px;" id="file" class="invisible">
			Optionally, you can provide a text file with the list of descriptors:<input type="file" name="file-with-descriptor-list"/>
			<br/><small>The file must contain one descriptor name per line; descriptor names are case-sensitive</small>
		</div>
		
		<br/><br/>
		<xsl:if test="not(//param[@key='no-normalisation'])">
		<b>Normalisation parameters</b>
		<div class="params normalisation">
		<table>
			<tr><td>Descriptors normalization</td><td>
				<select name="normx">
					<option value="" selected="true">Do not normalize</option>
					<option value="STANDARDIZE">Standardize (zero mean and unit variance)</option>
					<option value="RANGE">To range [0, 1] (e.g., for LibSVM)</option>	
					<option value="RANGE_MINUS1_PLUS1">To range [-1, 1] (e.g., for LibSVM)</option>
				</select>
			</td></tr>
			<tr><td>Values normalization</td><td>
				<select name="normy">
					<option value="" checked="true">Do not normalize</option>
					<option value="STANDARDIZE">Standardize (zero mean and unit variance)</option>
					<option value="RANGE">To range [0,1] (e.g., for LibSVM)</option>
				</select>
				</td></tr>
		</table>
		</div>
		</xsl:if>
		
		
		<script language="javascript">
			$(document).ready(function(){
				$("[name=manual]").change(function(){
					$("#file").setClass("invisible", !$(this).is(":checked"));
				});
				$("[name=manual]").change();
				
				
				var updateVisibility = function(){
					var details = $("#"+$(this).attr('name')+"-params");
				if (details.length > 0)
					if ($(this).is(':checked'))
						details.removeClass("invisible");
					else
						details.addClass("invisible");
						
				};
				$("input[details]").click(updateVisibility);
				$("input[details]").each(function() {updateVisibility.call(this);});
				
			});
		</script>
		
		<script language="javascript">
			<xsl:if test="//attachment/configuration/selection">
				setCheckbox("unique-values", '<xsl:value-of select="//attachment/configuration/selection/numDifferentValues &gt; 0"/>');
				setValue("unique-values-threshold", '<xsl:value-of select="//attachment/configuration/selection/numDifferentValues"/>');
				setCheckbox("decorrelation", '<xsl:value-of select="//attachment/configuration/selection/correlationThreshold &gt; 0"/>');
				setValue("decorrelation-threshold", '<xsl:value-of select="//attachment/configuration/cselection/orrelationThreshold"/>');
				setCheckbox("maximum", '<xsl:value-of select="//attachment/configuration/selection/maximumValueThreshold &lt; 999999999"/>');
				setValue("maximum-threshold", '<xsl:value-of select="//attachment/configuration/selection/maximumValueThreshold"/>');
				
				
				setCheckbox("maximum", '<xsl:value-of select="//attachment/configuration/selection/maximumValueThreshold &lt; 999999999"/>');
				setValue("maximum-threshold", '<xsl:value-of select="//attachment/configuration/selection/maximumValueThreshold"/>');
				
				setCheckbox("std", '<xsl:value-of select="//attachment/configuration/selection/stdThreshold &gt; 0"/>');
				setValue("std-threshold", '<xsl:value-of select="//attachment/configuration/selection/stdThreshold"/>');
				
				setCheckbox("ufs", '<xsl:value-of select="//attachment/configuration/selection/useUFS"/>');
				
				setCheckbox("pca", '<xsl:value-of select="//attachment/configuration/selection/pca/maximumComponents"/>');
				setValue("std-threshold", '<xsl:value-of select="//attachment/configuration/selection/pca/minimumVarianceThreshold"/>');
				setValue("max-pca-components", '<xsl:value-of select="//attachment/configuration/selection/pca/maximumComponents"/>');
				
				setValue("normx", '<xsl:value-of select="//attachment/configuration/modelConfiguration/scaleTypeX"/>');
				setValue("normy", '<xsl:value-of select="//attachment/configuration/modelConfiguration/scaleTypeY"/>');
			</xsl:if>
			
		</script>
		
	</xsl:template>
	
</xsl:stylesheet>