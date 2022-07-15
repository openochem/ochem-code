<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - Multi-Learning options</title>
		<h1>Property weighting options</h1>
		<style type="text/css">
			.weights TD, .weights TH {padding: 6px; border-right: 3px solid white;}
			.weights TD {background-color: #FFE;}
			.weights TH {background-color: #FFE; font-weight: bold;}
			
			.cost-matrix INPUT {width: 40px; border: 1px solid gray; text-align: right; margin: 3px; font-size: 14pt;}
			
		</style>

		
		<input type="checkbox" name="nomultilearning"/><label>Disable multi-property learning and develop individual model for each property</label>
		<br/>
		<br/>
		(Experimental for sparse data) Substitute implicit values (separated by ;) with this category:<input type="text" name="implicitValues" value=""/>
		
		<br/><br/>
		N.B.! Weighting is under development and <b>(currently available only for DNN, DEEPCHEM and NNF).</b><br></br> It allows to give more weight to a particular property or class while learning the model.
		This could be used for training that contain multiple properties and/or multiple classes (e.g., "active" or "inactive").<br/>
		The default option is to give equal weight to each property/class.
		For classification the property weight is a multiplication of the property and class weights.
		<br/>
		<br/>
		Weighting options :<br/>
		<input type="checkbox" name="global-weighting"/><label>Set weights inversionally proportional to class/property frequencies</label>
		<br/>
		<input type="checkbox" name="global-normalization" checked="checked"/><label>Normalize weights</label>
		<br/>
		<br/>
	
		
		<table class="weights">
		<tr style="border-bottom: 1px solid black;">
			<th>Property/class name</th>
			<th>Weight</th>
		</tr>
		
		<xsl:for-each select="//others/labelWeighting/property-weight">
			<tr style="padding-top: 5px;">
				<td><xsl:value-of select="@name"/></td><td><input type="text" value="{@weight}" name="prop-{@name}"/></td>
			</tr>
			<xsl:for-each select="property-class">
				<tr>	
					<td style="padding-left: 20px;">•   <xsl:value-of select="@name"/></td><td><input type="text" value="{@weight}" name="prop-{../@name}$$$$class-{@name}"/></td>	
				</tr>
			</xsl:for-each>
		</xsl:for-each>
		</table>
		<script language="javascript">
			$("[property]").change(function(){
				$("#cm-" + $(this).attr("property")).setClass("invisible", !$(this).is(":checked"));
			});
			$(document).ready(function(){
				setTimeout("$('[property]').change()", 100);
			});
			
			setCheckbox("global-weighting", '<xsl:value-of select="//attachment/configuration/modelConfiguration/labelWeighting/globalWeighting"/>');
			setCheckbox("global-normalization", '<xsl:value-of select="//attachment/configuration/modelConfiguration/labelWeighting/globalNormalization"/>');
			setCheckbox("nomultilearning", '<xsl:value-of select="//attachment/configuration/modelConfiguration/noMultiLearning"/>');
		</script>
		
		
	</xsl:template>
	
</xsl:stylesheet>