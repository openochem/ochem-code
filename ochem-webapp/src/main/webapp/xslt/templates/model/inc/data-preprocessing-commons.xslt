<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="data-preprocessing">
		<xsl:param name="ranges"/>
		<b><doc term="Molecule+preprocessing">Preprocessing of molecules</doc></b><br/>
		 
		<doc term="Molecule+preprocessing#Moleculepreprocessing-Standardization" hide="true">
		Choose your preferred standardizer: <select name="standardizer">
			<xsl:for-each select="standardizerOptions/standardizer">
				<xsl:if test="self::node()[text()!='NONE' and text()!='Chemaxon']">
					<option value="{.}"><xsl:value-of select="."/></option>
				</xsl:if>
			</xsl:for-each>
		</select>
		</doc> 
		<doc term="Molecule+preprocessing#Moleculepreprocessing-Standardization" hide="true"><input type="checkbox" checked="checked" name="standardize"/> Standardization</doc>
		<doc term="Molecule+preprocessing#Moleculepreprocessing-Neutralize" hide="true"><input type="checkbox" checked="checked" name="neutralize"/> Neutralize</doc>
		<doc term="Molecule+preprocessing#Moleculepreprocessing-Remove+salts" hide="true"><input type="checkbox" checked="checked" name="desalt" onchange="updateVisibility()"/> Remove salts</doc>
	
		<xsl:if test="//param[@key='use-mixtures-validation'] != 'false'">
			<div id="desalt-section"  style="margin-left: 30px;" class="invisible section">
				Additional validation option for mixtures (salts):<br/> 
				<input type="radio" name="mixtureValidation" value="component" checked="checked"/>Validation by each mixture component<br/>
				<input type="radio" name="mixtureValidation" value="anions"/>Validation by anions (ILs only)<br/>
				<input type="radio" name="mixtureValidation" value="cations"/>Validation by cations (ILs only)<br/>
				<!--
				<input type="radio" name="mixtureValidation" value="all"/>Validation by all components (most rigorous)<br/>
				<input type="radio" name="mixtureValidation" value="maxcomp"/>Validation by maximum mixture component (skip validation of counterions)<br/>
				-->
				<input type="radio" name="mixtureValidation" value="mixture"/>Validation by mixtures (each mixture is considered as a compound - likely to overfit)<br/>
			</div>
		</xsl:if>
		
		<br/>
		<doc term="Molecule+preprocessing#Moleculepreprocessing-Clean+structure" hide="true"><input type="checkbox" checked="checked" name="cleanstructure"/> Clean structure</doc><br/><br/>
		
		<xsl:if test="(param[@key='numIntervals'] != '0') or (param[@key='numAppequals'] != '0') or (param[@key='numGreaterless'] != '0') or $ranges">
			<b>Records with ranges</b><br/>
			
			<div style="margin-left: 20px;" id="ranges">
			Include following records:<br/>
			<xsl:if test="(param[@key='numAppequals'] != '0') or $ranges">
			<input type="checkbox" checked="checked" name="approximateequals"/>Include "approximately equals" records <xsl:if test="(param[@key='numAppequals'] != '0')">(<i><xsl:value-of select="param[@key='numAppequals']"/> records</i>)</xsl:if><br/>		
			</xsl:if>
			<xsl:if test="(param[@key='numIntervals'] != '0') or $ranges">
			<input type="checkbox" checked="checked" name="intervals"/>Include interval records <xsl:if test="(param[@key='numIntervals'] != '0')">(<i><xsl:value-of select="param[@key='numIntervals']"/> records</i>)</xsl:if><br/>
			</xsl:if>
			<xsl:if test="(param[@key='numGreaterless'] != '0') or $ranges">
			<input type="checkbox" checked="checked" name="greaterless"/>Include "greater" and "less" records <xsl:if test="(param[@key='numGreaterless'] != '0')">(<i><xsl:value-of select="param[@key='numGreaterless']"/> records</i>)</xsl:if><br/>
			</xsl:if>
			<div id="handleRanges">
				<br/>
				Handling of records with ranges:<br/> 
				<input type="radio" name="handleRanges" value="false"/>Use average values (for intervals) and boundary values (for greater and less ranges)<br/>
				<input type="radio" name="handleRanges" value="true" checked="checked"/>Handle ranges as ranges (experimental)<br/>
			</div>
			</div>
		</xsl:if>
		
		<script language="javascript">
			(function()
			{
				var updateVisibility;
				$(document).ready(function(){
				
					$("input[type=checkbox]").change(updateVisibility = function(){
						$("#handleRanges").setClass("invisible", $("#ranges input[type=checkbox]:checked").length == 0);
						$("#desalt-section").setClass("invisible", $("input[name=desalt]").is(":checked"));
					});
					
					updateVisibility();
				});
			}());
			
			
			<xsl:if test="//ochem-model">
				$("[name=standardize]").setChecked(<xsl:value-of select="count(//ochem-model/attachment/standartization/standardizeWith) &gt; 0"/>);
				$("[name=neutralize]").setChecked(<xsl:value-of select="count(//ochem-model/attachment/standartization/neutralizeWith) &gt; 0"/>);
				$("[name=desalt]").setChecked(<xsl:value-of select="count(//ochem-model/attachment/standartization/desaltWith) &gt; 0"/>);
				$("[name=cleanstructure]").setChecked(<xsl:value-of select="count(//ochem-model/attachment/standartization/cleanStructureWith) &gt; 0"/>);
				
				setCheckbox("approximateequals", '<xsl:value-of select="//ochem-model/attachment/datahandling/approximateequals = 'use'"/>');
				setCheckbox("intervals", '<xsl:value-of select="//ochem-model/attachment/datahandling/intervals = 'use'"/>');
				setCheckbox("greaterless", '<xsl:value-of select="//ochem-model/attachment/datahandling/greaterless = 'use'"/>');
				$("[name='mixtureValidation']").val(<xsl:value-of select="//ochem-model/attachment/protocol/baggingConfiguration/mixtureValidation"/>);
				
				<xsl:if test="//ochem-model/attachment/datahandling/handleRanges = 'false'">
					setRadio("handleRanges", "false");
				</xsl:if>
			</xsl:if>
		</script>
	</xsl:template>
	
</xsl:stylesheet>