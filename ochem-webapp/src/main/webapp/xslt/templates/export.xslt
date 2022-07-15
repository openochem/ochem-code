<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	<xsl:template name="content">
		<style type="text/css">
			input[type=checkbox] {margin-right: 5px;}
			form img {margin-left: 5px; margin-right: 5px;}
			form a {margin-right: 5px;}
			.params {margin-left: 20px; background-color: #FFC; padding: 10px;}
			#converted-unit-params TD {padding-right: 10px;}
			
			SMALL.restriction {
				padding: 10px;
				border: 1px solid #AAA;
				border-radius: 5px;
				display: block;
				width: 800px;
			}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up silver">
				<h1>Data export</h1>
			Export the selected data as an Excel, CSV or SDF file</td></tr>
			<tr>
				<td class="itunes-right">
					Please, select the items that you want to export:<br/><br/>
					<form method="post">
					
					<a href="javascript:selectAll()">[select all]</a><a class="restriction" href="javascript:selectUnrestricted()">[select unrestricted only]</a><a href="javascript:selectNone()">[select none]</a><br/>
					<xsl:for-each select="exportable-data/availableColumns">
							<input type="checkbox" name="{column}">
								<xsl:if test="checked = 'true' and restricted = 'false'">
									<xsl:attribute name="checked">checked</xsl:attribute>
								</xsl:if>
								<xsl:attribute name="restricted"><xsl:value-of select="restricted"/></xsl:attribute>
							</input>
							<xsl:value-of select="title"/>
							<xsl:if test="restricted = 'true'">
								<img class="restriction" src="img/icons/lock.png" title="Export of this column is limited. You can export only particular amount of such data per week."/>
							</xsl:if>
						<br/>					
					</xsl:for-each>			
					<br/>
					<xsl:if test="exportable-data/properties[@qualitive='false']">
						<div class="params" id="converted-unit-params">
							Select the units to which the exported values will be converted:<br/>
							<table>
							<xsl:for-each select="exportable-data/properties[@qualitive='false']">
								<tr>
									<td><xsl:value-of select="@name"/></td>
									<td>
										<select name="unit-{@id}">
											<xsl:for-each select="unitCategory/unit">
												<option value="{@id}">
													<xsl:if test="@id=../../selectedUnit/@id">
														<xsl:attribute name="selected">selected</xsl:attribute> 
													</xsl:if>
													<xsl:value-of select="@name"/>
												</option>
											</xsl:for-each>
										</select>
									</td>
								</tr>
							</xsl:for-each>
							</table>
						</div>
						<br/><br/>
					</xsl:if>						
					
					<input type="hidden" name="format"/>
					<input type="submit" name="submit" class="fancy-button" value="Get Excel file" format="xls"/>
					<input type="submit" name="submit" class="fancy-button" value="Get CSV file" format="csv"/>
					<input type="submit" name="submit" class="fancy-button" value="Get SDF file" format="sdf"/>
					<input type="submit" name="submit" class="fancy-button" value="Get R script" format="r"/>
					
					<br/><br/>
<!--
					<small class="restriction">
					Note:<br/>
					<img src="img/icons/lock.png"/>The columns marked with a "lock" icon are limited for download. The download of the other (unrestricted) columns is unlimited.<br/>
					<br/><b>This restriction does not apply to freely downloadable data, which are available under the CC-BY 4.0 license.</b>
					<br/>Download of the limited data requires <a href="https://docs.ochem.eu/display/MAN/Bonus+points+system" target="_blank">bonus points.</a> Some amount of free bonus points is restored on a weekly basis. <br/>
					Your currently have <b><xsl:value-of select="exportable-data/freeBonuses"/></b> free weekly bonus points.
					</small> 
					
					<br/>
					<a tab="Export history" href="export/show.do">View export history</a><br/>
					Export history allows to to re-download the previously exported files without any additional restrictions.
-->					
					</form>	
				</td>
			</tr>
		</table>
		<script language="javascript">
			var downloadLimit = <xsl:value-of select="exportable-data/freeBonuses"/>;
			if (downloadLimit &gt;= 1000000)
				$(".restriction").addClass("invisible");
			$("input[type=submit]").click(function(){

<!--			
				if (!$("input[type=checkbox]:checked").length)
				{
					window.alert("Please, select at least one column to export");
					return false;
				}
-->	
				$("input[name=format]").val($(this).attr('format'));
				return true;	
			});
			
			<xsl:for-each select="exportable-data/selectedColumns">
				var selectedColumnName = '<xsl:value-of select="."/>';
				$("[name='" + selectedColumnName + "']").attr("checked", "checked");
			</xsl:for-each>
			
			$("input[name='DESCRIPTORS']").change(function(){
				$("#DESCRIPTORS-params").setClass("invisible", !$(this).is(":checked"));
			}).change();
			
			
			function onConvertedUnitfieldChange()
			{
				invisible = true;
				$("input[name='PREDICTED_VALUE'],input[name='EXP_VALUE_CONVERTED']").each(function(){ 
					if ($(this).is(":checked")) 
					{ 
						invisible = false 
					} 
				});
				$("#converted-unit-params").setClass("invisible", invisible);
			}
			
			$("input[name='PREDICTED_VALUE']").change(onConvertedUnitfieldChange).change();			
			$("input[name='EXP_VALUE_CONVERTED']").change(onConvertedUnitfieldChange).change();
			
			$("input[name='use-shuffle-key']").change(function(){
				$("select[name='shuffle-key']").setClass("invisible", !$(this).is(":checked"));
			}).change();
			
			function selectAll()
			{
				$("input[type=checkbox]").attr("checked", "checked");
				$('input[type=checkbox]').parent(".invisible").removeClass("invisible");
				$(".params>input").removeAttr("checked");
			}
			
			function selectNone()
			{
				$("input[type=checkbox]").removeAttr("checked");
				$('input[type=checkbox]').parent(".params").addClass("invisible");
			}
			
			function selectUnrestricted()
			{
				selectNone();
				$("input[restricted='false']").attr("checked", "checked");
			}
			
			
		</script>
	</xsl:template>
	
</xsl:stylesheet>
