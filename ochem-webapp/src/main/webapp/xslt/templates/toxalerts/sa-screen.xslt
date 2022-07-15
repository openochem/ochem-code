<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-compounds.xslt" />
	<xsl:include href="../model/inc/data-preprocessing-commons.xslt" />
	<xsl:template name="content">
	
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<style type="text/css">
			.step {border: 1px solid #CCC; width: 600px; padding: 10px; margin-bottom: 5px;}
			.step TABLE TD {padding-bottom: 15px;}
		</style>
		<table height="100%" width="100%">
			<tr><td class="itunes-up">
				<img src="img/icons/screen.png"/>
				<h1><span class="toxalerts">ToxAlerts</span>: Compounds screening</h1>
				
				</td></tr>
			<tr>
				<td class="itunes-right">
				
				<xsl:if test="//message">
					<div class="warning">
						<xsl:value-of select="//message/message"/>
					</div>
				</xsl:if>
				
				
				<form method="post" enctype="multipart/form-data" action="alerts/screenSubmit.do">
					<div class="step">
						<b>Provide the compounds to screen against the structural alerts</b><br/><br/>
						<xsl:call-template name="select-compounds"/>
					</div>
					<div class="step">
						<xsl:call-template name="data-preprocessing"/>
						<!-- Aromatization: <select name="aromatize">
						<option value="none">None</option>
						<option value="Chemaxon-General">Chemaxon General</option>
						<option value="Chemaxon-Basic" selected="selected">Chemaxon Basic</option>
						</select-->
					</div>
					<div class="step">
					<b>Select the structural alerts</b><br/><br/>
					<table>
						<tr><td>Publication&#160;</td>
							<td>
								<select name="article" filter="1">
									<option value="">All articles</option>
									<xsl:for-each select="//available-alerts/articles">
										<option value="{@id}"><xsl:value-of select="publication-date/year"/>&#160;<xsl:value-of select="authors/author[1]/LastName"/></option>
									</xsl:for-each>
								</select>
							</td></tr>
						<tr><td>Endpoint&#160;</td>
							<td><select name="property" filter="1">
									<option value="">All endpoints</option>
									<xsl:for-each select="//available-alerts/endpoints">
										<option value="{@id}"><xsl:value-of select="@name"/></option>
									</xsl:for-each>
								</select>
							</td></tr>
						<tr class="invisible" id="alert-filter"><td>A single alert selected:&#160;</td>
							<td>
								<b><span id="alert-name"></span></b>
								<input type="hidden" name="alert-id"/>
							</td>
						</tr>
						<tr><td colspan="2"><input type="checkbox" name="approved-only"/>&#160;Only approved alerts</td></tr>
						<tr id="selected-filter"><td colspan="2"><input type="checkbox" name="selected-only"/>&#160;Only <b><xsl:value-of select="//available-alerts/selectionSize"/> selected alerts</b></td>
						
						</tr>
					</table>
					</div>
					<input type="hidden" name="article" send="1"/>
					<input type="hidden" name="article-title"/>
					<input type="hidden" name="property" send="1"/>
					<input type="hidden" name="property-title"/>
					
					<input type="submit" name="submit" class="fancy-button" value="Start screening"/><br/><br/>
				</form>
				</td>
			</tr>
		</table>
		
		<script language="javascript">
			var selectionSize = <xsl:value-of select="//available-alerts/selectionSize"/>;
			$("#selected-filter").setClass("invisible", selectionSize == 0);
			if (getParams["alert-id"])
			{
				$("#alert-filter").removeClass("invisible");
				$("[name='alert-id']").val(getParams["alert-id"]);
				$.ajax({
					url: webRoot + "alerts/list.do",
					data: "out=json&amp;id=" + getParams["alert-id"],
					method: "post", 
					dataType: "json",
					success: function(response){
						$("#alert-name").html(response.list["substructure-alert"].name);
					},
					error: function(){
						window.alert("error");
					}
					
				});
			}
				
		</script>
		
		</xsl:template>
</xsl:stylesheet>
