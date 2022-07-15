<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="select-sets">
		<style type="text/css">
			DIV.section {margin-left: 15px;}
		</style>
		<div id="select-sets">
		<div id="hidden-value">
		
			<h1>Select the training and validation sets:</h1><br/>
			<div class="section">
				<div id="training-set">
					Training set <i style='color:gray;'>(required)</i>:
					<a action="selectset" bindto="trainingsetid" storein="tset-title" title="Click to change">[...]</a>
					<span id="trainingsetid-more" class="invisible">   <a action="moreonset" set="trainingsetid">[details]</a></span><br/>
				</div>
				
				<div id="val-sets">
					<xsl:for-each select="//*[position() &lt;= 1]">
						<div class="invisible">
							Validation set #<xsl:value-of select="position()"/>: <a action="selectset" bindto="validationsetid{position()-1}" storein="vset-title{position()-1}" title="Click to change">[...]</a>
							<a action="delete" class="delete-link invisible">[x]</a>
							<span id="validationsetid{position()-1}-more" class="invisible"><a action="moreonset" set="validationsetid{position()-1}">[details]</a></span>
							<input type="hidden" name="validationsetid{position()-1}" class="invisible"/>
							<input type="hidden" name="vset-title{position()-1}" class="invisible"/>	
						</div>
					</xsl:for-each>
				</div>
				<a href="#" onclick="addValidationSet(); return false;">Add a validation set</a>
			</div>
			<br/>
				
			<div id="modelText"></div>
			<div id="propertyBrowser" class="section">
				<xsl:for-each select="//*[position() &lt;= 1]">
					<div id="property{position()-1}"  class="invisible"><a bindto="property{position()-1}-id" storein="property{position()-1}-title" action="propertyclick">[...]</a>
					using unit: <select id="unit{position()-1}" name="unit{position()-1}"></select></div>
				</xsl:for-each> 
			</div>
		
			<input name="upload" type="hidden" value="" />
			<input name="test" type="hidden" value="" />
			<input type="hidden" name="trainingsetid" class="invisible"/>
			
			<input type="hidden" name="tset-title" class="invisible"/>
			
			
			<xsl:for-each select="//*[position() &lt;= 1]">
				<input type="hidden" name="property{position()-1}-id" class="invisible" fromreset="true"/>
				<input type="hidden" name="unitcategory{position()-1}-id" class="invisible" fromreset="true"/>
				<input type="hidden" name="unit{position()-1}-id" class="invisible" fromreset="true"/>
				<input type="hidden" name="property{position()-1}-title" class="invisible" fromreset="true"/>
			</xsl:for-each>
			</div>
		</div>
		
		<script language="javascript">
				// Load the UI from template
				var setSetById = function(id, validationSet)
				{
					var url = 'basket/list.do?id='+id+'&amp;out=json&amp;public=1';
					$.ajax({
						url: url,
						dataType: "json",
						success: function(response)
						{
						console.log(response);
							var inputID = validationSet ? "validationsetid" : "trainingsetid";
							var div = validationSet ? $("#val-sets div.invisible").eq(0) : $("#training-set");
							div.removeClass("invisible");
							form.onSetSelected(response.list.basket, div.find("[bindto]").get(0));
						}
					});	
				}
				
				function getSetsQueryString()
				{
					var params = new Array();
					$("#select-sets").find("input,select").each(function(){
						if ($(this).val() &amp;&amp; $(this).attr("name").indexOf("title") == -1)
							params.push($(this).attr("name") + "=" + $(this).val());
					});
					
					return params.join("&amp;");
				}
				
				function addValidationSet()
				{
					var baskWin = openTab("Select a validation set", webRoot+"basket/show.do?render-mode=popup");
					baskWin.callback = function(basket)
					{
						var div = $("#val-sets div.invisible").eq(0);
						//div.removeClass("invisible");
						form.onSetSelected(basket, div.find("a[action=selectset]").get(0));
						baskWin.closeTab();
					}	
				}
				
				<xsl:if test="//ochem-model/trainingSetId">
					setSetById(<xsl:value-of select="//ochem-model/trainingSetId"/>, false);	
				</xsl:if>
				<xsl:for-each select="//ochem-model/validationSetId">
					setSetById(<xsl:value-of select="."/>, true);	
				</xsl:for-each>
				
				QSPR.defaultPropertyUnit = new Array();
				<xsl:for-each select="//ochem-model/properties/property">
					QSPR.defaultPropertyUnit['<xsl:value-of select="@name"/>'] = '<xsl:value-of select="@unit"/>';
				</xsl:for-each>
			</script> 
	</xsl:template>
	
</xsl:stylesheet>