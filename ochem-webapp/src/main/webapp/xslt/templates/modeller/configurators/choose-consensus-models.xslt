<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<style type="text/css">
			TABLE.torefactor TD, TABLE.torefactor TH {background-color: #FFD; padding: 10px; text-align: left; border-bottom: 1px solid white;}
			TABLE.torefactor TH {font-style: bold; border-bottom: 1px solid black;}
			#checkboxes {margin-top: 20px;}
			#checkboxes INPUT {margin-right: 5px;}
		</style>
		<title>Model builder - Consensus model</title>
		<h1>Choose the individual models for consensus</h1>
		In order to build a consensus model, you must select several (at least two) individual models based on the selected training set <b><xsl:value-of select="/model/model/training-set/@name"/></b>.
		<br/><br/>
		
		<a href="javascript:addModel()">[Add a model]</a><br/>
		<table id="models" class="invisible torefactor">
			<tr><th>Model name</th><th style="text-align: center">Method</th><th></th></tr>
			<tr class="template invisible">
				<td></td>
				<td style="text-align: center"></td>
				<td>
					<input type="hidden" value="" name="model-id"/>
					<a href="#" onclick="removeRow(this); return false;" class="delete-link">[x]</a>
				</td>
			</tr>
		</table>
			
		<br/><br/>Consensus type: <select name="type">
			<option value="OPTIMAL">Optimal combination of models</option>
			<option value="AVERAGE">Simple average</option>
			<option value="RMSE_WEIGHTED">Weighted by model accuracy</option>
			<option value="WEIGHTED_AVERAGE">Weighted by molecule prediction accuracy</option>
			<option value="BEST_MODEL">Best model for molecule</option>
		</select>
		
		<br><br></br></br><doc term="Consensus+model" hide="true"><input type="checkbox" checked="checked" name="allow-errors"/> Ignore errors in individual submodels</doc>
		
		<script language="javascript">
		
			function addModel()
			{
				var win = openTab("Select a model", "model/select.do?trainingSet=<xsl:value-of select="/model/model/training-set/@id"/>");
				win.callback = function(models)
				{
					var lmodels = array(models);
					for (var i = 0; i &lt; lmodels.length; i++)
						drawModel(lmodels[i]);
					win.closeTab();
				}	
			}
			
			function drawModel(model)
			{
				var row = $("#models TR.template").clone();
				row.removeClass("invisible");
				row.removeClass("template");
				row.find("td").eq(0).html(model.name);
				row.find("td").eq(1).html(model.template.name);
				row.find("input").val(model.id);
				$("#models").append(row);
				updateVisibility();
			}
			
			function addModelById(id)
			{
				$.ajax({
					url: 'model/list.do?query='+id+'&amp;out=json',
					dataType: "json",
					success: function(response)
					{
						drawModel(response.list.model);
					}
				});	
			}
			
			function updateVisibility()
			{
				var numOfModels = $("#models tr").length - 2;
				$("#models").setClass("invisible", numOfModels &lt; 1);
				//$("#checkboxes").setClass("invisible", numOfModels &lt; 2);
				if (numOfModels &lt; 2)
					$("input[name=next]").attr("disabled", 1);
				else
					$("input[name=next]").removeAttr("disabled");
			}
			
			function removeRow(link)
			{
				$(link).parents("tr").eq(0).remove();
				updateVisibility();
			}
			
			$(document).ready(function(){
				updateVisibility();
				$("input[type=radio]").change(function(){
					updateVisibility();
				});
				
				<xsl:if test="/model/others/ochem-model/attachment/configuration/individual-models">
				    setValue("type", '<xsl:value-of select="//attachment/configuration/modelConfiguration/type"/>');
					<xsl:for-each select="/model/others/ochem-model/attachment/configuration/individual-models/individual-model">
						addModelById(<xsl:value-of select="@id"/>);
					</xsl:for-each>
				</xsl:if>
			});
			
			// Uncomment when tested
			//$("form").submit(function(){
			//	frameLoading();
			//});
		</script>
		
	</xsl:template>
</xsl:stylesheet>