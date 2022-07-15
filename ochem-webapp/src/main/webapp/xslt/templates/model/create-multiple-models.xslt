<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-sets.xslt" />
	<xsl:include href="../model/inc/data-preprocessing-commons.xslt" />
	<xsl:include href="../model/inc/structure-optimisation-commons.xslt" />
	

	<xsl:template name="content">
		<style type="text/css">
			.templates TD {background-color: #EEF; border: 10px solid white; padding: 20px; vertical-align: top;}
			.templates TD INPUT {margin-right: 3px;}
			.templates TD H1 {font-family: Georgia; font-size: 12pt;}
			.templates TD A {font-size: 90%; margin-bottom: 5px;}
			#success-message {background-color: #FFE; padding: 10px;}
		</style>
		<title>Multiple Models Builder</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/browser.js?ver=1.8.6"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/blocks/model-template.js?ver=1.8.9"></script>
		<script language="javascript" src="js/blocks/select-model.js"></script>	
		
		<table width="100%">
			<tr><td class="itunes-up">
				<img src="img/comprehensive.png" id="page-icon"/>
				<h1 id="page-title">Comprehensive modeling</h1>
				Create multiple models simultaneously
			</td></tr>
			<tr><td class="itunes-right">
				<div id="start-models">
				The comprehensive modeling feature allows you to simultaneously run multiple models with different machine learning methods, molecular descriptors and validation protocols.<br/>
				Please note that running multiple models may require significant computational resources and time.
				<br/><br/>
				<xsl:call-template name="select-sets"/>
				<br/>
				
				<b>Select the methods you want to use for the modeling:</b>
				<table class="templates creator-scope">
					<tr>
						<td>
							<h1>Method</h1>
							<a href="#" onclick="selectAll(this, true); return false;">[all]</a><a href="#" onclick="selectAll(this, false); return false;">[none]</a><br/>
							
							<xsl:for-each select="//others/model-template[type='METHOD']">
								<div>
									<span><input type="checkbox" name="template" value="{@id}">
										<xsl:if test="recommended = 'true'">
											<xsl:attribute name="checked">true</xsl:attribute>
										</xsl:if>
									</input><xsl:value-of select="@name"/></span>
									<xsl:if test="/model/session/user/rank = 10 or /model/session/user/@id = introducer/@id">
										<a title="Edit this template" tab="Edit model template" href="modeltemplate/edit.do?id={@id}"> [edit]</a>
										<a title="Delete this template" action="delete" ajax-data="id={@id}" class="delete-link"> [x]</a>
									</xsl:if>
									<br/>
									<div class="invisible description">
										<xsl:value-of select="description"/>
									</div>
								</div>
							</xsl:for-each>
							<xsl:call-template name="add-custom">
								<xsl:with-param name="type">METHOD</xsl:with-param>
							</xsl:call-template>
						</td>
						<td>
							<h1>Descriptors</h1>
							<a href="#" onclick="selectAll(this, true); return false;">[all]</a><a href="#" onclick="selectAll(this, false); return false;">[none]</a><br/>
							<xsl:for-each select="//others/model-template[type='DESCRIPTORS']">
								<xsl:if test="@id != '573' or //session/user/@id = '5'">
								<div>
									<span><input type="checkbox" name="template" value="{@id}">
										<xsl:if test="recommended = 'true'">
											<xsl:attribute name="checked">true</xsl:attribute>
										</xsl:if>
									</input><xsl:value-of select="@name"/></span>
									<xsl:if test="/model/session/user/rank = 10 or /model/session/user/@id = introducer/@id">
										<a title="Edit this template" tab="Edit model template" href="modeltemplate/edit.do?id={@id}"> [edit]</a>
										<a title="Delete this template" action="delete" ajax-data="id={@id}" class="delete-link"> [x]</a>
									</xsl:if>
									<br/>
									<div class="invisible description">
										<xsl:value-of select="description"/>
									</div>
								</div>
								</xsl:if>
							</xsl:for-each>
							<xsl:call-template name="add-custom">
								<xsl:with-param name="type">DESCRIPTORS</xsl:with-param>
							</xsl:call-template>
						</td>
						<td>
							<h1>Descriptor selection</h1>
							<a href="#" onclick="selectAll(this, true); return false;">[all]</a><a href="#" onclick="selectAll(this, false); return false;">[none]</a><br/>
							<xsl:for-each select="//others/model-template[type='SELECTION']">
								<div>
									<span><input type="checkbox" name="template" value="{@id}">
										<xsl:if test="recommended = 'true'">
											<xsl:attribute name="checked">true</xsl:attribute>
										</xsl:if>
									</input><xsl:value-of select="@name"/></span>
									<xsl:if test="/model/session/user/rank = 10 or /model/session/user/@id = introducer/@id">
										<a title="Edit this template" tab="Edit model template" href="modeltemplate/edit.do?id={@id}"> [edit]</a>
										<a title="Delete this template" action="delete" ajax-data="id={@id}" class="delete-link"> [x]</a>
									</xsl:if>
									<br/>
									<div class="invisible description">
										<xsl:value-of select="description"/>
									</div>
								</div>
							</xsl:for-each>
							<xsl:call-template name="add-custom">
								<xsl:with-param name="type">SELECTION</xsl:with-param>
							</xsl:call-template>
						</td>
						<td>
							<h1>Model validation</h1>
							<a href="#" onclick="selectAll(this, true); return false;">[all]</a><a href="#" onclick="selectAll(this, false); return false;">[none]</a><br/>
							<xsl:for-each select="//others/model-template[type='VALIDATION']">
								<div>
									<span><input type="checkbox" name="template" value="{@id}">
										<xsl:if test="recommended = 'true'">
											<xsl:attribute name="checked">true</xsl:attribute>
										</xsl:if>
									</input><xsl:value-of select="@name"/></span>
									<xsl:if test="/model/session/user/rank = 10 or /model/session/user/@id = introducer/@id">
										<a title="Edit this template" tab="Edit model template" href="modeltemplate/edit.do?id={@id}"> [edit]</a>
										<a title="Delete this template" action="delete" ajax-data="id={@id}" class="delete-link"> [x]</a>
									</xsl:if>
									<br/>
									<div class="invisible description">
										<xsl:value-of select="description"/>
									</div>
								</div>
							</xsl:for-each>
							<xsl:call-template name="add-custom">
								<xsl:with-param name="type">VALIDATION</xsl:with-param>
							</xsl:call-template>
						</td>
					</tr>
				</table>
				<small><a href="#" onclick='toggleVisibility(this, "advanced"); return false;' alternate-text="Hide advanced options">Show advanced options&gt;&gt;</a><br/></small>
				<div id="advanced" class="invisible params">
				<input type="checkbox" name="skipDublicates" checked="checked"/> Skip duplicates (if a model with the same configuration already exists)<br/><br/>
					<xsl:call-template name="data-preprocessing">
						<xsl:with-param name="ranges">true</xsl:with-param>
					</xsl:call-template><br/><br/>

				<a href="#" onclick='toggleVisibility(this, "mixtures"); return false;' alternate-text="Show mixture options">Show mixture options (required if "Remove salts" is not checked)&gt;&gt;</a><br/>
					<div id="mixtures" class="mixtures invisible">
						<b>Mixture processing options:</b>
						<select name="mixtures">
				<option value="fraction">Descriptors calculated for the whole mixture plus molar fractions are added as descriptors</option>
				<option value="concatenate">Descriptors are concatenated (only binary mixtures)</option>			
				<option value="wsum"  checked="checked">Descriptors are averaged (weighted by molar fraction; any mixtures)</option>			
				<option value="wsumdiff">Sum and absolute differences of weighted by molar fraction descriptors (binary mixtures)</option>
						</select>
				<table style="margin-left: 30px;">
				<tr><td><b>Additional validation option for mixtures (salts):</b></td></tr>
				<tr>
				<td>
				<input type="radio" name="mixtureValidation" value="component" checked="checked"/>Validation by each mixture component<br/>
				<input type="radio" name="mixtureValidation" value="anions"/>Validation by anions (ILs only)<br/>
				<input type="radio" name="mixtureValidation" value="cations"/>Validation by cations (ILs only)<br/>
				<!--
				<input type="radio" name="mixtureValidation" value="all"/>Validation by all components (most rigorous)<br/>
				<input type="radio" name="mixtureValidation" value="maxcomp"/>Validation by maximum mixture component (skip validation of counterions)<br/>
				-->
				<input type="radio" name="mixtureValidation" value="mixture"/>Validation by mixtures (each mixture is considered as a compound - likely to overfit)<br/>
				</td>
				</tr>
				</table>
				<br></br>
				</div>
				<br></br>

				<a href="#" onclick="showConditions(); return false;">Show conditions&gt;&gt;</a><br/>
				<div id="conditions" class="invisible">
					<br/>Conditions of experiments<br/>
					<div class="params" id="Conditions-Browser">
					</div>
				 </div>	

				<xsl:call-template name="structure-optimisation">
					<xsl:with-param name="ranges">true</xsl:with-param>
				</xsl:call-template><br/><br/>

					<input type="checkbox" onclick="toggleVisibility(this, 'normalisation'); return true;" name="override-normalisation"/> Override default normalisation options
					<div class="invisible params" id="normalisation">
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
					<br class="cb"/>
					<input type="checkbox" name="save-models" checked="checked"/> <span title="Uncheck this option to skip saving the models. This will make your models unapplicable to new compounds.">Save models</span>
					<br/>
					
					<br class="cb"/>
					<input type="checkbox" name="singleLearningEnforce"/> <span title="Check to apply multi-task learning as single-task">Enforce Single Task Learning (STL) for MTL method</span>
				<br/>
					<input type="checkbox" name="featureNet"/> <span title="Check to use feature net approach (see DOI:10.1021/ci8002914)">Use Feature Net learning for multitask basket</span>
				<br/>
					<input type="checkbox" name="sanitize"/> <span title="Use sanitisation of SMILES for RDKit (no descriptors methods)">Use sanitisation</span>
				<br/>
					(Experimental) Merge with the uploaded descriptor type: <input type="text" name="additionalDescriptorsMerge" value=""/>
				<br/>
					(Experimental for sparse data) Substitute implicit values (separated by ;) with: <input type="text" name="implicitValues" value=""/>
				<br/>
					(Experimental for taxonomy) Use taxonomy: <input type="text" name="taxonomyValues" value=""/>
				<br/>
				<br/>
				</div>
				<br/>
				
				Considering the selection above, <b><span id="models-count"></span> models</b> are requested to be created. Notice that only maximum 25 models based on descriptors can be started simultaneously.<br/><br/>
				<input id="start-button" type="button" class="fancy-button" value="Create the models" onclick="start()"/>
				</div>
				<div class="finished invisible creator-scope">
					<h1>Success!</h1>
					<div id="success-message"></div><br/><br/>
					<p>
						<!-- You can track the status of your models in the <a action="basket_summary">basket summary page</a>. -->
						You can track the status of your models in the <input type="button" class="fancy-button" value="basket summary page" onclick="basketSummary()"/>
					</p>
				</div>
				
			</td></tr>
		</table>
		
		<div id="xmlDialog">
			<div class="hd">XML Configuration of the model</div> 
   		 	<div class="bd">
   		 		<div style="overflow: auto; height: 300px;">
					<pre id="configuration-xml">
					</pre>
				</div>
			</div>
		</div>
		
		<script language="javascript">
			var ajax = new QSPR.Ajax();
			var timer = 0;
			
			var creator = new Actionable();
			creator.actionURL = "modeltemplate/action.do";
			creator.scope = ".creator-scope";
			form.scope = "#select-sets";
		    
		    creator.xmlDialog = new YAHOO.widget.Dialog("xmlDialog", { 
			    	width:"700px", 
					fixedcenter:true, 
					modal:true, 
					visible:false 
			    });
			    
		    creator.xmlDialog.setContent = function(content)
		    {
		    	$("#configuration-xml").html(content);
		    }
		    
		    creator.doAddtemplate = function(link)
		    {
		    	var wnd = openTab("Select a model to create template from", webRoot+"model/select.do?single=1");
				wnd.callback = function(entity)
				{
					var ajax = new QSPR.Ajax();
					ajax.send({
						url: webRoot + "/multiplemodels/addCustomTemplate.do?type="+$(link).attr("type")+"&amp;model=" + entity.id,
						success: function()
						{
							window.alert("A new template has been successfully added");
							window.location.reload();
						}
					});
					wnd.closeTab();
				}
		    }
		    
		    creator.doShowxml = function()
		    {
		    	window.alert("xml!");
		    }
		    
		    creator.doBasket_summary = function()
		    {
		    	openTab("Basket summary", "multiplemodels/show.do?set=" + form.getValue("trainingsetid"))
		    }
		    
		    function basketSummary()
		    {
		    	openTab("Basket summary", "multiplemodels/show.do?set=" + form.getValue("trainingsetid"))
		    }
		    
		    
		    creator.beforeDelete = function(link)
		    {
		    	creator.deleteLink = link;
		    	return window.confirm("Are you sure you want to delete this template?");
		    }
		    
		    creator.onDeleteSuccess = function()
		    {
		    	$(creator.deleteLink).parent("div").remove();	
		    }
		    
		    creator.xmlDialog.render();
			
			function updateModelsCount()
			{
				var total = 1;
				$(".templates TD").each(function(){
					total *= $(this).find("input:checked").length;
				});
				$("#models-count").html(total);
				creator.totalModels = total;
				updateVisibility();
				
			}
			
			function updateVisibility()
			{
				$("#start-button").setDisabled(creator.totalModels == 0);
			}
			
			function getCheckedIds()
			{
				var ids = new Array();
				$(".templates INPUT:checked").each(function(){
					ids.push($(this).attr("value"));
				});
				return ids.join(",");
			}
			
			function getQueryString(scope)
			{
				var data = new Array();
				$(scope).find("input").each(function(){
					var elem = $(this);
					if (elem.attr("type") == "text")
						data.push(elem.attr("name") + "=" + elem.attr("value"));
					else if (elem.attr("type") == "checkbox")
					{
						if (elem.is(":checked")){
							if(elem.attr("name") == "condition")
								data.push(elem.attr("name") + "=" + elem.val());
									else
								data.push(elem.attr("name") + "=checked");
							}
					}
					else if (elem.attr("type") == "radio") 
					{
						if (elem.is(":checked"))
							data.push(elem.attr("name") + "=" + elem.val());
					}
				});
				$(scope).find("select").each(function(){
					var elem = $(this);
					data.push(elem.attr("name") + "=" + elem.val());
				});
				
				return data.join("&amp;");
			}
			
			function start()
			{
				if (!form.getValue("trainingsetid"))
					window.alert("Please select a training set!");
				else 
					new LongOperation({
						url: "multiplemodels/createSubmit.do",
						data: "out=json&amp;ids=" + getCheckedIds() + "&amp;" + getSetsQueryString() + "&amp;" + getQueryString("#advanced"),
						finished: function(status)
						{
							status = status.split(":")[1];
							$("#start-models").addClass("invisible");
							$(".finished").removeClass("invisible");
							$("#success-message").html(status);
						}
					}).start();
			}
			
			function toggleVisibility(link, id)
			{
				var ele = $("#" + id);
				if (ele.hasClass("invisible"))
					ele.removeClass("invisible");
				else
					ele.addClass("invisible");
				var alternate = $(link).attr("alternate-text");
				$(link).attr("alternate-text", $(link).html());
				$(link).html(alternate);
				
			}
			
			function selectAll(link, all)
			{
				$(link).parent().find("input").setChecked(all);
				updateModelsCount();
			}
		
			$(".templates INPUT").change(function(){
				updateModelsCount();
			});
			
			
			function showConditions()
			{
			var condBrowser = null;
			condBrowser = new Browser();
			condBrowser.url = "basket/listConditions.do";
			if(form.getValue("trainingsetid") == 0)return;
			condBrowser.filters.setValue("id", form.getValue("trainingsetid"));

			condBrowser.container = "Conditions-Browser";
			condBrowser.scope = "#Conditions-Browser";
			condBrowser.itemElement = "property";
			condBrowser.itemTemplate = "js/templates/models/conscriptor.ejs";
			
//			condBrowser.listenEvent("items_loaded", function(){
//				if (condBrowser.pager.totalNum &gt; 0){
//					$("#Conditions").removeClass("invisible");
//					}
//			});
		
			condBrowser.onItemDrawn = function()
			{	
				this.currentBlock.find(".regression").setClass("invisible", this.currentEntity.type != 0);
				this.currentBlock.find(".classification").setClass("invisible", this.currentEntity.type != 1);
				var dom = this.currentBlock.get(0);
				
				if (this.currentEntity.type == 1)
				{
					var optionsSelect = new DynamicSelect('option-' + this.currentEntity.id, 'properties/listoptions.do', 'option', this.currentBlock);
					optionsSelect.update("id=" + this.currentEntity.id + "&amp;basket="+form.getValue("trainingsetid"));
				}
				else
				{
					var unitsSelect = new UnitSelect('unit-' + this.currentEntity.id, 'unit/list.do', 'unit', this.currentBlock);
					unitsSelect.update("category=" + this.currentEntity.unitCategory.id, this.currentEntity.defaultUnit.id);
				}
				
				
				this.currentBlock.find("input[name=condition]").change(function(){
					$(this).siblings(".details").setClass("invisible", !$(this).is(":checked"));
				});
								
			}
			
			condBrowser.doDetails = function()
			{
					this.currentBlock.find(".mappings").removeClass("invisible");
			}
			
			condBrowser.doAddmapping = function()
			{
				var div = $("<div> = </div>");
				var selectRule = "id=" + this.currentEntity.id + "&amp;basket="+form.getValue("trainingsetid");
				div.prepend(sel1 = $("<select></select>"));
				div.append(sel2 = $("<select></select>"));
				sel1.attr("name", "mapping-" + this.currentEntity.id + "-1");
				sel2.attr("name", "mapping-" + this.currentEntity.id + "-2");
				this.currentBlock.find(".mappings").append(div);
				var dSel1 = new DynamicSelect(sel1, 'properties/listoptions.do', 'option', this.currentBlock);
				dSel1.update(selectRule);
				var dSel2 = new DynamicSelect(sel2, 'properties/listoptions.do', 'option', this.currentBlock);
				dSel2.update(selectRule);
			}
			
			condBrowser.listenEvent("items_loaded", function(){
				condBrowser.mainDiv.find("input[name=condition]").change();
			});
			
			condBrowser.initialize();

			toggleVisibility(this, "conditions")
			}
			
			
			$(document).ready(function(){
				$(".description").each(function(){
					$(this).parent("div").find("span").attr("title", $(this).html().replace(/\n/g, "<br/>"));
				});
				$("span[title]").tooltip({showURL: false});
								
			});
			
			
			updateModelsCount();
		</script>
	</xsl:template>
	
	<xsl:template name="add-custom">
		<xsl:param name="type" />
		<br/>
		<a action="addTemplate" type="{$type}">+add a custom template</a>
	</xsl:template>
	
		
</xsl:stylesheet>