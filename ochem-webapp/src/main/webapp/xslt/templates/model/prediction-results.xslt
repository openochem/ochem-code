<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:template name="content">		
		<style type="text/css">
			.floating {float: left; margin:10px;}
			.standalone {clear: both;}
			IMG.mol {width: 100px; height: 100px;}
			IMG.big.mol {width: 200px; height: 200px;}
			.cached, .in-ad, .out-ad {padding: 4px; font-size: 75%; color: white;
				-moz-border-radius: 5px;
				border-radius: 5px;
			}
			.cached, .in-ad {background-color: #070 !important; }
			.out-ad {background-color: #700 !important; }
			
			#stats {
				float: right; 
				width: 400px; 
				background-color: #F4F4F4; 
				padding: 10px; 
				font-size: 90%;
				margin-bottom: 10px;
				margin-left: 10px;
			}
			
			.red-error
			{
				fcolor: red;
			}
			
			.error {color: #A00 !important; font-style: normal;}
			
			.top-align 
			{
				vertical-align: top;
			}
			
			.result-piece
			{
				padding: 5px;
			}
			
			.actions {
				display: none;
			}
			
			.block-table TR:hover .actions {display: inline;}
		</style>
		<script language="javascript" type="text/javascript" src="js/lib/excanvas.pack.js"></script>
		<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/commons/plotting.js"></script>
	  	<title>Prediction results</title>
		<table width="100%">
			<tr><td class="itunes-up silver">
				<h1><doc term="Using+OCHEM+predictor">OCHEM predictor - results</doc></h1>
				Here you can browse the predictions for your compounds and export them in a variety of formats
			</td></tr>
			<tr><td class="itunes-right">
				<xsl:if test="//modelApplierSummary/results">
					<div id="stats">
						<b>Predicted accuracy for the set</b><br/>
						<xsl:for-each select="//modelApplierSummary/results">
							<a tab="Model profile" href="model/profile.do?mapping_id={modelMapping/@id}" title="Open model profile"><xsl:value-of select="modelMapping/property/@name"/></a> for <xsl:value-of select="simulatedStatistics/@matchedPoints"/> compounds
							<div style="padding-left: 10px;">
							<xsl:choose>
								<xsl:when test="modelMapping/property/@type = 0">	
									RMSE = <xsl:value-of select="simulatedStatistics/rmse"/> ± <xsl:value-of select="simulatedStatistics/rmsestd"/><br/>
									MAE = <xsl:value-of select="simulatedStatistics/mae"/> ± <xsl:value-of select="simulatedStatistics/maestd"/><br/>
								</xsl:when>
								<xsl:otherwise>
									Accuracy = <xsl:value-of select="simulatedStatistics/classificationSummary/accuracyTotal/@formatted-value"/>%<br/>
								</xsl:otherwise>
							</xsl:choose>
							</div>
						</xsl:for-each>
					</div>
				</xsl:if>
				<xsl:if test="//modelApplierSummary/@thereAreErrors">
						<b>Predictions for some molecules failed:</b><br/>
						<xsl:for-each select="//modelApplierSummary/models">
							<xsl:if test="predictionErrors != 0">
								<xsl:value-of select="@name"/> had <xsl:value-of select="predictionErrors"/> failed predictions<br/>
							</xsl:if>
						</xsl:for-each>		
						<br/>
				</xsl:if>
				
				<xsl:if test="//modelApplierSummary/@conditions">	
					<xsl:value-of select="//modelApplierSummary/@conditions"/><br/>
					<br/>
				</xsl:if>
								
				<a tab="Export predictions" href="modelapplier/exportPredictions.do"><img src="img/icons/xls.gif"/> Export results in a file (Excel, CSV or SDF)</a>
				<div class="mmp-link invisible">
				<a tab="MMP analysis of a predicted set" href="mmpqsar/prediction.do?prediction=1"><img src="img/icons/mmp-16.png"/> Analyse matched molecular pairs</a>
				</div>
				<br/>
				
				<xsl:for-each select="//modelApplierSummary/models">
					<xsl:if test="/model/param=./@id">
						<img src="img/icons/chart_dots.png"/><a href="model/profile.do?add-validation-set=1&amp;id={@id}" tab="Model profile">Add the results as a validation set for model <i><xsl:value-of select="@name"/></i></a><br/>
					</xsl:if>
				</xsl:for-each>
				<br/>
				
				<div id="DivExample">
					<!-- Applicability domainn -->
					<xsl:if test="//modelApplierSummary/ads">
						<a href="javascript:showAD()"><img src="img/icons/ad.gif"/> Advanced applicability domain charts</a>
						<div id="demo" class="yui-navset"> 
						    <ul class="yui-nav"> 
						    	<xsl:for-each select="//modelApplierSummary/ads">
							    	<li class="selected">
							    		<a href="#tab{position()}"><em><xsl:value-of select="model/model/@name"/></em></a>
							    	</li>
						    	</xsl:for-each>
						    </ul> 
						    
							<div class="yui-content">
								<xsl:for-each select="//modelApplierSummary/ads">
									<div style="overflow: auto;">
										Applicability domain for model 
										<a href="model/profile.do?mapping_id={model/@id}" tab="Model profile"><i><xsl:value-of select="model/model/@name"/></i></a> 
										for property <i><xsl:value-of select="model/property/@name"/></i>  
										<br/>
										<br/>
										<table>
											<tr>
												<td>
													<div id="placeholder{position()}" style="width:400px;height:250px;float: left;"></div>
												</td>
												<td valign="top">
													<div style="float: left; margin-left: 30px; font-size: 100%; max-width: 600px;">
													The applicability domain chart allow to estimate the expected prediction accuracy.<br/>
													The green dots indicate the predicted compounds, where its X-position is its "distance to model" and its Y-position is the expected predicton accuracy (for classification models) or the expected RMSE (for regression models).
													</div>
												</td>
											</tr>
											<tr>
												<td align="center">Distance to model (<xsl:value-of select="@name"/>)</td>
												<td></td>
											</tr>
										</table>
										
										
									</div>
								</xsl:for-each>
							</div>
						</div>
						<br/>
					</xsl:if>
					<br/>
					
					<nobr>
					<xsl:choose>
						<xsl:when test="//modelApplierSummary/@thereAreErrors">
							Display : 
							<select name="display-mode" filter="1">
								<option value="all" selected="selected">all predictions</option>
								<option value="errors">only errors</option>
								<xsl:if test="count(//modelApplierSummary/models) &gt; 1">		
									<xsl:for-each select="//modelApplierSummary/models">
										<xsl:if test="predictionErrors != 0">
											<option value="model-{@id}">only errors for model <xsl:value-of select="@name"/></option>
										</xsl:if>
									</xsl:for-each>						
								</xsl:if>
							</select>
						</xsl:when>
						<xsl:otherwise>
							<input type="hidden" name="display-mode" value="all" filter="1"/>
						</xsl:otherwise>
					</xsl:choose>
					Sorting
					<select name="ordering" filter="1">
						<option value="">none</option>
						<option value="prediction">by prediction value</option>
						<xsl:if test="//modelApplierSummary/ads">
							<option value="dm">by prediction accuracy</option>
						</xsl:if>
					</select>
					<input type="checkbox" name="asc" filter="1"/> Ascending
					</nobr>
					<br/>
					<div class="pager-strip">
						<span><b class="showed">none</b> of <b class="total">none</b></span>
						<div id="pager" class="pgr pager">
						</div>
					</div>
					
					<div id="Browser">
					</div>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
				</div>		
			<div id="DivExample2" style='visibility:hidden'>
			<img src='/img/long_green.gif' width='160' height='16'></img>									
			</div>
			</td></tr>
		</table>
		
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var resBrowser = new Browser();
			resBrowser.url = "modelapplier/reslist.do";
			resBrowser.itemElement="predictionRow";
			resBrowser.itemTemplate = 'js/templates/models/modelresult.ejs';
			resBrowser.doExpanderror = function()
			{
				$(this.currentBlock).find(".error").addClass("invisible");
				$(this.currentBlock).find(".fullerror").removeClass("invisible");
			}
			
			resBrowser.listenEvent("items_loaded", function(){
				if (resBrowser.pager.totalNum &gt; 100)
					$(".mmp-link").removeClass("invisible");
			});
			
			resBrowser.doZoom = function(link)
			{
				var size = Math.round(Math.min(window.innerHeight, window.innerWidth) * 0.9);
				var id = this.currentBlock.find('.block-image IMG').attr("id");
				var img = $(this.molViewDialog.body).find("IMG");
				img.attr("src","");
				img.attr("src","depiction.jsp?w="+size+"&amp;h="+size+"&amp;id="+id);
				img.attr("width",size);
				img.attr("height",size);
			    this.molViewDialog.render();
				this.molViewDialog.show();
			}
			
			function zoom(link)
			{
				$(link).find('img').toggleClass("big");
			}
			
			$(document).ready(function() {
				resBrowser.molViewDialog = new YAHOO.widget.Dialog("molViewDialog", { 
					fixedcenter:true, 
					modal:true, 
					visible:false,
					postmethod: "none" 
			    });
				resBrowser.initialize();
			});
			
			models = new Array();
			<xsl:for-each select="//modelApplierSummary/models">
				models.push({id: <xsl:value-of select="@id"/>, trainingSetId: <xsl:value-of select="training-set/@id"/>});
			</xsl:for-each>
			
			<xsl:if test="//modelApplierSummary/ads">
				var dt;
				var plot;
				var plots = new Array();
				<xsl:for-each select="//modelApplierSummary/ads">
					var qualitative = false;
					plot = new Plot();
					
					<xsl:if test="model/property/@qualitive != 'true'">
					dt = ([
					<xsl:for-each select="ms/set[1]/point">
						[<xsl:value-of select="ad/@dmValue"/>, <xsl:value-of select="ad/@error"/>]
						<xsl:if test="position() != last()">,</xsl:if>
					</xsl:for-each>
					]);
					plot.addData(dt, "#C00");
					</xsl:if>
					
					dt = ([
					<xsl:for-each select="ms/set[@setId = 'predicted']/point">
						[<xsl:value-of select="ad/@dmValue"/>, <xsl:value-of select="ad/@error"/>]
						<xsl:if test="position() != last()">,</xsl:if>
					</xsl:for-each>
					]);
					plot.addData(dt, "#0C0");
					
					<xsl:if test="model/property/@qualitive = 'true'">
						plot.minX = 0;
						plot.minY = 0;
						qualitative = true;
					</xsl:if>
					
					
					var AD = undefined;
		       		<xsl:apply-templates select="ad-configuration" mode="adData"/>
		       		if (AD)
		       		if (AD.intervals)
		       		{
		       			plot.minX = 0;
		       			// Draw MGD-steps
		       			plot.drawLine(plot.minX, AD.errors[0], AD.intervals[0], AD.errors[0]);
		       			var color = qualitative ? "#777" : "#000";
		       			plot.currentData.color = color;
		       			for (var i = 0; i &lt; AD.intervals.length - 1; i++)
		       			{
		       				var inc = qualitative ? 0 : 1;
		       				plot.drawLine(AD.intervals[i], AD.errors[i + inc], AD.intervals[i+1], AD.errors[i+1]);
		       				plot.currentData.color = color;	
		       			}
		       		}
		       		plots.push(plot);
					
				</xsl:for-each>
				
				$(document).ready(function()
	   			{
	   				for (var i = 0; i &lt; plots.length; i++)
	   					plots[i].render("#placeholder"+(i+1));
	   				tabView = new YAHOO.widget.TabView("demo");	
	   				$("#demo").addClass("invisible");
	   			});
	   			
	   			function showAD()
	   			{
	   				$("#demo").toggleClass("invisible");
	   			}
			</xsl:if>
		</script>
		
		<div id="molViewDialog"> 
		    <div class="hd">Molecule structure</div> 
		    <div class="bd">
		    	<a href="javascript:void()" width="300" height="300" onclick="compoundBrowser.molViewDialogClose(); return false;"><img src=""/></a> 
		    </div> 
		</div>
	</xsl:template>
	
	<xsl:template match="*" mode="adData">
		var AD = new Object();
		AD.intervals = [
			<xsl:for-each select="interval">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
		AD.errors = [
			<xsl:for-each select="error">
				<xsl:value-of select="."/><xsl:if test="position() != last()">, </xsl:if>
			</xsl:for-each>
			];
	</xsl:template>
</xsl:stylesheet>