<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:template name="content">		
		<style type="text/css">
			.floating {float: left; margin:10px;}
			.standalone {clear: both;}
			IMG.mol {width: 100px; height: 100px;}
			IMG.big.mol {width: 200px; height: 200px;}
		</style>
		<script language="javascript" type="text/javascript" src="js/lib/excanvas.pack.js"></script>
		<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/commons/plotting.js"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1><doc term="Using+model+browser">Models applier browser</doc></h1>
			</td></tr>
			<tr><td class="itunes-right">
				<h1>Prediction results</h1>
				<a tab="Export predictions" href="modelapplier/exportPredictions.do"><img src="img/icons/xls.gif"/> Export results in a file (Excel, CSV or SDF)</a>
				<br/>
				
				<xsl:for-each select="others/model">
					<xsl:if test="/model/param=./@id">
						<img src="img/icons/chart_dots.png"/><a href="model/profile.do?add-validation-set=1&amp;id={@id}" tab="Model profile">Add the results as a validation set for model <i><xsl:value-of select="@name"/></i></a><br/>
					</xsl:if>
				</xsl:for-each>
				<br/>
				
				<div id="DivExample">
					<!-- Applicability domainn -->
					<xsl:if test="others/ad">
						<a href="javascript:showAD()"><img src="img/icons/ad.gif"/> Show applicability domain &gt;&gt;</a>
						<div id="demo" class="yui-navset"> 
						    <ul class="yui-nav"> 
						    	<xsl:for-each select="others/ad">
							    	<li class="selected">
							    		<a href="#tab{position()}"><em><xsl:value-of select="model/model/@name"/></em></a>
							    	</li>
						    	</xsl:for-each>
						    </ul> 
						    
							<div class="yui-content">
								<xsl:for-each select="others/ad">
									<div style="overflow: auto;">
										Applicability domain for model 
										<a href="model/profile.do?mapping_id={model/@id}" tab="Model profile"><i><xsl:value-of select="model/model/@name"/></i></a> 
										for <xsl:value-of select="model/property/@name"/>  
										using dataset <i><xsl:value-of select="model/model/training-set/@name"/></i><br/>
										Distance to model: <i><xsl:value-of select="@name"/></i><br/>
										<div id="placeholder{position()}" style="width:400px;height:250px;float: left;"></div>
										<div style="float: left; margin-left: 20px; font-size: 70%;">Gray line shows percentage of molecules with DM less than a threshold<br/></div>
									</div>
								</xsl:for-each>
							</div>
						</div>
					</xsl:if>
					<br/>
					Sorting
					<select name="ordering" filter="1">
						<option value="">none</option>
						<option value="prediction">by prediction value</option>
						<xsl:if test="others/ad">
							<option value="dm">by prediction accuracy</option>
						</xsl:if>
					</select><br/>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
					<div id="pager">
					</div>
					<div id="Browser">
					</div>
					<div class="pager-strip">
						<b class="showed">none</b> of <b class="total">none</b>
					</div>
	      			<div class="standalone">	
		      			<div class="formsubmit">		
						<input type="button" name="submit" value="&lt;&lt;Back" onclick="location.href='modelapplier/apply.do';"/>
						</div>
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
			resBrowser.itemElement="datatable";
			resBrowser.itemTemplate = 'js/templates/models/modelresult.ejs';
			
			function zoom(link)
			{
				$(link).find('img').toggleClass("big");
			}
			
			$(document).ready(function() {resBrowser.initialize();});
			
			<xsl:if test="others/ad">
				var dt;
				var plot;
				var plots = new Array();
				<xsl:for-each select="others/ad">
					plot = new Plot();
					dt = ([
					<xsl:for-each select="ms/set[1]/point">
						[<xsl:value-of select="ad/@dmValue"/>, <xsl:value-of select="ad/@error"/>]
						<xsl:if test="position() != last()">,</xsl:if>
					</xsl:for-each>
					]);
					plot.addData(dt, "#C00");
					dt = ([
					<xsl:for-each select="ms/set[@setId = 'predicted']/point">
						[<xsl:value-of select="ad/@dmValue"/>, <xsl:value-of select="ad/@error"/>]
						<xsl:if test="position() != last()">,</xsl:if>
					</xsl:for-each>
					]);
					plot.addData(dt, "#0C0");
					
					// 7 Apr 2001 - This functionality is neither used nor supported anymore. Commented out. / Midnighter
					// Distribution of DM in the set
					//var percentages = new Array();
					//for (var i = 0; i &lt; dt.length; i++)
					//{
					//	percentages.push([dt[i][0], i / dt.length * 100]);
					//}
					//plot.addData(percentages, "#888");
					//plot.currentData.yaxis = 2;
					//plot.setPoints({line: true});
					
					var AD = undefined;
		       		<xsl:apply-templates select="ad-configuration" mode="adData"/>
		       		if (AD)
		       		if (AD.intervals)
		       		{
		       			// Draw MGD-steps
		       			plot.drawLine(plot.minX, AD.errors[0], AD.intervals[0], AD.errors[0]);
		       			plot.currentData.color = "#000";
		       			for (var i = 0; i &lt; AD.intervals.length - 1; i++)
		       			{
		       				plot.drawLine(AD.intervals[i], AD.errors[i+1], AD.intervals[i+1], AD.errors[i+1]);
		       				plot.currentData.color = "#000";	
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