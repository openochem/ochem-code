<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
		<script language="javascript" src="js/commons/plotting.js"></script>
		
		<style type="text/css">
			#area INPUT.threshold {width: 70px; text-align: center; background-color: #EFE; border: 1px solid #999;}
			#area INPUT.option {width: 100px; text-align: center;}
			#area {padding: 50px;}
			#area H1 {border-bottom: 1px solid black; font-size: 16pt; font-family: Georgia; font-weight: bold; margin-bottom: 15px;}
			.property {border-top: 1px solid black; padding-top: 10px;}
		</style>
		<div id="area">
		<h1>Basket discretization</h1>
		
		
		<span style="position: relative" class="threshold invisible">
				&lt;
				<input type="text" name="option" class="option" value="low"/> 
				&lt;
				<input type="text" name="threshold" value="" class="threshold"/>
		</span>
		
		<form action="basket/discretizeSubmit.do" method="get" id="discretize-form">
		<input type="hidden" name="id" value="{//basket/@id}"/>
		You are going to discretize the basket 
		<a tab="Basket profile" href="basket/edit.do?id={//basket/@id}"><xsl:value-of select="//basket/@name"/></a>. 
		The new basket name will be <input type="text" name="basket-name" value="{//basket/@name} (discretized)" size="40"/><br/>
		Please, provide the thresholds for discretization.<br/><br/>
		
		
		<xsl:for-each select="//basket/properties">
			<div class="property">
			<input type="hidden" name="property" value="{@id}"/>
			<b><xsl:value-of select="@name"/></b> (in <xsl:value-of select="defaultUnit/@name"/>) converted into <input type="text" name="property-name" value="{@name} (qualitative)" size="40"/><br/> 
			
			
			<div id="histogram" style="width: 600px; height: 300px;"></div>
			<xsl:apply-templates select="distribution"/>
			
			<br/>
			Thresholds: <input type="text" disabled="1" value="{distribution/min}" class="threshold"/> 
			<span class="thresholds">
			</span>
			
			&lt;
			<input type="text" name="option" class="option" value="high"/> &lt;
			<input type="text" disabled="1" value="{distribution/max}" class="threshold"/>
			<br/>
			<small>To add a threshold, click on the histogram. To remove a threshold, empty its value in the text box.</small>
			</div>
		</xsl:for-each>
		<br/><br/>
		<input type="button" value="Create the discretized basket" onclick="longOperation.start(); return false;"/>
		</form>
		</div>
		
		<script language="javascript">
			var longOperation = new LongOperation({
				url: "basket/discretizeSubmit.do",
				formScope: "#discretize-form",
				finished: function(){
					window.alert("The operation has successfully finished.");
				}
			});
		</script>
		
	</xsl:template>
	
	<xsl:template match="distribution">
		<script language="javascript">
			
			var plot = new Plot();
			plot.options.grid.autoHighlight = false;
			var bins = new Array();
			var frequencies = new Array();
			bins.push(<xsl:value-of select="min"/>);
			<xsl:for-each select="bins">
				bins.push(<xsl:value-of select="."/>);
			</xsl:for-each>
			<xsl:for-each select="frequencies">
				frequencies.push(<xsl:value-of select="."/>);
			</xsl:for-each>
			
			var dt = [];
			for (var i = 0; i &lt; bins.length; i++)
				dt.push([bins[i], frequencies[i]]);
				
			plot.addData(dt);
			plot.currentData.points.show = false;
			plot.currentData.bars = {show: true, barWidth: bins[1] - bins[0]};
			plot.defaultColor = "#000";
			
			plot.thresholds = [1.0 * <xsl:value-of select="middle"/>];
			
			function sortNumber(a,b)
			{
				return a - b;
			}
			
			plot.onClick = function(e, pos)
			{
				plot.thresholds.push(1.0 * pos.x);
				plot.drawThresholds();
			};
			
			plot.drawThresholds = function()
			{
				var container = $(".thresholds");
				
				// Save the labels
				var labels = [];
				container.find(".option").each(function(){
					labels.push($(this).val())
				});
				if (labels.length == 0)
					labels = ["low"];
			
				plot.thresholds = plot.thresholds.sort(sortNumber);
				
				container.html("");
				for (var key in plot.thresholds)
				{
					var span = $(".threshold.invisible").clone().removeClass("invisible");
					span.attr("pos", key);
					span.appendTo(container);
					span.find(".threshold").val(plot.thresholds[key].toPrecision(4));
				}
				
				var i = 0;
				container.find(".option").each(function(){
					$(this).val(labels[i]);
					if (i &lt; labels.length - 1)
						i++;
				});
			
				// Draw vertical lines
				while (plot.dataArray.length > 1)
					plot.dataArray.pop();
				for (var key in plot.thresholds)
					plot.verticalLine(plot.thresholds[key]);
					
				// Update the thresholds as the user types
				container.find("INPUT.threshold").change(function(){
					var pos = 1 * $(this).parent().attr("pos");
					if (!$(this).val())
						plot.thresholds.splice(pos, 1);
					else
						plot.thresholds[pos] = 1.0 * $(this).val();
					plot.drawThresholds();
				});
				
				plot.render("#histogram");
			}
			$(document).ready(function(){
				plot.drawThresholds();
			});
			
		</script>
	</xsl:template>
	
	
</xsl:stylesheet>