<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<style type="text/css">
			#area INPUT.threshold {width: 70px; text-align: center; background-color: #EFE; border: 1px solid #999;}
			#area INPUT.option {width: 100px; text-align: center;}
			#area {padding: 50px;}
			#area H1 {border-bottom: 1px solid black; font-size: 16pt; font-family: Georgia; font-weight: bold; margin-bottom: 15px;}
			.property {border-top: 1px solid black; padding-top: 10px;}
			TEXTAREA {align: center; height: 400px; font-family: Courier; float: left; margin-right: 10px;}
			.w600 {width: 600px;}
			.hints {float: left;}
			.examples {font-family: Courier; font-size: 9pt; overflow: auto;}
			.hints PRE {font-family: Verdana; font-size: 14pt;}
			.error {color: #800;}
			.success {color: #080;}
		</style>
		
		<form id="oscript-form">
		<div id="area">
			<h1>Basket transformation</h1>
			<table>
			<tr><td>Basket</td><td><a href="basket/edit.do?id={//basket/@id}&amp;render-mode=popup" tab="Basket profile"><xsl:value-of select="//basket/@name"/></a></td></tr>
			<tr><td>Properties in the basket</td><td><xsl:for-each select="//basket/propertyUsed/property"><a href="properties/edit.do?id={@id}" tab="Property profile"><xsl:value-of select="@name"/></a>, </xsl:for-each></td></tr>
			<xsl:if test="//conditions-used/condition">
				<tr><td>Conditions used in the basket&#160;</td><td><xsl:for-each select="//basket/conditions-used/condition"><a href="properties/edit.do?id={@id}" tab="Property profile"><xsl:value-of select="@name"/></a>, </xsl:for-each></td></tr>
			</xsl:if>
			</table>
			<br/>
			New basket name:<br/><input type="text" name="new-basket-name" value="{//basket/@name} (transformed with OScript)" class="w600"/><br/>
			<input type="hidden" name="id" value="{//basket/@id}"/>
			<br/>
			Transformation script: <a href="#" onclick="helpDialog.show(); return false;">[?]</a><br/>
			<textarea name="script" class="w600">// Enter the transformation rules below:
</textarea>	
			<div class="hints">
				<pre id="msg">
				</pre>
</div>
<br style="clear: both;"/><br/>
<input type="button" name="submit" value="Transform the basket" style="clear: both;" onclick="longOperation.start();"/>
		</div>
		
		</form>
		
<div id="helpDialog"> 
	    <div class="hd">OScript examples</div> 
	    <div class="bd"> 
	    	OScript transformation scripts allow you to transform your basket using a set of simple rules.
	    	<br/><br/>For example, you may treat two properties as one in your QSAR model or, vise versa, split one property into two using a condition value (e.g., different properties for different Species).
	    	<br/>N.B.! All properties and conditions should exist before OScript is used.
	    	<br/>Examples of simple transformation rules follow:
	      	<div class="examples">
<pre>
<b>Simple renaming of a property (or a merge of two properties)</b>
logKow-&gt;logPow

<b>Merge two numeric properties, distinguish by a condition</b>
LC50_Rat -> LC50[Species=Rat]
LC50_Fish -> LC50[Species=Fish]

<b>Split a property into two using a condition</b>
LC50[Species=Rat] -> logPow_Rat
LC50[Species=Fish] -> logPow_Fish

<b>Merge condition values</b>
logPow[Species=Rat] -> logPow[species=Rodent]
logPow[Species=Mouse] -> logPow[species=Rodent]

<b>Using intrevals with conditions</b>
"Boiling Point"[Pressure &gt; 100] -> Boiling_Point_Low_Pressure
"Boiling Point"[Pressure &lt; 100] -> Boiling_Point_High_Pressure

// .. or
"Boiling Point"[Pressure &gt; 100]-> Boiling Point[Pressure_type=High]
"Boiling Point"[Pressure &lt; 100]-> Boiling Point[Pressure_type=Low]
</pre>

			</div>  
	    </div> 
	</div>
		
		<script language="javascript">
			var longOperation = new LongOperation({
				url: "baskettransformer/submit.do",
				formScope: "#oscript-form",
				finished: function(status){
					status = status.split("~~~");
					$("#msg").html(status[0] + '   <a href="epbrowser/show.do?basket-select='+status[1]+'" tab="Transformed basket">view the transformed basket</a>').addClass("success").removeClass("error");
					$(document).trigger("DOM_updated");
				},
				error: function(error){
					$("#msg").html(error).addClass("error").removeClass("success");
					longOperation.waitingDialog.hide();
					$(document).trigger("DOM_updated");
				}
			});
			
			var helpDialog = new YAHOO.widget.Dialog("helpDialog", { 
	    		width:"600px", 
				fixedcenter:true, 
				modal:true, 
				visible:false
	    	});
	    	
	    	helpDialog.render();
		</script>
	
	
	
						
		</xsl:template>
</xsl:stylesheet>
