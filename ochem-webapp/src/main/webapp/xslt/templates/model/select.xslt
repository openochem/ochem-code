<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/browsers/model-browser.js"></script>
		<script language="javascript">
			include.plugins('view');
			var sampleBrowser = new ModelBrowser("model", "model");
			
			$(document).ready(function() {sampleBrowser.initialize();});
			
			function onNext()
			{
				if ($("img[src='img/icons/checked.gif']").size() == 0)
				{
					window.alert("You must select at least one model to continue");
					return false;
				}
				return true;
			}
		</script>
		<style type="text/css">
			 .command-panel {white-space: nowrap; width: 85px;}
             .compact-item {height: 65px;}
              A.apply-model {background-color: #EEE; padding: 4px; border-radius: 3px; color: #444;white-space: nowrap;margin-right: 10px; border: 1px solid #888;}
              A.apply-model:hover { background-color: #DDD; text-decoration: none;}
              
              .upper-filters NOBR {margin-top: 5px; padding-top: 5px;}

		</style>
		<title>Model applier</title>
		<table width="100%">
			<tr><td class="itunes-up silver">
				<h1><doc term="Using+model+browser">Models applier browser</doc></h1>
				The complete list of models at OCHEM available for you is displayed below. If you are new here, you can also switch to a simplified <a href="predictor/show.do?render-mode=popup">OCHEM predictor</a>
			</td></tr>
			<tr><td class="itunes-right">
				 <a action="selectsubmit" class="fancy-button invisible">Submit selected models</a><br/>
				<div class="upper-filters">
				<nobr>
				Model name or model ID:
				<input type="text" name="query" filter="1"/>
				and property name:
				<input type="text" name="proquery" filter="1"/>
				</nobr>
				<nobr>Models visibility:
				<select name="visibility" filter="1">
					<option value="private">Private</option>
					<option value="public">Public</option>
					<option value="all" selected="selected">Public and private</option>
				<xsl:if test="//user/superuser = 'true'" >
					<option value="publishedAll">All published</option>
					<option value="allModels" selected="selected">All tasks</option>
				</xsl:if>					
				</select>
				Order by:
				 <select name="order" filter="1">
				 	<option value="creation">creation time</option>
					<option value="access">last access time</option>
				 	<option value="modification">last modification time</option>
				 </select> 
				<xsl:if test="/model/session/user/group">
					<input type="checkbox" name="group" filter="1"/><xsl:value-of select="/model/session/user/group/name"/> members' models
				</xsl:if>
				<a href="javascript:sampleBrowser.request(true)" class="fb-button">
					refresh list
				</a>
				</nobr>
				</div>
				<div class="pager-strip">
					<span><b class="showed">none</b> of <b class="total">none</b></span>
					<div id="pager" class="pgr">
					</div>
				</div>
				<div id="aaa">&#160;</div>
				<input type="hidden" name="name" value="" filter="1"/>	
				<div id="Browser">
				</div>
				
				<div class="pager-strip">
					<span><b class="showed">none</b> of <b class="total">none</b></span>
					<div id="pager" class="pgr">
					</div>
				</div>
				
				<form action="modelapplier/apply.do" method="post" onsubmit="return onNext();">
					<div class="formsubmit">
						<input type="submit" id="next" name="next" value="Next&gt;&gt;"/>
					</div>	
				</form>			
			</td></tr>
		</table>
		
		<div id="xmlDialog">
			<div class="hd">XML Configuration of the model</div> 
   		 	<div class="bd">
   		 		<div style="overflow: auto; height: 300px;">
   		 			<pre></pre>
   		 		</div>
			</div>
		</div>
		
		<script language="javascript">
			var xmlDialog = new YAHOO.widget.Dialog("xmlDialog", { 
			    	width:"700px", 
					fixedcenter:true, 
					modal:true, 
					visible:false 
		    });
		    
		    xmlDialog.render();
		</script>
	</xsl:template>
</xsl:stylesheet>