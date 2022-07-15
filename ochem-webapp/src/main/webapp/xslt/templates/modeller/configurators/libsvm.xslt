<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<style type="text/css">
			.left{padding-left:15px;}
			.left input{padding-right:5px;}
			.font{font-size: 10.5pt; padding-top:3px;}
			TD {padding-top:2px;}
			
			#advanced-section TD {text-align: right; padding-left: 3px;}
			.section {
				background-color: #FFE;
				padding: 20px;
			}
		</style>
		<script language="javascript">
			function updateVisibility()
			{
				$("#advanced-section").setClass("invisible", !$("input[name=advanced]").is(":checked"));
				$("#grid-section").setClass("invisible", !$("input[name=grid]").is(":checked"));
				$("#one-class-class").setClass("invisible", $("select[name=use-nu]").val() != "2");
				$("#cwr1").setClass("invisible", !$("input[name=weight]").is(":checked"));
				$("#cwr2").setClass("invisible", !$("input[name=weight]").is(":checked"));
			}
		</script>
		<title>Model builder - SVM</title>
		<h1>Configure LibSVM method</h1>
		<b>SVM Method options</b><br/>
		<div class="section">
		<table>
			<tr>
				<td>SVM algorithm </td>
				<td><select name="use-nu" onchange="updateVisibility()">
					<option value="0">Classic (C-SVC, epsilon-SVR)</option>
					<option value="1">New (nu-SVC, nu-SVR)</option>
					<option value="2">One-class SVM (classification only)</option>
				</select></td>
			</tr>
			<tr id="one-class-class" class="invisible">
				<td>Select class </td>
				<td>
				<select name="one-class-class">
					<xsl:for-each select="others/option">
						<option value="{@id}"><xsl:value-of select="@name"/></option>
					</xsl:for-each>
				</select></td>
			</tr>			
			<tr>
				<td>Kernel type </td>
				<td><select name="kernel-type">
			<option value="0">Linear</option>
			<option value="1">Polynomial</option>
			<option value="2" selected="selected">Radial Basis Function</option>
			<option value="3">Sigmoid</option>
		</select></td>
			</tr>
		</table>
		</div>
		
		<br/>
		<b><input type="checkbox" name="weight" value="true" onchange="updateVisibility()"/> Use class weighting for classification tasks</b><br/>
		<b><input type="checkbox" name="grid" checked="true" value="true" onchange="updateVisibility()"/> Enable grid search</b><br/> 
		<div id="grid-section" class="section">
			<i>Parameters are assigned values 2^(current_step) where current_step = (min, min+step, ..., max)</i><br/>
			<table>
				<tr><td>Cost min, max, step:</td><td><input type="text" class="small" name="cost-min" value="-10"/>,<input type="text" class="small" name="cost-max" value="10"/>,<input type="text" class="small" name="cost-step" value="2"/></td></tr>
				<tr><td>Gamma min,max,step:</td><td><input type="text" class="small" name="gamma-min" value="-10"/>,<input type="text" class="small" name="gamma-max" value="10"/>,<input type="text" class="small" name="gamma-step" value="2"/></td></tr>
				<tr><td>Epsilon min,max,step:</td><td><input type="text" class="small" name="svr-epsilon-min" value="-16"/>,<input type="text" class="small" name="svr-epsilon-max" value="10"/>,<input type="text" class="small" name="svr-epsilon-step" value="2"/></td></tr>
				<tr id="cwr1"><td>Class weight min,max,step:<br/><small>Relevant for classification tasks only</small></td><td><input type="text" class="small" name="class-ratio-min" value="0.1"/>,<input type="text" class="small" name="class-ratio-max" value="1"/>,<input type="text" class="small" name="class-ratio-step" value="0.1"/></td></tr>
				<tr><td>Grid search set size<br/> (fraction of training set):</td><td><input type="text" class="small" name="grid-search-set-size" value="0.1"/></td></tr>
				<tr><td>PARALLEL:</td><td><input type="text" class="small" name="grid-search-parallel" value=""/></td></tr>
			</table>
		</div><br/>
		
		<b><input type="checkbox" name="advanced" onchange="updateVisibility()"/> Configure advanced options</b><br/>
		<div id="advanced-section" class="invisible section">
			<table>
				<tr id="cwr2">
					<td>Class weight ratio:</td><td><input type="text" class="small" name="class-ratio" value="1.0"/></td>
				</tr>
				<tr>
					<td>-d <i>(degree)</i>:</td><td><input type="text" class="small" name="degree" value="3"/></td>
					<td></td><td></td>
				</tr>
				<tr>
					<td>-g <i>(gamma, leave blank for 1/num_features)</i>:</td><td><input type="text" class="small" name="gamma" value=""/></td>
					<td>-c <i>(cost)</i></td><td><input type="text" class="small" name="cost" value="1"/></td>
					<td>-r <i>(coef0)</i>:</td><td><input type="text" class="small" name="coef0" value="0"/></td>
				</tr>
				<tr>
					<td>-n <i>(nu)</i>:</td><td><input type="text" class="small" name="nu" value="0.001"/></td>
					<td>-p <i>(epsilion)</i>:</td><td><input type="text" class="small" name="svr-epsilon" value="0.001"/></td>
					<td>-e <i>(eps)</i>:</td><td><input type="text" class="small" name="epsilon" value="0.001"/></td>
				</tr>
				
			</table>
		</div>
		
		<script language="javascript">
			<xsl:if test="//method = 'LIBSVM'">
				setValue("use-nu", '<xsl:value-of select="//attachment/configuration/modelConfiguration/type"/>');
				setValue("kernel-type", '<xsl:value-of select="//attachment/configuration/modelConfiguration/kernel_type"/>');
				setValue("cost-min", '<xsl:value-of select="//attachment/configuration/modelConfiguration/costMin"/>');
				setValue("cost-max", '<xsl:value-of select="//attachment/configuration/modelConfiguration/costMax"/>');
				setValue("cost-step", '<xsl:value-of select="//attachment/configuration/modelConfiguration/costStep"/>');
				setValue("gamma-min", '<xsl:value-of select="//attachment/configuration/modelConfiguration/gammaMin"/>');
				setValue("gamma-max", '<xsl:value-of select="//attachment/configuration/modelConfiguration/gammaMax"/>');
				setValue("gamma-step", '<xsl:value-of select="//attachment/configuration/modelConfiguration/gammaStep"/>');
				setValue("svr-epsilon-min", '<xsl:value-of select="//attachment/configuration/modelConfiguration/svrEpsilonMin"/>');
				setValue("svr-epsilon-max", '<xsl:value-of select="//attachment/configuration/modelConfiguration/svrEpsilonMax"/>');
				setValue("svr-epsilon-step", '<xsl:value-of select="//attachment/configuration/modelConfiguration/svrEpsilonStep"/>');
				setValue("grid-search-set-size", '<xsl:value-of select="//attachment/configuration/modelConfiguration/gridSearchSetSize"/>');
				setValue("degree", '<xsl:value-of select="//attachment/configuration/modelConfiguration/degree"/>');
				setValue("gamma", '<xsl:value-of select="//attachment/configuration/modelConfiguration/gamma"/>');
				setValue("cost", '<xsl:value-of select="//attachment/configuration/modelConfiguration/cost"/>');
				setValue("coef0", '<xsl:value-of select="//attachment/configuration/modelConfiguration/coef0"/>');
				setValue("nu", '<xsl:value-of select="//attachment/configuration/modelConfiguration/nu"/>');
				setValue("epsilion", '<xsl:value-of select="//attachment/configuration/modelConfiguration/epsilion"/>');
				setValue("eps", '<xsl:value-of select="//attachment/configuration/modelConfiguration/eps"/>');
				setValue("grid-search-parallel", '<xsl:value-of select="//attachment/configuration/modelConfiguration/PARALLEL"/>');
			</xsl:if>
			$(document).ready(function(){
				updateVisibility();	
			});
		</script>
	</xsl:template>
</xsl:stylesheet>