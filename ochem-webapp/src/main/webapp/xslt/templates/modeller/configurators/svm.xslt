<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<style type="text/css">
			.left{padding-left:15px;}
			.left input{padding-right:5px;}
			.font{font-size: 10.5pt; padding-top:3px;}
			TD {padding-top:2px;}
		</style>
		<script language="javascript">
			function showOptions(block){
				$("#"+block).removeClass("invisible");
			}
			function hideOptions(block){
				$("#"+block).addClass("invisible");
			}
		</script>
		<title>Model builder - SVM</title>
		<h1>Configure SVM method</h1>
		<div class="font">
			<table>
				<tr><td><b>svm scale options</b></td></tr>
				<tr><td class="left"><label>scaling lower limit:</label><input type="text" class="small" name="lower" value="-1"/></td></tr>
				<tr><td class="left"><label>scaling upper limit:</label><input type="text" class="small" name="upper" value="1"/></td></tr>
				<tr><td></td></tr>
			</table>
		</div>
		<div class="font">
			<table>
				<tr><td><b>svm train options</b></td></tr>
				<tr>
					<td  class="left"><label>svm type:</label>
					<select name="svm-type">
						<option value="3">epsilon-SVR</option>
						<option value="4">nu-SVR</option>
					</select></td>
				</tr>
				<tr>
					<td  class="left"><label>kernel type:</label>
					<select name="kernel-type">
						<option value="0">linear</option>
						<option value="1">polynomial</option>
						<option value="2" SELECTED="SELECTED">radial basis function</option>
						<option value="3">sigmoid</option>
					</select></td>
				</tr>
			</table>
		</div>
		<div class="font">
			<b><input type="checkbox" name="grid" checked="true" value="true"/> Grid search options</b> 
				<a href="javascript:showOptions('gridOptions')">[show options]</a>
				<a href="javascript:hideOptions('gridOptions')">[hide options]</a>
			<table id="gridOptions" class="invisible">
				<tr>
				<td  class="left">
					<label> c lower</label><input type="text" class="small" name="c-lower" value="-5"/>
					<label> c upper</label><input type="text" class="small" name="c-upper" value="15"/>
					<label> c increment</label><input type="text" class="small" name="c-increment" value="2"/>
				</td>
				</tr>
				<tr>
					<td class="left">
						<label> g lower</label><input type="text" class="small" name="g-lower" value="-15"/>
						<label> g upper</label><input type="text" class="small" name="g-upper" value="3"/>
						<label> g increment</label><input type="text" class="small" name="g-increment" value="2"/>
					</td>
				</tr>
				<tr>
					<td class="left">
						<label> p lower</label><input type="text" class="small" name="p-lower" value="0.0001"/>
						<label> p upper</label><input type="text" class="small" name="p-upper" value="10"/>
						<label> p multiplication</label><input type="text" class="small" name="p-increment" value="15"/>
					</td>
				</tr>
				<tr><td class="left"><label>n fold</label><input type="text" class="small" name="n-fold" value="5"/></td></tr>
			</table>
		</div>
		<div class="font"><b>Advance options</b><a href="javascript:showOptions('adOptions')">[show options]</a>
			<a href="javascript:hideOptions('adOptions')">[hide options]</a>
			<div id="adOptions" class="invisible left">
				<label>-d <i>(degree)</i>:</label><input type="text" class="small" name="degree" value="3"/>
				<label>-g <i>(gamma)</i>:</label><input type="text" class="small" name="gamma" value="1"/>
				<label>-r <i>(coef0)</i>:</label><input type="text" class="small" name="coef0" value="0"/><br/>
				<label>-c <i>(cost)</i>:</label><input type="text" class="small" name="cost" value="1"/>
				<label>-n <i>(nu)</i>:</label><input type="text" class="small" name="nu" value="0.5"/>
				<label>-p <i>(epsilion)</i>:</label><input type="text" class="small" name="p-epsilion" value="0.1"/><br/>
				<label>-e <i>(eps)</i>:</label><input type="text" class="small" name="eps" value="0.001"/>
				<label>-h <i>(shrinking)</i>:</label><input type="text" class="small" name="shrinking" value="1"/>
			</div>
		</div>
	</xsl:template>
</xsl:stylesheet>