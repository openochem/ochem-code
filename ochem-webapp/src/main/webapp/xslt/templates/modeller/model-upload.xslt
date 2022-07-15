<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-sets.xslt" />
	
	<xsl:template name="content">
	<style type="text/css">
		DIV.validation {margin-top: 5px; margin-left: 5px;}
		#template-configuration {padding: 10px; background-color: #FFE;}
	</style>
	<title>Model Builder</title>
		<table width="100%">
			<tr><td class="itunes-up silver">
				<h1><doc term="Uploading+a+stub+QSAR+model">Upload a stub model</doc></h1>
				<p>Select the training and validation sets, and upload the predicted values from an external file</p>
			</td></tr>
			<tr><td class="itunes-right">
				<form name="modelform" id="modelform" action="modelupload/submit.do" method="post" enctype="multipart/form-data">
					<br/>
					<xsl:call-template name="select-sets"/>
					<h1>Provide the uploaded model information</h1>
					<b>Upload a file with predicted values</b> (you can take a look at a sample <a href="/documents/model-upload.xls">here</a>)<br/><input type="file" name="file"/><br/><br/>
					<b>Provide brief model description</b><br/>
					<textarea name="description" onkeyup="return checkModelDescription(this)" style="width: 500px; height: 80px;">
					</textarea>
					<div id="modelDesc" class="message" style="color:red">model description is empty</div>
					<div class="formsubmit"> 
					<input type="button" action="submit" disabled="disabled" value="Upload" name="next"/>
					</div>
					
				</form>
			</td></tr>
		</table>
		<script language="javascript">
			function checkModelDescription(obj)
			{
				if (obj.value.length &lt; 1)
				{
					$("#modelDesc").css("color","red").html("model description is empty");
					$("[name=next]").attr("disabled","disabled");
				}else
				if (obj.value.length &gt;= 1 &amp;&amp; obj.value.length &lt; 100)
				{
					$("#modelDesc").css("color","#FA58F4").html("model description length is "+obj.value.length+" and may not be informative");
					$("[name=next]").removeAttr("disabled");
				}
				else if(obj.value.length &gt;= 100)
				{
					$("#modelDesc").css("color","green").html("model property description is good");
					$("[name=next]").removeAttr("disabled");
				}
			}
		</script>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/blocks/model-template.js?ver=1.8.9"></script>
	</xsl:template>
</xsl:stylesheet>