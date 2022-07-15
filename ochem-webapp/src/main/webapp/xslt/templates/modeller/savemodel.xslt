<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

	<xsl:include href="configurators/configurators-common.xslt" />

	<xsl:template name="configuration-content">
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" type="text/javascript" src="js/lib/excanvas.pack.js"></script>
	  	<script language="javascript" type="text/javascript" src="js/lib/jquery.flot-0.6.min.js"></script>
	  	<script language="javascript" src="js/commons/plotting.js"></script>
	  	<script language="javascript">
	  		$("#wizardform").submit(function(){
	  			if (QSPR.discard)
	  				return true;
	  			if ($("input[name=modelname]").val() == "")
	  			{
	  				window.alert("Please enter the name of the model");
	  				return false;
	  			}
	  			return true;
	  		});
	  	</script>
		<h1>Save the model</h1>
		Please enter your model's name: <input type="text" send="1" name="modelname" value="{model/@name}"/><br/><br/>
		<iframe frameborder="0" src="model/profile.do?id={model/@id}&amp;save=1&amp;render-mode=popup" width="100%" height="550px"></iframe>
	</xsl:template>
	
</xsl:stylesheet>