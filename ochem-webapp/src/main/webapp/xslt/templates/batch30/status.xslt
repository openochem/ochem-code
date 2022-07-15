<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<title>Uploaded data browser</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/blocks/batch-upload.js"></script>
		<style type="text/css">
			#process 
			{
				align: center; text-align: center; width: 100%;
				font-family: Tahoma;
				font-size: 14pt;
				position: relative;
				top: 200px;
			}
		</style>
		<div id="process">
		<img src="img/roller.gif"/><br/>
		<span id="status">Initializing...</span><br/>
		<a href="#" onclick="setInterrupt(); return false;">[interrupt]</a>
		</div>
		
		<script language="javascript">
			var ajax = new QSPR.Ajax();
			var interval = 100;
			var timer = setInterval("checkStatus()", interval);
			var interrupt = 0;
			
			function setInterrupt()
			{
				interrupt = 1;
			}
						
			function checkStatus()
			{
				clearInterval(timer);
				interval *= 2;
				
				if (interval > 5000)
					interval = 5000;
					
				timer = setInterval("checkStatus()", interval);
				
				ajax.url = 'batchupload30/status.do';
				ajax.send({
					data: (interrupt == 1) ? 'nodb=1&amp;interrupt=1' : 'nodb=1',
					success: function(json)
					{
						var message = json.message.message;
						var url = json.message.title;
						if (message != "null")
							$("#status").html(message);
						
						if (url != '') 
						if (url != 'status')
						{
							location.replace(webRoot + "batchupload30/"+url+".do?render-mode=popup");
						}	
					},
					error: function(msg)
					{
						$("#status").html('QSPR server is not available');
					}
				});
			}
		</script>
	</xsl:template>
</xsl:stylesheet>