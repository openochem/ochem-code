<?xml version="1.0" encoding="UTF-8"?>
<!-- Universal template to track the status of long-running server side operations -->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<div style="font-size: 150%; text-align: center; margin-top: 100px;">
			<img src="img/roller.gif"/>
			<span id="status"></span>
		</div>
		
		<script language="javascript">
			var operation = new LongOperation({
				tracker: $("#status"),
				finished: function(){
					window.location.replace(webRoot + "longoperations/operationSuccess.do?render-mode=popup&amp;operation-id=" + operation.operationId);
				}
			});
			
			LongOperation.currentOperation = operation;
			operation.operationId = getParams['operation-id'];
			operation.startCheckingStatus();
		</script>
	</xsl:template>
</xsl:stylesheet>
