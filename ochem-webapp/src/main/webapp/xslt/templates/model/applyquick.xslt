<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/midnighter.js"></script>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1><doc term="Using+model+browser">Models applier browser</doc></h1>
			</td></tr>
			<tr><td class="itunes-right">
				<div id="process">
					<img src="img/roller.gif" id="progress-img"/><br/>
					<span id="status">Initializing...</span>
				</div>
				
				<form action="modelquick/results.do" method="post">
					<div class="formsubmit">
						<input type="button" name="submit" value="&lt;&lt;Back" onclick="location.href='modelquick/select.do';"/>
						<input type="submit" id="next" name="next" value="Next&gt;&gt;" disabled="disabled"/>
					</div>	
					
				</form>
			</td></tr>
		</table>
		
		<script language="javascript">
			var ajax = new QSPR.Ajax();
			var timer;
			
			function checkStatus()
			{
				ajax.url = 'modelquick/status.do';
				ajax.send({
					data:'',
					success: function(xml)
					{
						var status = xml.message.message;
						
						if (status == "Finished")
						{
							$("[name='next']").removeAttr("disabled");		
							clearInterval(timer);	
							
							// To delete browser history entry
							window.location.replace(webRoot + "modelquick/results.do?render-mode=popup");				
						}
						else 
						{
							if (status.substring(0, 5) == "Error")
							{
								$("#progress-img").attr('src', 'img/icons/error.jpg');		
								clearInterval(timer);	
							}
							$("#status").html(status.replace(/\_\$\$\_/g, "<br/>"));
						}
					},
					error: function(msg)
					{
						$("#status").html('QSPR server is not available');
					}
				});
			}
			
			function start()
			{
				$("#status").html('Starting...');
				ajax.url = 'modelapplier/start.do';
				ajax.send({
					success: function()
					{
						$("#status").html('Requesting status...');
						timer = setInterval("checkStatus()", 3000);
					}
				});
			}
			
			function onResultsSuccess() {
				var stop = "";
			
			}
			
			$(document).ready(function() {start();});
		</script>
	</xsl:template>
</xsl:stylesheet>