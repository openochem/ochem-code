<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript">
			var form = new AjaxForm();
			var self = form;
			
			
			form.doSubmit = function()
			{
				var nm = self.getValue("moldata");
				if (nm == "") 
					return window.alert("Please provide valid set for the model");
		
				jQuery('#progress').removeClass('invisible');	
				
				start();
				document.modeldata.submit();
				
			}
			
			form.doEditmolecule = function()
			{
				
				jQuery('#depiction').addClass('invisible');
				jQuery('#jme').removeClass('invisible');
			}
			
			form.doSelect = function()
			{
				var jmeMol2 = document.JME.molFile();
				jQuery('#depiction').attr('alt', "whereIsTheMolecule");
				
				var mol = URLEncode(jmeMol2);
				jQuery('#depiction').attr('src', 'depiction.jsp?mol=' + mol);
				
				self.setValue("moldata", jmeMol2);
				self.doCancel();
				
				
			}
			
			form.doCancel = function()
			{
				jQuery('#depiction').removeClass('invisible');
				jQuery('#jme').addClass('invisible');
			}
			
			function onResultsSuccess() {
				var stop = "";
			
			}
			
			
			
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
			
			
			
		</script>
		<title>Record editor</title>
		<style type="text/css">
			INPUT {border: 1px solid black; padding: 2px 2px 2px 2px;}
			.options TD {padding-right: 5px; padding-bottom: 20px;}
		</style>
		<table width="100%">
			<tr><td class="itunes-up">
				<h1>Apply the model</h1>
			</td></tr>
			<tr><td class="itunes-right">
				<xsl:if test="//message">
						<div class="warning">
							<xsl:value-of select="//message/message"/>
						</div>
					</xsl:if>
				<h1>Provide the compound(s)</h1>

				Please provide compounds for which you want to predict the target property<br/>Several options are available:<br/><br/>
				
				<form name="modeldata" id="modeldata" action="modelquick/apply.do" method="post" enctype="multipart/form-data">
				<div>
					<input type="hidden" name="moldata" send="1" value=""/>
					<table class="options">
						<tr>
							<td><input type="radio" name="data" value="externalfile"/></td>
							<td>Upload compounds from a file<br/> <i>(SDF/MOL2/SMILES/Excel sheet)</i></td>
							<td><input type="file" name="externalfile"/></td>
						</tr>
						<tr>			
							<td><input type="radio" name="data" value="name"/></td>
							<td>Provide a Name/CAS-RN/SMILES</td>
							<td><input type="text" name="name" value=""/></td>
						</tr>
						<tr>			
							<td><input type="radio" name="data" value="moldata" checked="checked"/></td>
							<td>Draw Molecule<br/><i>(click on depiction to the right to draw)</i></td>
							<td id="molecule"><a action="editmolecule"><img id="depiction" src="depiction.jsp?id=947" width="100px" height="100px"/></a></td>
							<td id="jme" valign="top" width="1px" height="1px" class="invisible">
								<applet code="JME.class" name="JME" archive="{@web-root}applets/jme.jar" width="400px" height="300px">
								<div id="javaexc">
									<p>Because you don't have Java in you machine you can't use JME editor. Please Install java to use JME editor but you can use 
									<b>submit smiles and upload file option</b>
									</p>
								</div>
								</applet>
								<div>
									<input id="selectButton" type="button" action="select" value="Select" />
									<input id="cancelButton" type="button" action="cancel" value="Cancel" />
								</div>	
							</td>
						</tr>
					</table>
					<br/>
				</div>
				<div id="nextButton" class="formsubmit">
					<input type="button" action="submit" value="Next&gt;&gt;"/>
				</div>
				<div id="progress" class='invisible'>
					<img src='img/roller.gif' id='progress-img'/><br/>
					<span id='status'>Starting...</span><br/>				
				</div>
				</form>			
			</td></tr>
		</table>
	</xsl:template>
</xsl:stylesheet>