<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />
	
	<xsl:template name="content">
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/browser.js"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/blocks/message.js"></script>
		<style type="text/css">
			BODY {background-image: url(/img/envelope.jpg); background-repeat: no-repeat;}
			.content TABLE {padding: 10px; width: 90%; margin: 10px; }
			TD {padding: 2px 2px 2px 2px;}
			.compact-item TD {border: none;}
			INPUT, TEXTAREA {border: 1px solid black; padding: 2px 2px 2px 2px; margin: 1px 1px 1px 1px;}
			.right {align:right;}
			.content {padding: 30px;}
		</style>
		<title>Message</title>
		<div class="content" style="min-height: 450px;">
		<h1>Message via internal post</h1>
		<xsl:if test="not(//session/user/@login)">
			Please note, that you are logged in as a guest and we dont have your contact info.<br/>If you want to be contacted back by the receiver, please leave your contact details in your message. <br/><br/>
		</xsl:if>
		<table>
		<tr><td width="150">
		<input type="hidden" name="org_msg_id" value="{//message/@org_message_id}" send="1"/>
		From:</td><td><input disabled="disabled" type="text" name="sender" value="{//message/sender}"/>
		</td></tr>
		<tr><td>
		To:</td><td><input disabled="disabled" type="text" name="receiver" value="{//message/receiver}" send="1"/><xsl:if test="//message/receivername = //message/sendername">
			Are you sending a message to yourself? Nice.
		</xsl:if>
		</td></tr>
		<tr><td>
		Subject:</td><td><input type="text" name="subject" value="{//message/subject}" style="width: 100%;" send="1"/>
		</td></tr>
		<tr><td colspan="2">
			<br/>Dear <xsl:value-of select="//message/receivername"/>,
			<textarea id="body" name="body" send="1" style="width: 100%; height: 120px;">
			</textarea>
			Sincerely,<br/>
			<xsl:value-of select="//message/sendername"/>
		</td></tr>
		</table>
		<div class="formscope popup-footer"><a action="send">send</a><a href="javascript:window.closeTab();">cancel</a></div>
		<div id="waitingDialog"> 
	   		<div class="hd">Please wait</div> 
		    <div class="bd" style="text-align: center;"> 
		        Please wait until message is send.<br/>
		        <img src="img/roller_small.gif"/> 
		    </div> 
		</div>
		
		<div id="alertDialog"> 
		    <div class="hd">Message</div> 
	   		 <div class="bd"> 
	           <div id="alertMessage">message</div>
		    </div> 
		</div>
		</div>
	</xsl:template>
</xsl:stylesheet>