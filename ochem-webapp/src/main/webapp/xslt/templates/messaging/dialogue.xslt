<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />

	<xsl:template name="content">
	<title>Intenral messaging</title>
	<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<style type="text/css">
			.dialogue, BODY, HTML {
				height: 100%;
			}
			.dialogue .messages{
				width: 800px;
				border: 1px solid gray;
				vertical-align: top;
				padding-bottom: 170px;
				padding-top: 100px;
			}
			
			.header {
				color: #555;
				background-color: #F2F2F2;
				width: 800px;
				border-bottom: 1px solid gray;
			}
			.dialogue .header {
				position: fixed;
				top: 0px;
			}
			.header > DIV {padding: 20px;}
			
			.messages DIV {
				padding: 10px;
				color: #222;
			}
			
			.dialogue TEXTAREA {
				width: 720px;
				height: 100px;
				margin-bottom: 5px;
				border: 1px solid gray;
				 resize: none;
				 color: #222;
				 font-family: Tahoma, Arial;
			}
			
			.dialogue .input > DIV {
				margin: 15px 45px;	
			}
			
			.dialogue .input {
				background-color: #F2F2F2;
				position: fixed;
				bottom: -2px;
				width: 800px;
				border-top: 1px solid gray;
				color: #888;
			}
			
			.dialogue .header {
				
			}
			
			.dialogue A, .dialogues A {color: #2B587A;}
			
			.messages TD {padding-right: 15px;}
			.messages .user {
				text-align: right;
				vertical-align: top;
			}
			
			.dialogues > DIV {
				width: 800px;
				height: 80px;
				border: 1px solid gray;
				border-top: none;
				vertical-align: middle;
				clear: right;
				padding: 0px 20px;
				overflow: hidden;
			}
			
			.dialogues {
				color: #333;
				font-family: Tahoma, Arial;
				
			}
			.dialogues DIV:hover {
				background-color: #EDF1F5;
				cursor: pointer;
			}
			
			.new {
				float: left; clear: both; background-color: #500; color: white; font-size: 9pt; margin-right: 10px; padding: 2px;
			}
			
			.read-false .message-text {
				font-weight: bold;
			}
			.own.read-false .message-text {
				font-weight: normal;
			}
			
			.preview {
				height: 80px;
				overflow: hidden;
			}
			.preview > DIV {
				overflow: hidden;
				max-height: 70px;
			}
			
			.time {text-align: right; color: #999;}
			.hint-blur {color: #888 !important;}
		</style>
		
		<table class="all-dialogues section invisible">
			<tr>	
				<td class="header"><div>Click on any dialog below to open it</div></td>
			</tr>
			<tr>
				<td class="dialogues"></td>
			</tr>
		</table>
		
		<table class="dialogue section invisible" height="100%">
			<tr>
				<td class="header">
					<div>
					Dialogie with <span></span>
					<br/>
					<a href="dialogue/dialogue.do?render-mode=popup">&#8592; View all dialogues</a>
					</div>	
				</td>
			</tr>
			<tr height="100%">
				<td class="messages">
					
				</td>
			</tr>
			<tr>
				<td class="input">
					<div>
					<textarea id="message-text" hint="Enter your message here"></textarea><br/>
					<a href="javascript:void()" type="button" class="fb-button" id="send" value="Send">Send</a> (Ctrl + Enter)
					</div>
				</td>
			</tr>
		</table>
		<div>
			<div>
			</div>
			
		</div>
		
		
		<script language="javascript">
			include.plugins('view');
			Chat = function() {
				Actionable.call(this);
				var self = this;
				this.dialogueUser = "novserj";
				
				var updateDialogues = function() {
					self.ajax.send({
						url: "dialogue/dialogues.do",
						success: function(response) {
							console.log(response);
							$(".dialogues").html("");
							var dialogues = array(response.others.dialogue);
							for (var i = 0; i &lt; dialogues.length; i++)
							{
								var msg = getTemplate("dialogue-template").render(dialogues[i]);
								$(".dialogues").append(msg);
							}
							
							$(".dialogues > DIV").click(function(){
								location.href = "dialogue/dialogue.do?user=" + $(this).attr("opponent");
							});
							
							if (dialogues.length == 0)
								$(".dialogues").append(getTemplate("no-dialogues").render({}));
							
							$(document).trigger("DOM_updated", $(".messages"));
						}
					});	
				}
				
				var update = function() {
					$(".dialogue .header SPAN").html(self.dialogueUser);
					self.ajax.send({
						url: "dialogue/messages.do",
						data: "user=" + self.dialogueUser,
						success: function(response) {
							$(".dialogue .messages").html("");
							var messages = array(response.dialogue.messages);
							for (var i = 0; i &lt; messages.length; i++)
							{
								var msg = $(getTemplate("message-template").render(messages[i]));
								if (messages[i].sender != self.dialogueUser)
									msg.addClass("own");
								else
									msg.addClass("incoming")
								$(".dialogue .messages").append(msg);
							}
							
							$("html, body").animate({ scrollTop: $(document).height() }, "slow");
							$(document).trigger("DOM_updated", $(".messages"));
							
							setTimeout(function(){
								self.ajax.send({
									url: "dialogue/markRead.do?user=" + self.dialogueUser,
									success: function() {
										$(".incoming.read-false").removeClass("read-false");
									}
								});
							}, 2000);
						}
					});	
				}
				
				var send = function() {
					self.ajax.send({
								url: "dialogue/send.do",
								data: "receiver=" + self.dialogueUser + "&amp;text=" + URLEncode($("#message-text").val()),
								success: function() {
									$("#message-text").val("");
									update();
								}
							});
				}
				
				this.openDialogue = function(user) {
					
				}
				
				this.showAllDialogues = function(user) {
				
				}
				
				this.selectSection = function(section) {
					$(".section").addClass("invisible");
					$(".section." + section).removeClass("invisible");
				}
				
				$(function(){
				
					// Message sending
					$("#message-text").keydown(function(event) {
						if (event.which == 13 &amp;&amp; event.ctrlKey)
							send();	
					});
					$("#send").click(send);
					$("textarea").hint();
					
					if (getParams["user"])
					{
						chat.dialogueUser = getParams["user"];
						update();
						chat.selectSection("dialogue");
					}
					else
					{
						updateDialogues();
						chat.selectSection("all-dialogues");
					}
				});
				
			}
			
			var templatesMap = new Object();
			var getTemplate = function(id) {
				if (!templatesMap[id])
					return templatesMap[id] = new View({element: id});
				else
					return templatesMap[id];
			}
			
			
			var chat = new Chat();
			
			
			
		</script>
		
		<script type="text/template" id="message-template">
			<div message="[%=id %]" class="read-[%=isRead %]">
				<table width="100%">
					<tr>
						<td width="130" class="user">	
							<img src="img/icons/user-16.png"/><a tab="User profile" href="user/profile.do?login=[%=sender%]">[%=sender%]</a>
						</td>
						<td class="message-text">
							[%=replaceHTMLTags(text) %]
						</td>
						<td class="time">
							[%=time %]
						</td>
					</tr>
				</table>
				
			</div>
		</script>
		
		<script type="text/template" id="dialogue-template">
			<div opponent="[%=opponent.login %]">
				<table width="100%" height="100%">
					<tr>
						<td width="120" class="user">	
							<img src="img/icons/user-16.png"/><a tab="User profile" href="user/profile.do?login=[%=opponent.login %]">[%=opponent.login %]</a>
						</td>
						<td class="preview">
							<DIV>
							[% if (unread == "true") { %]<div class="new">NEW</div>[% } %]
							[%=replaceHTMLTags(lastMessage.text) %]
							</DIV>
						</td>
						<td class="time">
							[%=lastMessage.time %]
						</td>
					</tr>
				</table>
			</div>
		</script>
		
		<script type="text/template" id="no-dialogues">
			You dont have any incoming or outgoing messages yet.<br/>
		</script>
		
	</xsl:template>
</xsl:stylesheet>