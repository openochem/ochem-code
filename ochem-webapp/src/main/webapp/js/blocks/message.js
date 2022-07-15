function Message()
{
	AjaxForm.call(this);
	this.controller = "mail";
	this.actionURL = "mail/action.do";
	this.itemElement = "message";
	var self = this;
	
	this.beforeSend = function()
	{	
		if($("#body").attr("value").length < 5)
		{
			window.alert("Empty mail cannot be sent. Please, write the message.");
			return false;
		}
		this.waitingDialog.show();
		return true;
	}

	this.onSendSuccess = function(entity)
	{
		window.closeTab();
	}
	
	this.onSendError = function(entity){
		this.waitingDialog.hide();
		$("#alertMessage").html(entity);
		this.recordDialog.show();
	}
	
	var handleSubmit = function() { 
		this.cancel();
	};
	
	this.initialize = function()
	{
		this.recordDialog = new YAHOO.widget.Dialog("alertDialog", { 
	    	width:"425px", 
			fixedcenter:true, 
			modal:true, 
			visible:false 
	    });
	    
		this.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", { 
	    	width:"325px", 
			fixedcenter:true, 
			modal:true, 
			visible:false, 
			close: false
	    });
	    
	    var myButtons = [{ text:"Ok", handler:handleSubmit}]; 
		this.recordDialog.cfg.queueProperty("buttons", myButtons);
	    
	    this.recordDialog.render();
	    this.waitingDialog.render();
	}
}
var message = new Message();
$(document).ready(function()
{
	message.initialize();	
});