
/*
 * A client-side class to track the status of long operations 
 * ...e.g., batch upload, discretize basket, create multiple models, apply models, delete records...
 * TODO Refactoring: All the long operations should use this class in future
 * 
 * Alpha version
 * 
 * @option url - URL to start the operation
 * @option formScope - the data from this selector will be appended to the request
 * @option data - additional data to send
 * @option finished - callback function
 * @option error - callback function
 * @option tracker = (dialog | empty | jQuery object)
 * 
 * by Midnighter
 */

function LongOperation(options)
{
	this.options = options;
	var self = this;
	
	this.operationId = Math.floor(Math.random() * 100000000);
	
	// Start an operation with a randomly generated identifier
	this.start = function()
	{
		LongOperation.currentOperation = this;
		
		// Generate a new operation ID
		this.operationId = Math.floor(Math.random() * 100000000);
		
		LongOperation.map[this.operationId] = this;
		
		var data = ["operation-id=" + this.operationId, "out=json"];
		if (this.options.formScope)
		{
			$(this.options.formScope).find("input[name], textarea[name]").each(function(){
				data.push($(this).attr("name") + "=" + URLEncode($(this).val()));
			});
		}
		var queryString = data.join("&");
		
		if (this.options.data)
			queryString += "&" + this.options.data.replace("operation-id", "tmp").replace("start-operation", "tmp");
		
		console.log(queryString);
		
		$.ajax({
			type: "POST",
			url: this.options.url,
			data: queryString,
			success: function()
			{
				self.waitingDialog.start();
				self.startCheckingStatus();
			}
		});
	}
	
	this.startCheckingStatus = function()
	{
		self.waitingDialog.setStatus("Starting...");
		LongOperation.currentOperation = this;
		LongOperation.map[this.operationId] = this;
		self.timer = setInterval('LongOperation.map[\''+self.operationId+'\'].checkStatus()', 1000);
		self.waitingDialog.show();
	}
	
	this.stopCheckingStatus = function()
	{
		clearInterval(self.timer);
		self.waitingDialog.hide();
	}
	
	// Retrieve and display the operation status
	this.checkStatus = function()
	{
		$.ajax({
			type: "POST",
			url: "longoperations/operationStatus.do",
			data: "out=json&nodb=1&operation-id=" + this.operationId,
			dataType: "json",
			success: function(response)
			{
				var status = response.message.message;
				if (!status)
					return;
				self.waitingDialog.setStatus(status);
				if (status.indexOf("Error") == 0)
				{
					clearInterval(self.timer);
					self.waitingDialog.stop();
					self.waitingDialog.error();
					if (self.options.error)
						self.options.error(status.replace("Error: ", ""));
				} 
				else if (status.indexOf("Finished") == 0)
				{
					clearInterval(self.timer);
					self.waitingDialog.hide();
					if (self.options.finished)
						self.options.finished(status);
				}
			},
			error: function(response)
			{
				// Could not request status. Do nothing. Not a big deal here
			}
		});
	}
	
	// Request cancellation
	this.cancel = function()
	{
		$.ajax({
			url: "longoperations/cancelOperation.do",
			data: "out=json&nodb=1&operation-id=" + this.operationId,
			success: function()
			{
				self.waitingDialog.setStatus("Cancel requested...");
				$("#cancel-button").setDisabled(true);
			},
			error: function()
			{
				
			}
		});
	}

	if (!this.options.tracker || this.options.tracker == "dialog")
		this.waitingDialog = new WaitingDialog();
	else
		this.waitingDialog = new SimpleTracker(this.options.tracker);
}

LongOperation.map = new Object();


function WaitingDialog() // "implements" Tracker
{
	// The waiting dialog HTML. Works, but might have to be refactored
	$("body").prepend($('<div id="operationWaitingDialog"><div class="hd">Please wait</div><div class="bd" style="text-align: center;"><span class="status">Starting...</span><br/><img src="img/roller_small.gif"/> <br/><input type="button" value="OK" class="invisible" id="ok-button" onclick="LongOperation.currentOperation.waitingDialog.hide()"/><input type="button" id="cancel-button" value="Cancel" onclick="LongOperation.currentOperation.cancel();"/></div></div>'));

	this.yuiDialog = new YAHOO.widget.Dialog("operationWaitingDialog", { 
    	width:"400px", 
		fixedcenter:true, 
		modal:true, 
		visible:false, 
		close: false
    });
	this.yuiDialog.render();
	
	this.setStatus = function(status)
    {
		$("#operationWaitingDialog").find(".status").html(status.replace(/\n/g, "<br/>"));
    }

	this.stop = function(status)
    {
    	$("#operationWaitingDialog").find("img").addClass("invisible");
    	$("#operationWaitingDialog").find("#cancel-button").addClass("invisible");
    	$("#operationWaitingDialog").find("#ok-button").removeClass("invisible");
    }
	
	this.start = function()
	{
		$("#operationWaitingDialog").find("img").removeClass("invisible");
    	$("#operationWaitingDialog").find("#cancel-button").removeClass("invisible");
    	$("#operationWaitingDialog").find("#ok-button").addClass("invisible");
	}
	
	this.error = function()
	{
		
	}
	
	this.show = function()
	{
		this.yuiDialog.show();
	}
	
	this.hide = function()
	{
		this.yuiDialog.hide();
	}
}

function SimpleTracker(selector)
{
	this.setStatus = function(status)
    {
		selector.html(status.replace(/\n/g, "<br/>"));
    }

	this.stop = function(status)
    {
		selector.parent().find(".lo-progress-bar").addClass("invisible");
    }
	
	this.start = function(status)
	{
		selector.parent().removeClass("lo-error");
		selector.parent().find(".lo-progress-bar").removeClass("invisible");
	}
	
	this.error = function()
	{
		selector.parent().addClass("lo-error");
	}
	
	this.show = function()
	{
		selector.parent().removeClass("invisible");
	}
	
	this.hide = function()
	{
		selector.parent().addClass("invisible");
	}
}
