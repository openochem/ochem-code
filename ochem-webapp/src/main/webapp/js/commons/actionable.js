/* A class to attach action to HTML links
 * Actions can be both local and remote 
 * */
function Actionable()
{
	this.ajax = new QSPR.Ajax();
	this.ajax.owner = this;
	this.restriction = false;
	var self = this;
	
	$(document).ready
	(
		function()
		{
			self.scopeBlock = (self.scope) ? $(self.scope) : $(document);
			try {
			  self.attachActions(self.scopeBlock);
			}
			catch(err) {
			  // do nothing
			}
		}
	);

	this.callAction = function(actionName, link, params)
	{
		var beforeMethod = eval('this.before'+actionName.ucFirst());  // "beforeDelete"
		var actualMethod = eval('this.do'+actionName.ucFirst());  // "doDelete"
		if (!actualMethod)
			actualMethos = eval('this.do'+capitaliseFirstLetter(actionName));
		var onSuccess = eval('this.on'+actionName.ucFirst()+'Success'); // "onDeleteSucess"
		var onError = eval('this.on'+actionName.ucFirst()+'Error'); // "onDeleteError"
		
		if (!params)
			params = new Object();
		
		if (!params.success)
			params.success = onSuccess;
		if (!params.error)
			params.error = onError;
		
		if (beforeMethod)
			if (!beforeMethod.call(this, link))
				return false;

		// Either direct action or AJAX request
		if (actualMethod != undefined)
			actualMethod.call(this, link);
		else
			this.doAction(actionName, params, link);
		
		if (params.event)
		{
			params.event.preventDefault();
			console.log("PD");
			return false;
		}
		return true;
	}

	this.attachActions = function(block)
	{
		// Bind only to relevant actions
		var restriction = (self.restriction) ? "[restrict='"+self.restriction+"']" : ":not([restrict])";
		block.find('[action]'+restriction).filter(":not(form)").each(function()
		{
			if ($(this).attr('restrict') != undefined && $(this).attr('restrict') != self.restriction)
				return;
			
			if (this.actionable)
				throw("Double action attachment attempt for "+$(this).attr('action')+". Rearrange scopes or put restrictions");
			this.actionable = self;
			$(this).attr('href', 'javascript:void(0)');
			$(this).click(
				function(e)
				{
					this.actionable.setPosition(this);
					var action = $(this).attr('action');
					return this.actionable.callAction(action, this, {event: e});
				});	
		});
	}
	
	// Ajax request action
	this.doAction = function(action, params, link)
	{
		var data = "out=json&action="+action+"&"+this.getActionQuery(action);
		if (link && $(link).attr('ajax-data'))
			data += "&" + $(link).attr('ajax-data');
		if (params.data)
			data += "&" + params.data;
		data = data.replace(/&&/, "&");
		
		var ajaxOptions = { 
				data: data, 
				success: params.success,
				error: params.error
			};
		
		if (this.actionURL)
			ajaxOptions.url = this.actionURL.replace("*", action);
		
		this.ajax.send(ajaxOptions);
	}
	
	this.getActionQuery = function(action)
	{
		// stub to override
		return "";
	}
	
	this.setPosition = function(element)
	{
		// stub to override
	}
}

// Event dispatcher
function EventDispatcher()
{
	var self = this;
	this.eventScope = $('<div id="'+Math.random()+'" class="invisible"/>');
	$(document).append(this.eventScope);

	this.listenEvent = function(event, funk)
	{
		$(this.eventScope).bind(event, null, funk);
	}
	
	this.fireEvent = function(event, data)
	{
		$(this.eventScope).trigger(event, data);
	}
}

function FormValueHandler() {
	
	var self = this;
	this.fields = new Array();
	
	this.getValue = function(name)
	{
		var elem = $("[name='"+name+"']");
		if (elem.length == 0)
			return this.fields[name];
		else if (elem.attr('type') == "checkbox")
			return elem.get(0).checked;
		else if (elem.attr('type') == "radio")
        {
            for (i=0; i<elem.length; i++)
                if ($(elem[i]).is(':checked'))
                    return $(elem[i]).val();
            return null;
        }
        else
			return elem.val();
	}
	
	this.setValue = function(name, val, title)
	{
		// Set value, set title, store title in hidden variable
		// Possibility to store title is provided for restoring values after hitting browser back button. (C) Midnighter.
		var targetElement = $("[name='"+name+"']");
		if (targetElement.length > 0)
		{
			if (targetElement.attr("type") == "radio")
				$("[name='"+name+"']").filter("[value="+val+"]").prop('checked', true);
			else
			{
				targetElement.attr('value', val);
				if (title != undefined)
				{
					var titleElement = $("[bindto='"+name+"']");
					if (titleElement.length > 0)
					{
						$("[bindto='"+name+"']").html(title);
						if (titleElement.attr('storein'))
							$("[name='"+titleElement.attr('storein')+"']").val(title);
					}
				}
			}
		}
		else
			this.fields[name] = val;
	}
}

function createWaitingDialog()
{
	var dialog = new YAHOO.widget.Dialog("waitingDialog", { 
    	width:"325px", 
		fixedcenter:true, 
		modal:true, 
		visible:false, 
		close: false
    });
	
	dialog.setStatus = function(status)
    {
    	$("#waitingDialog").find(".status").html(status);
    }
	
	dialog.render();
	return dialog;
}

String.prototype.ucFirst = function()
{
   // split string
   var firstChar = this.substring(0,1);
   var remainChar = this.substring(1);

   // convert case
   firstChar = firstChar.toUpperCase(); 
   remainChar = remainChar.toLowerCase();

   return firstChar + remainChar;
}

function capitaliseFirstLetter(string)
{
    return string.charAt(0).toUpperCase() + string.slice(1);
}