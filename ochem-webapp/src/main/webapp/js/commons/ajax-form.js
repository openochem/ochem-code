// AjaxForm : Basic code
// Please think it over before modifying anything here
// Midnighter

function AjaxForm()
{
	var self = this;
	Actionable.call(this);
	FormValueHandler.call(this);

	this.getActionQuery = function()
	{
		var data = "";
		this.scopeBlock.find('*[send]').each
		(
			function()
			{
				if ($(this).attr('type') == "checkbox" && !this.checked)
					return;
				var value = $(this).attr('value');
				if (value != undefined && value != $(this).attr('prompt'))
					data += '&'+$(this).attr('name')+'='+URLEncode(value);
			}
		);
		for (var key in this.fields)
			data += '&'+key+'='+URLEncode(this.fields[key]);
		return data;
	}
	
	this.onEditSuccess = function()
	{
	}
	
	$(document).ready(function(){
		// Restore field titles from hidden fields
		$("[storein]").each(function(){
			var storageElement = $("[name='"+$(this).attr('storein')+"']");
			if (storageElement.val())
				$(this).html(storageElement.val());
		});
	});
}

function EditForm()
{
	AjaxForm.call(this);
	EventDispatcher.call(this);
	
	this.changed = false;
	var self = this;
	
	$(document).ready(function()
	{
		$('input[type="text"][send]').change(function()
		{
			if (self.changed != true)
			{
				self.changed = true;
				self.fireEvent('formchanged');
			}
		});
	});
	
	this.onEditSuccess = function(data)
	{
		this.entity = data;
		this.onItemSaved.call(this, data);
	};
	
	this.onItemSaved = function(entity)
	{
		if (window.callback)
			window.callback(entity);
		else
			window.alert('Item has been successfully saved');
	}
}