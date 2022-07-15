function DynamicSelect(container, url, itemElement, scope)
{
	this.ajax = new QSPR.Ajax(url);
	var self = this;
	this.ajax.owner = this;
	EventDispatcher.call(this);
	
	this.itemElement = itemElement;
	
	
	if (container instanceof jQuery)
		this.select = container.get(0);
	else
		this.selectName = container;
	this.scopeBlock = (scope) ? scope : $(document);
	this.fixedItemsCount = 0;
	
	this.update = function(data, defaultValue, event_name)
	{
		var dfd = new jQuery.Deferred();
		if (!self.select)
		{
			self.scopeBlock = (self.scopeBlock) ? self.scopeBlock : $(document);
			self.select = self.scopeBlock.find("select[name='" + self.selectName + "']").get(0);
		}
		
		if (self.select)
		{
			// Show "loading..."
			self.select.options.length = self.fixedItemsCount;
			self.select.options[self.select.options.length]
			                    = new Option("Loading...", "", false, false);
		}
		
		
		self.ajax.send(
			{
				data: data,
				success:
				function(xml)
				{
					self.scopeBlock = (self.scopeBlock) ? self.scopeBlock : $(document);
					if (!self.select)
						self.select = self.scopeBlock.find("select[name='" + self.selectName + "']").get(0);
					if (!self.select)
						window.alert("Dynamic select error: could not select "+self.selectName + " - reload form.");
					// Remove old items
					self.select.options.length = self.fixedItemsCount;
    				
    				var elems = array(xml.list[this.itemElement]);
    					
					for (var i=0; i < elems.length; i++)
					{
						var opt = new Option(elems[i].name || "", elems[i].id, false, false);
						opt.elem = elems[i];
						self.select.options[self.select.options.length] = opt;
					};
					
					if (defaultValue)
						$(self.select).val(defaultValue);
					
					if (event_name)
						self.fireEvent(event_name);
					else
						self.fireEvent("update");
					dfd.resolve();
				},
				error: function()
				{
					dfd.reject();
					self.select.options.length = 0;
					self.select.options[0] = new Option("Error loading options", "", false, false);
				}
			}
		);
		
		return dfd.promise();
	}
	
	$(document).ready(function(){
		if (!self.select)
			self.select = self.scopeBlock.find("select[name='" + self.selectName + "']").get(0);
		self.fireEvent("update");
	});
}

function UnitSelect(selectName, url, itemElement, scope)
{
	var self = this;
	DynamicSelect.call(this, selectName, url, itemElement, scope);
	this.listenEvent("update", function(){
		if (self.select.options.length == 0 || (self.select.options.length == 1 && self.select.options[0].text == "") || (self.select.options.length == 2 && self.select.options[0].value == "-1" && self.select.options[1].value == "9"))
			$(self.select).addClass('invisible');
		else
			$(self.select).removeClass('invisible');
	});
}