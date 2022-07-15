/**
 * The parent class for all DHTML client side "browsers" of items
 * 
 * @param url (e.g., 'browser/list.do')
 * @param controller (e.g. browser)
 * @param itemElement
 * @param itemTemplate 
 * 
 * @returns
 */
function Browser()
{
	Actionable.call(this); // superclass(es)
	EventDispatcher.call(this);
	
	this.container = "Browser";
	this.itemClass = "browser-item";
	this.htmlWaiting = '<img src="img/long_green.gif"/>';
	this.page = 1;
	this.itemElement = "items";
	this.pager = new DirectPager(11, this.useHistory);
	this.resetFiltersOnLoad = true;
	if (this.scope)
		this.pager.selectors.scope = this.scope;
	
	this.additionalRequestData = "";
	
	this.filters = new Filters();
	this.ready = false;
	this.useTableRows = false;
	this.options = {"editwindow_scrollbars":1, "editwindow_resizable":"yes", "highlight_row": "yes"};
	var self = this;

	if (this.controller)
	{
		this.url = this.controller+"/list.do";
		this.actionURL = this.controller+"/action.do";
	}
	
	this.loadItemsAndPagerFromJSONList = function(list)
	{
		this.mainDiv.html('');
		var items;
		
		if (list)
			items = array(list[this.itemElement]);
		
		// Cool drawing of items
		if (items)
			for (var i = 0; i < items.length; i++)
				self.drawFromJSON(items[i]);
		
		// Draw pager
		var pager;
		if (list)
		{
			self.pager.load(list);
			self.pager.render();
			if (list.sortingDiscarded == "true")
				$("#sorting-discarded").html("Sorting has been discarded because your query was too broad. Make a more specific query to get sorted results");
			$("#sorting-discarded").setClass("invisible", list.sortingDiscarded != "true");
		}
	}
	
	this.request = function(resetPage, onSuccess)
	{
		this.fireEvent("items_load");
		
		if (resetPage)
			this.page = 1;
			
		this.mainDiv
			.html(this.htmlWaiting);
			
		// Ignore previously sent requests
		if (this.ajaxRequest)
			this.ajaxRequest.success = this.ajaxRequest.error = function(){};
			
		if (this.operation)
			this.operation.stopCheckingStatus();
		
		if ($("#query-status").length > 0 && typeof LongOperation === 'function')
		{
			this.operation = new LongOperation({tracker: $("#query-status")});
			this.filters.setValue("operation-id", this.operation.operationId);
			this.filters.setValue("start-operation", 1);
		}
		
		var _data = "pagenum=" + this.page;
		_data += this.filters.getQueryString();
		_data += "&" + this.getAdditionalData();
		this.ready = false;
		
		
		
		this.ajax.send
		(
			this.ajaxRequest = {
				url: this.url, 
				data: _data, 
				success: function(response)
				{
					this.pager.setVisible(true);
					if (this.ajax.dataType == "xml")
					{
						response = xml2json(response, "\t");
						response = eval('('+response+')');
						response = response.model;
					}
					
					// Cool drawing of Filters
					if (response.filters)
					{
						var filters = array(response.filters.filter);
						if (this.resetFiltersOnLoad)
							this.filters.reset();
						for (var i = 0; i < filters.length; i++)
							this.filters.setValue(filters[i]['name'], filters[i]['value'], filters[i]['title']);
					}
					
					this.responseExtras = response.others || {};
					if (response.param) {
						var params = array(response.param);
						for (var i = 0; i < params.length; i++)
							this.responseExtras[params[i].key] = params[i].content;
					}
					
					this.list = response.list;
					
					this.fireEvent("before_draw");
					this.loadItemsAndPagerFromJSONList(response.list);
					
					this.selectionSize = response.selectionSize;
					this.setPositionById(getParams['selected'] ? getParams['selected'] : undefined);
					this.ready = true;
					if (onSuccess != undefined)
						onSuccess.call(this);
					this.fireEvent("items_loaded");
					$(document).trigger("DOM_updated", $(document));
				},
	
				error: function(e)
				{
					this.mainDiv
					.html('<div style="padding: 10px 10px; margin-top: 10px; color: #550000; font-size: 14pt; background-color: #FEE;">Error: '+e+'</div>');
					this.pager.setVisible(false);
				}
			}
		);
		
		if (this.operation)
			this.operation.startCheckingStatus();
	}
	
	
	this.getAdditionalData = function() //A stub, if a browser needs to add some data
	{
		return self.additionalRequestData;
	}
	
	
	this.drawFromXML = function(entity, options)
	{
		var json = xml2json(entity.get(0), "\t");
		json = eval("("+json+")");
		var rendered = this.drawFromJSON(json[this.itemElement], options);
	}
	
	this.drawFromJSON = function(JSONEntity, options)
	{
		var id = this.getRecordIdentifier(JSONEntity);
		var block = this.useTableRows ?
			$('<tr class="'+this.itemClass+'" rec-id="' + id + '"/>') : 
			$('<div class="'+this.itemClass+'" rec-id="' + id + '"/>');
		var container = (options && options.container ? options.container : this.mainDiv);
			
		var rendered = this.draw(JSONEntity);
		block.append(rendered);
		block.append("\r\n");
		block.get(0).entity = JSONEntity;
	
		this.attachActions(block);
		options = options || {};
		var xpend = (options.method == 'prepend') ? container.prepend : container.append;
		if (options.replaceId)
			this.mainDiv.find('[rec-id="'+options.replaceId+'"]').replaceWith(block);
		else
			xpend.call(container, block);
		
		this.setPositionById(id);
		if (this.onItemDrawn)
		{
			this.onItemDrawn();
		}
		if (!window.callback)
			this.currentBlock.find('[action="select"]').remove();
		else
			this.currentBlock.find('[action="select"]').removeClass("invisible");	
		
		if (this.options.highlight_row == "yes")
		{
			this.currentBlock
				.mouseover(function(){$(this).addClass("softhighlight");})
				.mouseout(function(){$(this).removeClass("softhighlight");});
		}
		
	}
	
	this.draw = function(entity)
	{
		if (!this.view)
			this.view = new View({url: this.itemTemplate});
		return this.view.render(entity);
	}
	
	this.getActionQuery = function()
	{
		if (this.currentRecordId)
		{
			var query = "id="+this.currentRecordId+this.filters.getQueryString();
			this.currentBlock.find("*[send]").each(function()
			{
				var value =  ($(this).attr('value') == undefined) ? "" : URLEncode($(this).attr('value'));
				query += "&" + $(this).attr('name') + "=" + value;
			});
			return query;
		}
		else
			return this.filters.getQueryString();
	}
	
	this.setPosition = function(element)
	{
		this.currentBlock = $(element).parents('[rec-id]');
		this.setPositionById(this.currentBlock.attr('rec-id'));
	}
	
	this.setPositionById = function(id)
	{
		this.currentRecordId = undefined;
		this.currentEntity = undefined;
		this.currentBlock = undefined;
		this.mainDiv.find("[rec-id]").removeClass("highlighted");
		
		if (id != undefined)
		{
			this.currentBlock = this.mainDiv.find('[rec-id="'+id+'"]');
			this.currentRecordId = id;
			if ((id != undefined) && this.currentBlock.get(0))
			{
				this.currentRecordId = this.currentBlock.attr('rec-id');
				this.currentEntity = this.currentBlock.get(0).entity;
			}
			
			// Highlight current block
			if (this.currentBlock != undefined)
			{
				this.currentBlock.addClass("highlighted");
			}
		}
	}
	
	
	this.getRecordIdentifier = function(entity)
	{
		if (entity.attr)
			return entity.attr('id');
		else if (entity.id)
			return entity.id;
		else
			return entity;
	}
	
	this.deleteRecord = function(id)
	{
		if (!id)
			id = this.currentRecordId;
		this.mainDiv.find('[rec-id="'+id+'"]').remove();
	}
	
	this.initialize = function(dontRunRequest)
	{
		console.log("Initializing the browser " + this.container);
		this.mainDiv = $('#'+this.container);
		
		// Surprisingly, "div#id"  selector doesn't work with Chrome and Safari
		// So we use "#id" insetad of "div#id" / Midnighter
		//this.mainDiv = $('div#'+this.container);
		
		// Initialize filters
		//if (this.filters != undefined)
		//	Filters.fill(this.filters);
		
		Autocomplete.updateFromHTML();
		if (!dontRunRequest)
			if (!this.pager.useHistory)
				this.request();
	}
	
	this.doSelect = function()
	{
		window.callback(this.currentEntity);
	}
	
	this.doReset = function()
	{
		this.filters.reset();
		this.page = 1;
		this.request();
	}
	
	this.doPopup = function(link)
	{
		var win = openTab($(link).attr('tab'), webRoot+$(link).attr('link'));
	}
	
	this.doEdit = function(link)
	{
		var id = (this.currentRecordId) ? this.currentRecordId : "-1";
		var editQueryString = this.editQueryString || "";
		
		//some time its undefined and because of that their was error in article browser
		if($(link).attr("query"))
			editQueryString += $(link).attr("query");
		var editController = this.editController ? this.editController : this.controller;
		var printedItemName = this.printedItemName || this.itemElement;
		
		var win = openTab("Edit "+printedItemName, webRoot + editController + "/edit.do?render-mode=popup&id=" + id + "&" + editQueryString);
		win.callback = function(entity) 
		{
			if (self.onItemSaved)
			{
				self.onItemSaved(entity);
			}
				
			win.closeTab();
			self.request(false, function()
			{
				var realEntity = entity[this.itemElement] || entity;
				self.setPositionById(realEntity.id);
			
			});
			
			
		};
	}	
	
	this.filters.onFilterChange = this.doRefresh = function()
	{
		self.request(true);
	}
	
	// Override
	this.pager.onStateChanged = function() 
	{
		self.filters.values.pagesize = self.pager.itemsPerPage;
		self.page = self.pager.currentPage;
		self.request();
	};
}

// UnifiedBrowser

function UnifiedBrowser()
{
	Browser.call(this);
	
	this.parentActionQuery = this.getActionQuery;
	
	this.getActionQuery = function(actionName)
	{
		var query = "";
		if (actionName == "saveall" ) //|| actionName == "savebasket")
		{
			query = this.filters.getQueryString();
			this.mainDiv.find("*[send]").each(function()
			{
				var value =  ($(this).attr('value') == undefined) ? "" : URLEncode($(this).attr('value'));
				query += "&" + $(this).attr('name') + "=" + value;
			});
		}
		else
			query = this.parentActionQuery();
		return query;
	}
	
	this.onSaveallSuccess = function()
	{
		this.request();
	}
	
	this.setValue = function(name, value)
	{
		this.currentBlock.find('[name="'+name+'"]').attr('value', value);
	}

	this.getValue = function(name)
	{
		return this.currentBlock.find('[name="'+name+'"]').attr('value');
	}
	
	this.getItemsCount = function()
	{
		return this.mainDiv.find('[rec-id]').length;
	}
}

function Filters() 
{
	var self = this;
	this.scope = document;
	this.firstTime = true;
	this.useUrlParameters = true;
	
	// Internal non-DOM filters
	this.values = new Object();
	this.titles = Object();
	
	this.setFromUrl = function() {
		for (key in getParams)
			if (key != "include")
				self.setValue(key, getParams[key]);
		this.useUrlParameters = false;
	}
	
	this.getQueryString = function()
	{
		_data = ""; 
		
		if (this.firstTime)
		{
			this.attachHandlers();
			
			if (this.useUrlParameters)
			for (key in getParams)
				if (key != "include")
				{
					_data += '&'+key+'='+getParams[key];
				}
			//this.reset();
			this.firstTime = false;
		}

		// Internal non-DOM filters
		for (var key in this.values)
			_data += '&'+key+'='+URLEncode(self.values[key]);
		
		$(this.scope).find('*[filter=1]').each
			(
				function()
				{
					var allowEmpty = $(this).attr('allowempty');
					
					var value;
					var type = $(this).attr('type');
					
					if (type == "checkbox")
						if (!this.checked)
							return;
						else
							value = "checked";
					else
						value = $(this).attr('value');
					
					if (type == "radio" && !this.checked)
						return; // do not sent unchecked radio-buttons
					
					
						
					if (value != undefined && value != $(this).attr('prompt') && (value != "" || allowEmpty=="1"))
					{
						if ( value != undefined ) {
							//value = value.replace(/^\s+|\s+$/g,"");
						}
						_data += '&'+$(this).attr('name')+'='+URLEncode(value.trim());
						$(this).addClass("highlighted");
					}
					else
						$(this).removeClass("highlighted");
				}
			);
		return _data;
	}

	this.getValue = function(name)
	{
		var target = $(this.scope).find('[filter=1][name="'+name+'"]');
		if (target.length > 0)
			return target.attr('value');
		else
			return this.values[name];
	}
	
	this.isSet = function(name)
	{
		var val = this.getValue(name);
		return val != undefined && val != "" && val != -1;
	}
	
	this.reset = function()
	{
		$(this.scope).find('input[filter][type="text"],select[filter]').each(function() { //Experimental
				var oldValue = $(this).attr('value');
				$(this).attr('value', '');
				$(this).attr('oldValue', oldValue);
		});
		
		$(this.scope).find('[filter] input').each(function () {
			var oldValue = $(this).attr('value');
			$(this).attr('value', '');
			$(this).attr('oldValue', oldValue);
		});
		
		$(this.scope).find('input[filter][type="checkbox"]').removeAttr('checked');
		$(this.scope).find('[filter]').removeClass("highlighted");
	}

	this.setValue = function(name, value, title)
	{
		var filter = $(this.scope).find('[filter][name="'+name+'"]');
		var filterDOM = filter.get(0);
		
		if (title)
			this.titles[name] = title;
		
		if (filterDOM && filterDOM != 'undefined')
		{
			// This is a regular filter in DOM
			if (filter.attr('type') == "radio")
			{
				filter.removeAttr("checked");
				filter.filter("[value='"+value+"']").attr("checked", true);
			}
			else if (filter.attr('type') == "checkbox")
				filter.attr('checked', true);
			else
			{
				filter.val(value);
				if (filterDOM.autocomplete)
				{
					filterDOM.autocomplete.setValue(title);
					filterDOM.autocomplete.setKey(value); //Experimental
				}
				if (filterDOM.tagName.toLowerCase() == "a")
					filterDOM.innerHTML = title;
				
				filter.addClass("highlighted");
			}
		}
		else
			this.values[name] = value;
		
	}
	
	this.attachHandlers = function()
	{
		$(this.scope).find('input[filter]').keypress(
			function(e)
			{
				if ((e.which == 13) && (self.onFilterChange))
				{
					self.onFilterChange();
					e.stopPropagation();
					e.preventDefault();
					return false;
				}
			}
		);
		$(this.scope).find('select[filter], input[filter=1][type="checkbox"], input[filter=1][type="radio"]').change(
			function(e)
			{
				self.onFilterChange(this);	
			}
		);
		
	}
}