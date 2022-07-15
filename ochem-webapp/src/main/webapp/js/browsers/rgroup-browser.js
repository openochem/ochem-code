function BrowserListHandler()
{
	this.addItemToList = function(JSONItem, list, itemElement)
	{
		if (!list[itemElement])
			list[itemElement] = new Array();
		
		var index = -1;
		
		for (i=0; i<list[itemElement].length; i++)
			if (JSONItem.id == list[itemElement][i].id)
				index = i;
		
		if (index != -1)
			return;
		
		list.lastResult++;
		list.size++;
		list[itemElement].push(JSONItem);
	}
	
	this.removeItemFromList = function(JSONItem, list, itemElement)
	{
		if (!list[itemElement])
			list[itemElement] = new Array();
		
		
		for (i=0; i<list[itemElement].length; i++)
			if (JSONItem.id == list[itemElement][i].id)
				index = i;

		if (index == -1)
			return;
		
		list.lastResult--;
		list.size--;
		list[itemElement].splice(index, 1);	
	}
}

function RGroupBrowser()
{
	var self = this;
	this.controller = "rgroup";
	Browser.call(this);

	this.scope = ".panelFull";
	this.pager.selectors.scope = this.scope;
	this.filters.scope = this.scope;
	this.filters.useUrlParameters = false;
	this.container = "BrowserFull";	
	this.itemElement = "rgroup";
	this.itemClass = "rgroupitem";
	this.itemTemplate = (getParams["edit"] == 1) ? 'js/templates/rgroup-edit.ejs' : 'js/templates/rgroup.ejs';
	this.pager.selectors.pager = ".pager";
	
	TooltipMenu.call(this);
	this.menuSelector = "#ttmenu-f";
	
	this.basketBrowser = null; // <-- Will be initialized outside
	this.activeBlock = null;
	this.parentParentInitialize = this.initialize;

	this.doRenamegroupdialog = function()
	{
		var text = $(self.activeBlock).find("div[class=groupname]").find("i").text();
		$("#renamegroup").find("input[type=text]").val(text);
		$("#renamegroup").dialog("open");
	}
	
	this.doRenamegroupSubmit = function(input)
	{
		var name = URLEncode($(input).val());
		var id = 	$(self.activeBlock).attr("rec-id");
		var params = new Object();
		params.data = "name="+name+"&groupid="+id;
		
		self.callAction("renamegroup", null, params);
		$("#renamegroup").dialog("close");
	}
	
	this.onRenamegroupSuccess = function(data)
	{
		self.request();
	}
	
	this.onBasketaddSuccess = function(data)
	{
		new BrowserListHandler().addItemToList(data[self.basketBrowser.itemElement], self.basketBrowser.list, self.basketBrowser.itemElement);
		self.basketBrowser.loadItemsAndPagerFromJSONList(self.basketBrowser.list);	
		$(this.menuSelector).addClass("invisible").removeClass("softhighlight");
	}
	
	this.doBsktadd = function()
	{
		var params = new Object();
		params.data = "groupid="+$(self.activeBlock).attr("rec-id")+"&basketid="+self.basketBrowser.filters.getValue("basketid");
		self.callAction("basketadd", null, params);
	}
	
	
	this.initialize = function()
	{
		$("#renamegroup").dialog({
			autoOpen: false,
			height: 220,
			width: 280,
			modal: true,
			buttons: {
				"OK": function() 
					{
						self.doRenamegroupSubmit($(this).find("input[type=text]").get(0));					
					},
				Cancel: function() 
					{
						$(this).dialog("close");
					}
			}
		});
		$("#renamegroup").keypress(function(e) {
		    if (e.keyCode == $.ui.keyCode.ENTER)
		    	self.doRenamegroupSubmit($(this).find("input[type=text]").get(0));
		});
		
		this.parentParentInitialize();
	}

}

function RGroupBasketBrowser()
{
	var self = this;
	this.controller = "rgroup";
		
	Browser.call(this);
	this.scope = ".panelBasket";
	this.pager.selectors.scope = this.scope;
	this.filters.scope = this.scope;
	this.filters.useUrlParameters = false;
	this.container = "BrowserBasket";
	this.itemElement = "rgroup";
	this.itemClass = "rgroupitem";
	this.itemTemplate = 'js/templates/rgroup.ejs';
	this.pager.selectors.pager = ".pager";
	this.basketFilter = new DynamicSelect("basketid", "rgroup/basket_list.do", "rgroupbasket");
	
	TooltipMenu.call(this);
	this.menuSelector = "#ttmenu-b";

	this.selectedBasketId = -1;
	
	this.parentParentInitialize = this.initialize;
	
	this.onItemsLoaded = function()
	{
		if (self.pager.totalNum == 0)
			self.mainDiv.html("<div class='content-message'>Empty R-Group basket</div>");
	}
	
	this.doCreatebasketdialog = function()
	{
		$("#createbasket").find("input[type=text]").val("");
		$("#createbasket").dialog("open");
	}
	
	this.doRenamebasketdialog = function()
	{
		$("#renamebasket").find("input[type=text]").val($(this.basketFilter.select).find("option:selected").html());
		$("#renamebasket").dialog("open");
	}
	
	this.doRenamebasketSubmit = function(input)
	{
		var name = URLEncode($(input).val());
		var id = 	$(this.basketFilter.select).val();
		var params = new Object();
		params.data = "basket-name="+name+"&basketid="+id;
		
		self.callAction("renamebasket", null, params);
		$("#renamebasket").dialog("close");
	}
	
	this.doCreatebasketSubmit = function(input)
	{
		var name = URLEncode($(input).val());
		var params = new Object();
		params.data = "basket-name="+name;
		self.callAction("createbasket", null, params);
		$("#createbasket").dialog("close");
	}
	
	
	
	this.doDeletebasketdialog = function()
	{
		var id = 	$(this.basketFilter.select).val();
		var params = new Object();
		params.data = "basketid="+id;
		self.callAction("deletebasket", null, params);		
	}
	
	this.onCreatebasketSuccess = this.onDeletebasketSuccess = this.onRenamebasketSuccess = function(item)
	{
		if (item == undefined || item.rgroupbasket == undefined)
			self.selectedBasketId = -1;
		else
			self.selectedBasketId = item.rgroupbasket.id;
		
		self.basketFilter.update("pagenum=1&pagesize=500",null,"basketlist_update");
	}
	
	this.onBasketremoveSuccess = function(data)
	{
		new BrowserListHandler().removeItemFromList(data[self.itemElement], self.list, self.itemElement);
		self.loadItemsAndPagerFromJSONList(self.list);
		$(this.menuSelector).addClass("invisible").removeClass("softhighlight");		
	}
	
	this.doBsktremove = function()
	{
		var params = new Object();
		params.data = "groupid="+$(self.activeBlock).attr("rec-id");
		self.callAction("basketremove", null, params);
	}
	

	this.onBasketFilterLoaded = this.onBasketFilterChange = function()
	{
	    for (i=0; i<self.basketFilter.select.options.length; i++)
	    {	
	    	if (self.basketFilter.select.options[i].value == self.selectedBasketId)
	    	{
	    		self.basketFilter.select.selectedIndex = i;
	    		break;
	    	}
	    }
	    self.selectedBasketId = -1;
		self.request();
	}
	
	this.basketFilter.listenEvent("basketlist_update", this.onBasketFilterLoaded);
	this.listenEvent("items_loaded", this.onItemsLoaded);
	
	this.initialize = function()
	{
		$(this.basketFilter.select).change(this.onBasketFilterChange);
		this.basketFilter.update("pagenum=1&pagesize=500",null,"basketlist_update");
		
		$("#createbasket").dialog({
			autoOpen: false,
			height: 220,
			width: 280,
			modal: true,
			buttons: {
				"OK": function() 
					{
						self.doCreatebasketSubmit($(this).find("input[type=text]").get(0));					
					},
				Cancel: function() 
					{
						$(this).dialog("close");
					}
			}
		});
		$("#createbasket").keypress(function(e) {
		    if (e.keyCode == $.ui.keyCode.ENTER)
		    	self.doCreatebasketSubmit($(this).find("input[type=text]").get(0));
		});
		
		
		$("#renamebasket").dialog({
			autoOpen: false,
			height: 220,
			width: 280,
			modal: true,
			buttons: {
				"OK": function() 
					{
						self.doRenamebasketSubmit($(this).find("input[type=text]").get(0));					
					},
				Cancel: function() 
					{
						$(this).dialog("close");
					}
			}
		});
		$("#renamebasket").keypress(function(e) {
		    if (e.keyCode == $.ui.keyCode.ENTER)
		    	self.doRenamebasketSubmit($(this).find("input[type=text]").get(0));
		});
			
		this.parentParentInitialize(true);
	}
}



