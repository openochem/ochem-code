function BatchBrowser()
{
	this.controller = "batchupload";
	Browser.call(this);
	this.itemElement = "exp-property";	
	this.itemTemplate = "js/templates/batchupload-record.ejs";
	var self = this;
	//self.filters.values["show"] =  "all";
	this.parentActionQuery = this.getActionQuery;
	this.pager.selectors.pager = ".pgr";
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('.block-image IMG').click(function(){
			$(this).toggleClass('big');
		});
		
		this.currentBlock.find('a').tooltip({showURL: false});
		
		var recStatus = this.currentEntity.ep_status;
		
		if(recStatus =="0")
			this.currentBlock.addClass("highlighterror");
		if(recStatus =="1")
			this.currentBlock.addClass("highlightnon-v");
	}
	
	this.filters.onFilterChange = function()
	{
	}
	
	this.doButton = function(link) 
	{
		//self.filters.values = new Object();
		self.filters.setValue("show", $(link).attr("value"));
		self.request();
	}
	
	this.listenEvent('items_loaded', function(){
		$('.openable h1').click(function(){
			$(this).parent('.openable').toggleClass('opened');
		});
	});
	
	this.doToggle = function(link)
	{
		var showtext = $(link).attr("showtext");
		var hidetexttext = $(link).attr("hidetext");
		var blocknum = $(link).attr("blocknum");
		
		var table = $(link).parents("div[rec-id]").find(".block-table").get(blocknum);
		$(table).toggleClass("invisible");
		if ($(table).hasClass("invisible"))
			$(link).html(showtext);
		else
			$(link).html(hidetext);
	}
	
	this.getActionQuery = function(action)
	{
		query = self.parentActionQuery();
		if (action=="state")
		{
			var radioButton = this.currentBlock.find("[action='state']").each(function()
			{
				if ($(this).attr("checked") == true) //doesn't work as a normal selector
				{
					var value =  ($(this).attr('value') == undefined) ? "" : URLEncode($(this).attr('value'));
					query += "&" + $(this).attr('name') + "=" + value;
				}
			});
		}
		return query;
	}
	
	this.onStateSuccess = function()
	{
	}
	
	this.onStateoverSuccess = this.onStatenoSuccess = this.onStatenewSuccess = function()
	{
		self.request();
	}
	
}

