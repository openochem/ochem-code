function BasketBrowser()
{
	this.controller = "basket";
	this.url = "basket/list.do";
	this.actionURL = "basket/action.do";
	this.itemTemplate = "js/templates/basket.ejs";
	Browser.call(this);
	this.pager.selectors.pager = ".pgr";
	this.itemElement = "basket";
	this.selectedSets = new Array();
	
	
	var self = this;
	
	this.listenEvent('items_loaded', function(){
		if (window.callback)
			$("[action='check']").remove();
	});
	
	this.parentSetPosition = this.setPosition;
	this.setPosition = function(link)
	{
		// Context menu is outside of the record context
		// So dont change position while clicking menu item link
		// Otherwise you lose focus
		if (!$(link).hasClass("yuimenuitemlabel"))
			this.parentSetPosition(link);
	}
	
	this.doCheck = function(link)
	{
		var id = $(link).attr('id');
		
		var text = $('img[name="'+id+'"]').attr("src");
		if(text=="img/icons/unchecked.gif")
		{
			this.selectedSets.push(id);
			$('img[name="'+id+'"]').attr("src","img/icons/checked.gif");
		}
		else
		{
			for(var i=0; i < this.selectedSets.length; i++)
			{
				if(this.selectedSets[i] == id)
				{
					this.selectedSets.splice(i,1);
					$('img[name="'+id+'"]').attr("src","img/icons/unchecked.gif");
				}
			}
		}
		
		if(this.selectedSets.length > 1)
		{
			$("[action='combine']").removeClass("invisible");
			$("[action='compare']").removeClass("invisible");
		}
		else
		{
			$("[action='combine']").addClass("invisible");
			$("[action='compare']").addClass("invisible");
		}
		
		// window.alert(this.selectedSets);
		// sampleBrowser.selectedItems.remove or .push();
	}
	
	this.doCompare = function()
	{
		// construct URL and redirect to it.
		var mol_set = this.selectedSets.join(",");
		url = webRoot+"epbrowser/show.do?basket-select="+mol_set+"&compare-baskets=true";
		window.location.href = url;
	}
	
	this.doCombine = function()
	{
		var mol_set = this.selectedSets.join(",");
		url = webRoot+"epbrowser/show.do?basket-select="+mol_set+"&join-baskets=true";
		window.location.href = url;
	}

	this.doUpload = function()
	{
		var modelWin = openTab("Upload basket", webRoot+"basket/getbasket.do?render-mode=popup&id="+this.currentRecordId);
	}
	
	this.doIndices = function()
	{
		var modelWin = openTab("Calculate descriptors", webRoot+"modelconfigurator/choose.do?descriptors=1&basketId="+this.currentRecordId);
	}	
	
	this.doBasket = function()
	{
		openTab("Export the basket", "basket/exportBasket.do?render-mode=popup&id="+this.currentRecordId);
	}
	
	this.doModel = function()
	{
		var modelWin = openTab("Build fast model", webRoot+"basket/model.do?render-mode=popup&id="+this.currentRecordId);
	}
	
	this.doRecordmenu = function(link)
	{
		this.recordMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
		this.recordMenu.show();
	}
	
	this.doBasketmenu = function(link)
	{
		this.basketMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
		this.basketMenu.show();
	}
	
	this.beforeAdd = function()
	{
		var newName = prompt("Enter new basket name", "");
		if (newName)
		{
			this.filters.setValue("name", newName);
			return true;
		}
		return false;
	}
	
	this.beforeDelete = function()
	{
		return window.confirm("You are going to delete the basket: \"" + this.currentEntity.name + "\"\nAre you sure you want to continue?");
	}
	
	this.beforeRename = function()
	{
		var newName = prompt("Enter new basket name", "");
		if (newName)
		{
			this.filters.setValue("newname", newName);
			return true;
		}
		return false;
	}
	
	this.onDeleteSuccess = this.onRenameSuccess = function() {this.request()};
	
	this.doNameclick = function()
	{
		if (window.callback)
			this.callAction("select");
	}

}

include.plugins('view');
var sampleBrowser = new BasketBrowser();
$(document).ready(function() {
	sampleBrowser.initialize();
	sampleBrowser.recordMenu = new YAHOO.widget.Menu("basicxlsmenu");
	sampleBrowser.recordMenu.render(); 
	sampleBrowser.basketMenu = new YAHOO.widget.Menu("basicbasketmenu");
	sampleBrowser.basketMenu.render(); 
});