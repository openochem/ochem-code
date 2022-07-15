function MoleculeBrowser()
{
	var self = this;
	this.controller = "molbrowser";
	this.scope = ".browserScope";
	Browser.call(this);
	EventDispatcher.call(this);
	this.itemTemplate = "js/templates/mol-browser.ejs";
	this.itemElement = "mapping";
	this.parentActionQuery = this.getActionQuery;
	
	this.doZoom = function(link)
	{
		var img = $(link).find('IMG');
		var id = img.attr("id");
		if (img.hasClass('big'))
		{
			img.removeClass('big');
			img.attr("src","depiction.jsp?id="+id);
		} else
		{
			img.attr("src","depiction.jsp?w=300&h=300&id="+id);
			img.addClass('big');
		}
	}
	
	
	this.fragSearch = function(fragId)
	{
		if (fragId == undefined)
		{
			$("#moleculeSearch").html("<a href='javascript:void()' onclick='moleculeBrowser.fragSearch(-1); return false;'>[search by fragment]</a>");
		}	
		else
		{
			var win = openTab("Draw molecule fragment", webRoot+"molecule/edit.do?render-mode=popup&id="+fragId);
			win.callback = function(newId)
			{
				if (newId != fragId)
				{
					self.fragBlock = $("#moleculeSearch").html();
					$("#moleculeSearch").html(new View({url: 'js/templates/molfragsearch.ejs'}).render(newId));
					self.page = 1;
					self.request();
				}
			}
		}
	}
	////
	////
	////
	this.beforeToggleselect = function()
	{
		// Changing checkbox status even before sending request to server
		var newValue = (this.currentEntity.selected == "true") ? "false" : "true";
		this.currentEntity.selected = newValue;
		this.currentBlock.find('img[name="checked"]').attr('src', 
			(newValue == "true") ? "img/icons/checked.gif" : "img/icons/unchecked.gif");
		return true;
	}
	
	this.onToggleselectSuccess = function(xml)
	{
	}
	
	this.beforeAddselect = this.beforeRemoveselect = this.beforeSelectpage = function()
	{
		this.waitingDialog.show();
		return true;
	}
	
	this.doAddtagselect = function()
	{
		var tagWin = openTab("Select tag", webRoot + "tags/show.do?render-mode=popup&type=molecule");

		tagWin.callback = function(tag)
		{
			self.callAction("addtag", "", {data: "new-tag=" + tag.id});
			tagWin.closeTab();
		}
	}
	
	this.onAddtagSuccess = this.onRemovetagSuccess = function()
	{
		this.request();
	}
	
	this.doRemovetagselect = function()
	{
		var tagWin = openTab("Select tag", webRoot + "tags/show.do?render-mode=popup&type=molecule");

		tagWin.callback = function(tag)
		{
			self.callAction("removetag", "", {data: "new-tag=" + tag.id});
			tagWin.closeTab();
		}
	}
		
	
	this.onItemsLoaded = function()
	{
		$(self.scope).find('a[action="addselect"]')
			.setClass("invisible", self.pager.totalNum > 5000);
	}
	
	this.listenEvent("items_loaded", this.onItemsLoaded);
	
	this.onItemDrawn = function()
	{
		this.currentBlock.find('a[title]').tooltip({showURL: false});
	}
		
	this.parentSetPosition = this.setPosition;
	this.setPosition = function(link)
	{
		// Context menu is outside of the record context
		// So dont change position while clicking menu item link
		// Otherwise you lose focus
		if (!$(link).hasClass("yuimenuitemlabel"))
			this.parentSetPosition(link);
	}
	
	this.parentInitialize = this.initialize;
		
	this.initialize = function()
	{
		this.parentInitialize();
	
		var dlgOptions;
		
		this.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", { 
	    	width:"325px", 
			fixedcenter:true, 
			modal:true, 
			visible:false, 
			close: false
	    });
		
		this.waitingDialog.render();
	}
	
	
	this.getActionQuery = function(name)
	{
		var data = this.parentActionQuery();

		if (name == "addbasket" || name == "removebasket")
		{
			var optionIndex = self.basketSelect.select.selectedIndex;
			var option = self.basketSelect.select[optionIndex];
			
			if (option.value == -10)
				data += "&basket-name=" + self.basketDialog.getData()["basket-name"];
			else
				data += "&basket-name=" + option.text;
		}
		else if (name == "selectpage")
		{
			// Page-targeted action. Probably move this to browser.js, since it looks pretty universal / Midnighter
			$("div[rec-id]").each(function(){
				data += "&id=" + $(this).attr('rec-id');
			});
		}
		
		return data;
	}
	
	
	this.onAddselectSuccess = this.onRemoveselectSuccess = this.onSelectpageSuccess = function()
	{
		// Unfreeze interface and reload page
		this.waitingDialog.cancel();
		self.request();
	}
	
	
	this.onRemoveselectedError = this.onSelectpageError = this.onAddselectedError = function(msg)
	{
		// On error, unfreeze interface and alert this error
		this.waitingDialog.cancel();
		window.alert(msg);
	}
	
	this.onAddbasketSuccess = function()
	{
		this.fireEvent("basket_ready");
	}
	
	this.onRemovebasketSuccess = function()
	{
		this.fireEvent("basket_ready");
	}
}


