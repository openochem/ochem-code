function FragmentBrowser()
{
	this.controller = "fragmentbrowser";
	//this.url = "fragmentbrowser/list.do";
	//this.actionURL = "fragmentbrowser/action.do";
	this.itemTemplate = "js/templates/fragment.ejs";
	Browser.call(this);
	this.pager.selectors.pager = ".pgr";
	this.itemElement = "mapping1Fragment";
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
	
	this.doZoom = function(link)
	{
		var img = $(link).find('IMG');
		var id = img.attr("id");
		if (img.hasClass('big'))
		{
			img.removeClass('big');
			img.attr("src","depiction.jsp?w=150&h=150&frag_id="+id);
		} else
		{
			img.attr("src","depiction.jsp?w=300&h=300&frag_id="+id);
			img.addClass('big');
		}
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
		url = webRoot+"epbrowser/show.do?fragment-select="+mol_set+"&compare-fragments=true";
		window.location.href = url;
	}
	
	this.doCombine = function()
	{
		var mol_set = this.selectedSets.join(",");
		url = webRoot+"epbrowser/show.do?fragment-select="+mol_set+"&join-fragments=true";
		window.location.href = url;
	}

	this.doUpload = function()
	{
		var modelWin = openTab("Upload fragment", webRoot+"fragment/getfragment.do?render-mode=popup&id="+this.currentRecordId);
	}
	
	this.doIndices = function()
	{
		var modelWin = openTab("Calculate e-state indices", webRoot+"fragment/indices.do?render-mode=popup&id="+this.currentRecordId);
	}	
	
	this.doFragment = function()
	{
		window.location = webRoot + "fragment/exportFragment.do?render-mode=popup&id="+this.currentRecordId;
	}
	
	this.doModel = function()
	{
		var modelWin = openTab("Build fast model", webRoot+"fragment/model.do?render-mode=popup&id="+this.currentRecordId);
	}
	
	this.doRecordmenu = function(link)
	{
		this.recordMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
		this.recordMenu.show();
	}
	
	this.doFragmentmenu = function(link)
	{
		this.fragmentMenu.cfg.setProperty("context", [$(link).attr('id'), "tl", "tr"]); 
		this.fragmentMenu.show();
	}
	
	this.beforeAdd = function()
	{
		var newName = prompt("Enter new fragment name", "");
		if (newName)
		{
			this.filters.setValue("name", newName);
			return true;
		}
		return false;
	}
	
	this.beforeRename = function()
	{
		var newName = prompt("Enter new fragment name", "");
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
var sampleBrowser = new FragmentBrowser();
$(document).ready(function() {
	sampleBrowser.initialize();
	//sampleBrowser.recordMenu = new YAHOO.widget.Menu("basicxlsmenu");
	//sampleBrowser.recordMenu.render(); 
	//sampleBrowser.fragmentMenu = new YAHOO.widget.Menu("basicfragmentmenu");
	//sampleBrowser.fragmentMenu.render(); 
	$("#original-molecule").attr('src', 'depiction.jsp?w=150&h=150&id=' + getParams["id"]);
});