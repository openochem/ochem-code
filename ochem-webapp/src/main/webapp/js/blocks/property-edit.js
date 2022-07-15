include.plugins("view");

var propertyChanged = false;

$(document).ready(function() 
{
	var select = $('select[name="category"]');
	var typeSelect = $("select[name='type']");
	select.change(function() {
		dynSelect.update("category="+ajaxForm.getValue('category'), ajaxForm.getValue("initial-unit"));
	});
	if ($('select[name="unit"]').get(0).options.length == 0)
		select.change();
	typeSelect.change(
		function(){
			$("#Options").setClass("invisible", this.value != 1);
			$("#Units").setClass("invisible", this.value != 0);
	});
	
	if (ajaxForm.getValue("id"))
	{
		$("[onlynew]").attr("disabled", "disabled");
	}
	
	typeSelect.val(typeSelect.attr("selected-value"));
	typeSelect.change();
	if (ajaxForm.embeddedOptionsBrowser())
		optionBrowser.initialize();
	conditionsBrowser.initialize(true);
	
	if (window.openerTab != undefined && window.openerTab.sampleBrowser != undefined && window.openerTab.sampleBrowser.filters.getValue('condition'))
	{
		// Its a condition
		ajaxForm.fields["condition"] = "true";
		$("#Conditions").addClass("invisible");
		$("#MiscConditions").addClass("invisible");
	}
	else
		$("select[name='type'] option[value=2]").remove();
	
	if (getParams["directory"])
		ajaxForm.fields["is-directory"] = "true";
	
	ajaxForm.confirmDialog = new YAHOO.widget.Dialog("confirm-dialog", { 
    	width: "500px", 
		fixedcenter: true, 
		modal: true, 
		visible: false 
    });
	
	ajaxForm.confirmDialog.cfg.queueProperty("buttons", 
			[{ 
				text:"Confirm changes", handler: function() 
				{
					ajaxForm.callAction("edit", null, {data: "confirmed=true"});
					this.cancel();
				} 
			}, 
	         { 
				text:"Cancel", handler: function()
				{
					this.cancel();
	         	} 
	         }]);
	
	ajaxForm.confirmDialog.render();
	
	ajaxForm.waitingDialog = new YAHOO.widget.Dialog("waiting-dialog", { 
    	width:"325px", 
		fixedcenter:true, 
		modal:true, 
		visible:false, 
		close: false
    });
	
	ajaxForm.waitingDialog.render();
});


var ajaxForm = new EditForm();

ajaxForm.scope = ".EditForm";
ajaxForm.actionURL = "properties/action.do";
ajaxForm.itemElement = 'property';
ajaxForm.restriction = "form";

ajaxForm.ajax.beforeRequest = function()
{
	this.waitingDialog.show();
}

ajaxForm.ajax.afterRequest = function()
{
	this.waitingDialog.cancel();
}

ajaxForm.isQualitive = function()
{
	return this.getValue("type") == 1;
}

ajaxForm.embeddedOptionsBrowser = function()
{
	return this.getValue("options-count") < 100;
}

ajaxForm.beforeEdit = function()
{
	if(ajaxForm.getValue('name').length < 3 || ajaxForm.getValue('name').length > 100)
	{
		window.alert("Property name is too short or too long and can not be added!");
		return false;
	}
	if(ajaxForm.getValue('description').trim().length < 50)
	{
		window.alert("Property description is too short and can not be added");
		return false;
	}
	return true;
}

ajaxForm.onEditSuccessOverrided = ajaxForm.onEditSuccess;
ajaxForm.onEditSuccess = function(entity)
{
	this.entity = entity;
	optionBrowser.filters.setValue('id', entity.property.id);
	tagsBrowser.filters.setValue('id', entity.property.id);
	tagsBrowser.onSaveallSuccess = function()
	{
		if (ajaxForm.embeddedOptionsBrowser())
			optionBrowser.callAction("saveall");
		else
			ajaxForm.onEditSuccessOverrided(entity);
	}
	tagsBrowser.callAction("saveall");
}

ajaxForm.onEditError = function(message)
{
	if (message.startsWith("Dear user"))
	{
		$("#confirm-message").html(message.replace(/\n/g, "<br/>"));
		this.confirmDialog.show();
	}
	else
		window.alert(message);
}

ajaxForm.doSelectparent = function()
{
	var selectWin = openTab("Select a parent property", webRoot + "properties/show.do?render-mode=popup&directories=true");
	selectWin.callback = function(selectedItem)
	{  
		selectWin.closeTab();
	}
}

function OptionBrowser()
{
	this.controller = "properties";
    this.container = "Browser";
	UnifiedBrowser.call(this);
    this.filters.useUrlParameters = false;
	this.url = "properties/listoptions.do"; // Do not load items through AJAX
	this.actionURL = "properties/saveoptions.do";
	this.scope = "#Options";
	this.itemElement = "option";
	this.itemTemplate = "js/templates/conditionoption.ejs";
	this.minusId = -1;
	this.changed = false;
	
	this.doDelete = function()
	{
		propertyChanged = true;
		this.deleteRecord();
	}
	
	this.doEditoption = function()
	{
		this.currentBlock.find(".option-name").addClass("invisible");
		this.currentBlock.find("input[name='co-name']").remove();
		this.currentBlock.find("a[action=editoption]").remove();
		this.currentBlock.find("input[type=text]").removeClass("invisible").attr("name", "co-name");
		this.currentBlock.find('input[name="co-id"]').val(--this.minusId); // When editing, delete the record and recreate it
		this.changed = true;
	}
	
	this.onSaveallSuccess = function()
	{
		window.callback(ajaxForm.entity);
	}
	
	
	this.doAdd = function()
	{
		var entity = {name: "unnamed", id: this.minusId--};
		this.drawFromJSON(entity, {method: 'prepend'});
	}
}

var ConditionsBrowser = function()
{
	var self = this;
	UnifiedBrowser.call(this);
	this.container = "ConditionsBrowser";
	this.scope = "#Conditions";
	this.itemTemplate = "js/templates/condition-inline.ejs";
	this.minusId = 0;
	
	this.onItemDrawn = function()
	{
		this.currentBlock.css({float: "left", width: "auto", padding: "0px 5px 0px 0px"});
		this.mainDiv.find('span').remove();
	}
	
	this.doDelete = function()
	{
		this.deleteRecord();
		propertyChanged = true;
	}
	
	this.doAdd = function()
	{
		//this.drawFromJSON({id: --this.minusId, name: "[...]"});
		var selectWin = openTab("Select a condition", webRoot + "properties/show.do?render-mode=popup&condition=true");
		selectWin.callback = function(selectedItem)
		{  
			self.drawFromJSON({id: selectedItem.id, name: selectedItem.name});
			propertyChanged = true;
			selectWin.closeTab();
		}
	}
	
	this.listenEvent("items_loaded", function() {
		if (self.pager.totalNum == 0)
			self.mainDiv.html("<span>No tags selected</span>");
	})
	
	this.doConditionclick = function()
	{
		if (this.getValue('condition-id') < 0)
		{
			var selectWin = openTab("Select a condition", webRoot + "properties/show.do?render-mode=popup&condition=true");
			selectWin.callback = function(selectedItem)
			{
				self.setValue('condition-id', selectedItem.id);
				self.currentBlock.find('a[name="condition-link"]').html(selectedItem.name);
				selectWin.closeTab();
			}
		}
//		else
//			openTab("Wiki", webRoot + 'wikipage/action.do?entities=property&id='+this.currentEntity.id);
	}
}

var ConditionsActionable = function()
{
	var self = this;
	Actionable.call(this);
	this.scope = "#MiscConditions";
	
	this.doConditionclick = function(link)
	{
		openTab("Condition profile", webRoot + 'properties/edit.do?id='+$(link).attr("id"));
	}
	
	this.doCountclick = function(link)
	{
		openTab("Compounds browser", webRoot + 'epbrowser/show.do?property='+$(link).attr("propid")+'&cond-id-1='+$(link).attr("id"));
	}
	
}

var dynSelect = new UnitSelect('unit', 'unit/list.do?lightweight=1', 'unit');
var optionBrowser = new OptionBrowser();
var conditionsBrowser = new ConditionsBrowser();
var conditionsActionable = new ConditionsActionable();
tagsBrowser.scope = "#Tags";

function onTabClose()
{
	if (propertyChanged || tagsBrowser.changed)
		return confirm("There are unsaved changes in the property. To save the changes, please use SAVE button.\n\nIf you close the tab, your changes will be lost. \nDo you want to close the tab?");
	return true;
}

function checkPropName(obj)
{
	if(obj.value.length <= 1)
	{
		$("#propName").css("color","red").html("error: property name length is "+obj.value.length+" and can not be saved (min. 2 characters)");
	}else
	if (obj.value.length > 1 && obj.value.length < 5)
	{
		$("#propName").css("color","#FA58F4").html("warning: property name length is "+obj.value.length+" and may not be informative");
	}
	else if(obj.value.length >= 5 && obj.value.length < 40)
	{
		$("#propName").css("color","green").html("property name is good and length is "+obj.value.length);
	}
	else if(obj.value.length >= 40)
	{
		$("#propName").css("color","red").html("error: property name length is "+obj.value.length+", its should be less than 40 characters and can not be saved (max. 40 characters)");
	}
}

function checkPropDesc(obj)
{
	if(obj.value.length <= 50)
	{
		$("#propDesc").css("color","red").html("error: property description length is "+obj.value.length+" and can not be saved (min. 50 characters)");
	}else
	if(obj.value.length > 50 && obj.value.length < 200)
	{
		$("#propDesc").css("color","#FA58F4").html("warning: property description length is "+obj.value.length+" and may not be informative");
	}
	else if(obj.value.length > 200)
	{
		$("#propDesc").css("color","green").html("property description is good");
	}
}