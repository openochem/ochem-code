
function PropertyForm()
{
	this.scope = ".formscope";
	var self = this;
	EditForm.call(this);
	
	this.actionURL = "eprecord/action.do";
	this.unitsSelect = new UnitSelect('n-unit', 'unit/list.do?lightweight=1', 'unit');
	this.optionsSelector = new OptionsSelector("n-option");
	//this.optionsSelect = new DynamicSelect('n-option', 'properties/listoptions.do', 'option');
	this.itemElement = "exp-property";
	
	this.doEditmolecule = function()
	{
		var molWin = openTab("Edit molecule", webRoot + "molecule/edit.do?render-mode=popup&id="+propertyForm.getValue('n-molecule'));
		molWin.moveBy(300, 50);
		molWin.callback = function(newId)
		{
			var count = propertyForm.getValue('count-similar-ep');
			if(propertyForm.getValue('n-molecule') != newId && count > 1)
			{
				$("#modify").removeClass("invisible");
			}
			propertyForm.setValue('n-molecule', newId);
			$('img#depiction').attr('src', 'depiction.jsp?id=' + newId);
			propertyForm.setValue('mol-id', newId);
			molWin.closeTab();
				
			setTimeout('synonymBrowser.request()', 1);
			setTimeout('nameBrowser.request()', 1);
		}
	}
	
	this.doEditproperty = function()
	{
		var propWin = openTab("Select property", webRoot + "properties/show.do?render-mode=popup&selected="+this.getValue("n-property"));

		propWin.callback = function(newProperty)
		{
			//propWin.document.write('callback!!!');
			propertyForm.property = newProperty;
			propertyForm.setValue('n-property', newProperty["id"]);
			$('[name="property-link"]').html(newProperty["name"]);
			propWin.closeTab();
			setTimeout('propertyForm.unitsSelect.update("category='+newProperty.unitCategory["id"]+'")', 1);
			setTimeout('propertyForm.onPropertyChanged()', 1);
			$("#Property-info").removeClass("invisible");
		}
	}
	
	this.unitsSelect.listenEvent("update", function(){
		if (propertyForm.property)
			$(propertyForm.unitsSelect.select).val(propertyForm.property.defaultUnit.id);
	});
	
	this.onPropertyChanged = function()
	{
		if (this.property.qualitive == "true")
		{
			$("#Qualitive").removeClass("invisible");
			$("#Quantitive").addClass("invisible");
		}
		else
		{
			$("#Quantitive").removeClass("invisible");
			$("#Qualitive").addClass("invisible");
		}
		this.optionsSelector.update(this.property);
	}
	
	this.doEditarticle = function()
	{
		var articleWin = openTab("Select article", webRoot + "article/show.do?render-mode=popup&selected="+this.getValue("n-article"));
		articleWin.callback = function(newArticle)
		{
			propertyForm.setValue('n-article', newArticle["id"]);
			$('[name="article-link"]').html(newArticle.title);
			articleWin.closeTab();
		}
	}
	
	this.doEditreference = function()
	{
		var appendix = (this.getValue("id") > 0) ? 
			"similarto="+this.getValue('id') :
			"similarmol="+this.getValue('n-molecule')+"&property="+this.getValue('n-property');
		var win = openTab("Select reference", webRoot + "epbrowser/show.do?render-mode=popup&experimental=1&render-mode=popup&"+appendix+"&selected="+this.getValue('n-reference'));
		win.moveBy(100, 100);
		win.callback = function(newValue)
		{
			$("a[name='reference-link']").html(newValue.article.title);
			self.setValue('n-reference', newValue.id);
			win.closeTab();
		}	
	}
	
	this.doDescribeproperty = function()
	{
		openTab("Property description", "properties/edit.do?id="+this.getValue("n-property"));
	}
	
//	this.doWikiproperty = function()
//	{
//		openTab("Property description", "wikipage/action.do?entities=property&id="+this.getValue("n-property"));
//	}
	
	this.doViewrecord = function()
	{
		var simiWin = openTab("View similar records", webRoot + "epbrowser/action.do?render-mode=popup&action=similar-record&id="+this.getValue('id'));
		simiWin.moveBy(50, 50);
	}
	
	this.beforeEdit = function()
	{
		this.waitingDialog.show();
		return true;
	}
	
	this.onItemSaved = function(entity)
	{
		conditionsBrowser.filters.setValue('id', entity["exp-property"].id);
		nameBrowser.filters.setValue('id', entity["exp-property"].id);
		this.changed = false;
		nameBrowser.callAction("saveall");
		
		//window.callback(entity);
	}
	
	this.listenEvent('formchanged', function()
	{
		document.title = document.title + " (modified)";
	});
	
	window.onbeforeunload = function(){
		if (self.changed)
			return "The data was changed. If you press OK changes will be lost";
	} 
	
	$(document).ready(function()
	{
		$("select[name='evidence']").change(function(){
			if (self.getValue("evidence") == 2)
				$("tr[name='reference']").removeClass("invisible");
			else
				$("tr[name='reference']").addClass("invisible");
		});	
		var slPredicate = $("select[name='n-predicate']");
		slPredicate.change(function(){
			var option = $(this.options[this.selectedIndex]);
			if (option.attr("shortname") == "+-" || option.attr("shortname") == "-")
			{
				$("#second-value").removeClass("invisible");
				$("#predicate").html(option.attr("name"));
			}
			else
			{
				$("#second-value").addClass("invisible");
			}
		});
		
		slPredicate.val(slPredicate.attr("selected"));
		
		$("select[name='evidence']").change();
		slPredicate.change();
		
		
		if (propertyForm.getValue("initial-category") != "")
		{
			propertyForm.property = {"defaultUnit": {"id": propertyForm.getValue("initial-unit")}};
			propertyForm.unitsSelect.update("category="+propertyForm.getValue("initial-category"));
		}
	});
	
	this.initialize = function()
	{
		if (this.getValue("options-count") > 0)
			this.optionsSelector.update(this.getValue("n-property"), this.getValue("original-option"), this.getValue("original-option-name"));
		
		this.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", { 
	    	width:"325px", 
			fixedcenter:true, 
			modal:true, 
			visible:false, 
			close: false
	    });
		this.waitingDialog.render();
		
		this.recordDialog = new YAHOO.widget.Dialog("recordDialog", { 
	    	width:"525px", 
			fixedcenter:true, 
			modal:true, 
			visible:false 
	    });
	    
	    var myButtons = [{ text:"Show record", handler:handleSubmit}, 
						{ text:"Cancel", handler:handleCancel } ]; 
		this.recordDialog.cfg.queueProperty("buttons", myButtons);
		this.recordDialog.render();
	}
	
	var handleSubmit = function() { 
		var dupRec = $("input[name='rec-id']").val();
		var recordWin = openTab("Duplicate record", webRoot + "epbrowser/show.do?render-mode=popup&id="+dupRec);
		this.cancel();
	};
	
	var handleCancel = function() { 
		this.cancel();
	};
	
	this.onEditError = function(msg, entity){
		// On error, unfreeze interface and alert this error
		this.waitingDialog.cancel();
		
		if (entity.attachment)
		{
			$("input[name='rec-id']").val(entity.attachment.content);
			$("#dupRecords").html(msg);
			this.recordDialog.show();
		}
		else
			window.alert(msg);
	}
	
	//this.ajax.showError = this.ajax.showErrorExtended;
}

function ConditionValuesBrowser()
{
	this.controller = "properties";
	this.scope = ".conditionsscope";
	UnifiedBrowser.call(this);
	this.itemElement = "property-value";
	this.itemTemplate = "js/templates/conditionvalue.ejs";
	this.filters.scope = this.scope;
	this.url = "properties/listvalues.do";
	this.actionURL = "properties/savevalues.do";
	this.container = "ConditionsBrowser";
	this.minusCounter = -1;
	
	var self = this;
	
	this.doEditcondition = function()
	{
		var propWin = openTab("Select condition", webRoot + "properties/show.do?render-mode=popup&condition=true");
		propWin.moveBy(50, 200);
		propWin.callback = function(newProperty)
		{
			self.currentEntity.property = newProperty;
			self.setValue('cond-id', newProperty["id"]);
			self.currentBlock.find('[name="condition-link"]').html(newProperty.name);
			propWin.closeTab();
			
			if (newProperty.qualitive == "true")
				// Update options dropbox
				setTimeout('conditionsBrowser.currentBlock.get(0).optionsSelect.update({id: "'+newProperty.id+'", "options-count": "'+newProperty["options-count"]+'"})', 1);
			else
			{
				// Update units dropbox
				setTimeout('conditionsBrowser.currentBlock.get(0).unitsSelect.update("category='+newProperty.unitCategory.id+'&lightweight=1")', 1);
				setTimeout('conditionsBrowser.currentBlock.get(0).predSelect.update()', 1);
			}
			self.updateVisibility();				
		}
	}
	
	this.onItemDrawn = function()
	{
		var dom = this.currentBlock.get(0);
		
		// Create dynamic select for units of condition
		var dynSelect = new UnitSelect('cond-unit', 'unit/list.do?lightweight=1', 'unit', this.currentBlock);
		dom.unitsSelect = dynSelect;
		
		// Create dynamic select for options for this condition ("qualitive conditions")
		var optSelect = new OptionsSelector('cond-option', this.currentBlock);
		console.log(this.currentEntity);
		if (this.currentEntity.option)
			optSelect.update(this.currentEntity.property, this.currentEntity.option.id, this.currentEntity.option.name, true);
		dom.optionsSelect = optSelect;
		
		// Dynamic select for predicates...
		var predSelect = new DynamicSelect('cond-pred', 'properties/listpredicates.do', 'predicate', this.currentBlock);
		dom.predSelect = predSelect;
		
		dom.unitsSelect.listenEvent("update", function()
		{
			$(conditionsBrowser.currentBlock.get(0).unitsSelect.select).val(conditionsBrowser.currentEntity.property.defaultUnit.name);
		});
		
		// Show/hide appropriate sections
		this.updateVisibility();
		
		var currentBlock = this.currentBlock;
		var slPredicate = currentBlock.find("select[name='cond-pred']");
		slPredicate.change(function(){
			var option = this.options[this.selectedIndex];
			if (option.value == 8 || option.value == 9) //Special - and +- predicates
			{
				currentBlock.find("span[name='second-value']").removeClass("invisible");
				currentBlock.find("span[name='predicate']").removeClass("invisible");
				currentBlock.find("span[name='predicate']").html(option.text);
			}
			else
			{
				currentBlock.find("span[name='second-value']").addClass("invisible");
			}
		});
		
		slPredicate.val(slPredicate.attr("selected"));
		
		$("select[name='evidence']").change();
		slPredicate.change();
	}
	
	this.updateVisibility = function()
	{
		this.currentBlock.find('span[name]').addClass("invisible");
		this.currentBlock.find('span[name="type-' + this.currentEntity.property.type + '"]').removeClass("invisible");
	}
	
	this.doNew = function()
	{
		var entity = {"value":0.0, "property": {"id": this.minusCounter--, "name":"[...]", "unitCategory":{"unit": new Array()}, "predicates":{"predicate": new Array()}}, "unit": {"name": "Click", "id": -1}, "textualValue":""};
		this.drawFromJSON(entity);
	}
	
	this.doDelete = function()
	{
		this.deleteRecord();
	}
	
	this.getRecordIdentifier = function(entity)
	{
		return entity.property.id;
	}
	
	this.onSaveallSuccess = function(xml)
	{
		// After conditions are saved, close the window
		window.callback(propertyForm.entity);
	}
}


function SynonymBrowser()
{
	var self = this;
	
	this.controller = "eprecord";
	this.scope  = ".synonymscope";
	this.useHistory = false;
	
	UnifiedBrowser.call(this);
	
	this.pager = new DirectPagerLite(11, true);
	this.pager.selectors.scope = this.scope;

	this.pager.onStateChanged = function() 
	{
		self.filters.values.pagesize = self.pager.itemsPerPage;
		self.page = self.pager.currentPage;
		self.request();
	};

	
	this.filters.scope = this.scope;
	this.itemElement = "moleculename";
	this.container = "SynonymBrowser";
	this.itemTemplate = "js/templates/molsynonym.ejs";
	this.url = "eprecord/getsynonyms.do";
	this.filters.firstTime = false;
	this.minusCounter = -1;
}

function NameBrowser()
{
	this.controller = "eprecord";
	this.scope = ".namescope";
	
	UnifiedBrowser.call(this);
	
	this.itemElement = "moleculename";
	this.itemTemplate = "js/templates/molname.ejs";
	this.filters.scope = this.scope;
	this.url = "eprecord/listnames.do";
	this.actionURL = "eprecord/nameactions.do";
	this.container = "NameBrowser";
	this.filters.firstTime = false;
	this.minusCounter = -1;
	
	var self = this;
	
	this.entityCache = new Array();
	
	this.listenEvent('items_load', function()
	{
		$(self.scope).find("div[rec-id^=-]").each(
			function (i) {
				self.currentBlock = $(this);
				this.entity.name = self.getValue('moleculename').trim();
				self.entityCache[i] = this.entity;
			}	
      	);
	});

	this.listenEvent('items_loaded', function()
	{
		var i = 0;
		for (i = 0; i<self.entityCache.length; i++)
		{
			self.drawFromJSON(self.entityCache[i]);
		}
		this.entityCache= new Array();
	});
	
	this.drawFromArray = function()
	{	
		if (self.items)
			if (self.items instanceof Array){
				for (var i = 0; i < self.items.length; i++){
					self.drawFromJSON(self.items[i]);
				}
			}
			else
				self.drawFromJSON(self.items);					
	}
	
	this.doNew = function()
	{
		this.minusCounter--;
		var entity = {"id": this.minusCounter, "name": ""};
		this.drawFromJSON(entity);
	}
	
	
	this.doDelete = function()
	{
		this.deleteRecord();
	}
	
	this.getRecordIdentifier = function(entity)
	{
		return entity["id"];
	}
	
	this.onSaveallSuccess = function()
	{
		//conditionsBrowser.callAction("saveall"); / Conditions are now saved simultaneously with property, in the same request
		window.callback(propertyForm.entity);
	}
	
	this.doSearch = function()
	{
		var searchstring = "";
		if (undefined != this.currentBlock) 
			searchstring = escape(this.getValue('moleculename').trim());
		
		var propWin = openTab("NCBI Search", webRoot + "ncbisearch/searchBrowser.do?render-mode=popup&search-string=" + searchstring
							+ "&existmol-id=" + propertyForm.getValue('n-molecule')
							+ "&mapping1-id=" + propertyForm.getValue('n-molecule-mapping1')
							+ "&mapping2-id=" + propertyForm.getValue('n-molecule-mapping2'));
		
		propWin.moveBy(50, 50);
		propWin.callback = function(molecule)
		{

			self.setValue('moleculename', molecule.searchedBy.name);
			// maybe to redraw the item here to see the new name, if it was an saved name
			//window.alert(molecule.id);
			var count = propertyForm.getValue('count-similar-ep');
			if((propertyForm.getValue('n-molecule') != molecule.id) && (count > 1))
			{
				$("#modify").removeClass("invisible");
			}
			propertyForm.setValue('n-molecule', molecule.id);
			$('img#depiction').attr('src', 'depiction.jsp?id=' + molecule.id);
			propertyForm.setValue('mol-id', molecule.id);
			propWin.closeTab();
			
			
			setTimeout('synonymBrowser.request()', 1);
			setTimeout('nameBrowser.request()', 1);
		}
	}
	
	this.onValidateSuccess = function()
	{
		self.request();
	}
	
	this.onInvalidateSuccess = function()
	{
		self.request();
	}
	
	this.beforeChecknames = function()
	{
		$("#NameBrowser").html('<img src="img/long_green.gif"/>');
		return true;
	}
	
	this.onChecknamesSuccess = function()
	{
		self.request();
	}
	
}

/**
 * The UI element for selecting a property option.
 * Works in two modes:
 * 	1) Using a simple drop-down for small number of options
 *  2) Using a browser of options for large number of options
 */
OptionsSelector = function(selectName, scope)
{
	var scope = scope || $(document);
	var self = this;
	this.dynamicSelect = new DynamicSelect(selectName, 'propertyoptions/list.do', 'option', scope);
	
	this.update = function(propertyID, currentValue, currentTitle, onlyIfRequired)
	{
		
		if (typeof propertyID === 'object')
		{
			var property = propertyID;
			var select = scope.find("select[name="+selectName+"]");
			var span = select.next("span");
			if (property["options-count"] < 50)
			{
				// Select box
				select.attr("send", 1);
				span.html("");
				
				onlyIfRequired = onlyIfRequired || false;
				if (!onlyIfRequired)
				{
					self.dynamicSelect.update("property=" + property.id, currentValue);
					select.removeClass("invisible");
					span.addClass("invisible");
				}
			}
			else
			{
				// Browser link
				select.removeAttr("send");
				span.removeClass("invisible");
				select.addClass("invisible");
				var html = $('<a href="#" title="Click to choose the property class">[...]</a><input type="hidden" name="' + selectName + '" value="' + currentValue + '" send="1"/>"');
				span.html(html);
				if (currentTitle)
					span.find("a").html(currentTitle);
				$("select[name="+selectName+"]").next("span").find("a").click(function(){
					var win = openTab("Select an option", webRoot + "propertyoptions/show.do?property=" + property.id);
					win.callback = function(option)
					{
						span.find("a").html(option.name);
						span.find("input").val(option.id);
						win.closeTab();
					}
					return false;
				});
			}
		}
		else
		{
			// Fetch the property fom DB by ID
			var ajax = new QSPR.Ajax();
			ajax.send({
				url: "properties/edit.do?id=" + propertyID + "&out=json",
				success: function(model)
				{
					self.update(model.property, currentValue, currentTitle);
				}
			});
		}
	}
}

include.plugins('view');
var propertyForm = new PropertyForm();
var conditionsBrowser = new ConditionValuesBrowser();
var synonymBrowser = new SynonymBrowser();
var nameBrowser = new NameBrowser();

$(document).ready(function()
{
	conditionsBrowser.initialize(); 
	synonymBrowser.initialize();
	nameBrowser.initialize(true);
	nameBrowser.drawFromArray();
	propertyForm.initialize();	
});
