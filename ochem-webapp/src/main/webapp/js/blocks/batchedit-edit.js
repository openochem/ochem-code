function BatcheditEditForm()
{
	this.scope = ".formscope";
	var self = this;
	EditForm.call(this);
	
	this.actionURL = "batchedit/action.do";
	//this.propunitsSelect = new DynamicSelect('n-newunit', 'unit/list.do', 'unit');
	
	this.propunitsSelect = new UnitSelect('newunit', 'unit/list.do?lightweight=1', 'unit');
	this.optionsSelect = new DynamicSelect('newoption', 'properties/listoptions.do', 'option');
	
	this.itemElement = "exp-property";
	
	this.doEditproperty = function()
	{
		var propWin = openTab("Select property", webRoot + "properties/show.do?render-mode=popup&selected=" + this.getValue("newproperty"));
		
		propWin.callback = function(newProperty)
		{
			batcheditForm.property = newProperty;
			batcheditForm.setValue('newproperty', newProperty.id);
			$('[name="property-link"]').html(newProperty.name);
			propWin.closeTab();
			//alert("new " + newProperty.id + " - " + newProperty["id"] + " # " + batcheditForm.getValue("n-category") + " - " + newProperty.unitCategory.id);
			
//			this.waitingDialog.show();
			
			//var selectDefault = newProperty.unitCategory.id != batcheditForm.getValue("n-category");
//			if (selectDefault){
			setTimeout('batcheditForm.propunitsSelect.update("category=' + newProperty.unitCategory.id + '")', 1);
//			}
//			setTimeout('batcheditForm.onPropertyChanged("' + newProperty.unitCategory.id != batcheditForm.getValue("n-category") + '")', 1);
//			setTimeout('batcheditForm.onPropertyChanged(' + selectDefault + ')', 1);
			setTimeout('batcheditForm.onPropertyChanged()', 1);
		}
	}
	
	this.propunitsSelect.listenEvent("update", function()
	{
		// property has not changed or property has the same unit category --> keep the existing unit
		if ( ! batcheditForm.property || batcheditForm.getValue("unitcategory") == batcheditForm.property.unitCategory.id) {
			var elSel = self.propunitsSelect.select;
			var elOptOld = elSel.options[0];
			
			//var elOptNew1 = new Option("do not change -- keep: \"" + self.getValue("n-unitname") + "\"", 0, false, false);
			var elOptNew1 = new Option("do not change", 0, false, false);
			try
			{
				elSel.add(elOptNew1, elOptOld);
			}
			catch(ex) 
		    {
				elSel.add(elOptNew1, 0); // IE only
		    }
			elSel.selectedIndex = 0;
		}
		
		// property has changed and the category is not the same anymore --> select default unit
		if (batcheditForm.property && batcheditForm.property.unitCategory.id != batcheditForm.getValue("unitcategory")) {
			$(batcheditForm.propunitsSelect.select).val(batcheditForm.property.defaultUnit.name);
			batcheditForm.setValue("unitcategory", batcheditForm.property.unitCategory.id);
			batcheditForm.setValue("newunit", propertyForm.property.defaultUnit.id);
		}
		
		//this.waitingDialog.cancel();
		//if (batcheditForm.property)
		//	$(batcheditForm.propunitsSelect.select).val(batcheditForm.property.defaultUnit.id);
	});
	
	this.onPropertyChanged = function()
	{
		if (this.property.qualitive == "true")
		{
			$("#Qualitive").removeClass("invisible");
			$("#Quantitive").addClass("invisible");
			
//			if (selectDefault)
//				this.optionsSelect.update("id=" + this.property.id);
		}
		else
		{
			$("#Quantitive").removeClass("invisible");
			$("#Qualitive").addClass("invisible");
			
		}
		
//		if (selectDefault) {
			this.optionsSelect.update("id=" + this.property.id);
//			$(batcheditForm.propunitsSelect.select).val(batcheditForm.property.defaultUnit.name);
//		}
		
	}
	
	this.doEditarticle = function()
	{
		var articleWin = openTab("Select article", webRoot + "article/show.do?render-mode=popup&selected="+this.getValue("newarticle"));
		articleWin.callback = function(newArticle)
		{
			batcheditForm.setValue('newarticle', newArticle["id"]);
			$('[name="article-link"]').html(newArticle.title);
			articleWin.closeTab();
		}
	}
	
	this.beforeEdit = function()
	{
		this.waitingDialog.show();
		return true;
	}
	
	this.onEditError = function(msg)
	{
		// On error, unfreeze interface and alert this error
		this.waitingDialog.cancel();
		window.alert(msg);
	}
	
	this.onItemSaved = function(entity)
	{
		
		
		if (   this.getValue("apply_prop")    == false 
			&& this.getValue("apply_page")    == false 
			&& this.getValue("apply_line")    == false 
			&& this.getValue("apply_table")   == false 
			&& this.getValue("apply_evi")     == false 
			&& this.getValue("apply_cond")    == false 
			&& this.getValue("apply_article") == false
		)
		{
			window.closeTab();
		}
		
		else
		{
			this.changed = false;
			window.callback(entity);
			// Unfreeze interface and reload page
			this.waitingDialog.cancel();
		}
	}
		
	this.listenEvent('formchanged', function()
	{
		document.title = document.title + " (modified)";
	});
	
	window.onbeforeunload = function(){
		if (self.changed)
			return "The data was changed. If you press OK changes will be lost";
	} 
	
	this.initialize = function()
	{
		this.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", 
		{ 
	    	width:"325px", 
			fixedcenter:true, 
			modal:true, 
			visible:false, 
			close: false
	    });
		
		this.waitingDialog.render();
		
		//batcheditForm.property = {"defaultUnit": {"id": "1"}};
		
//		batcheditForm.property = {"defaultUnit": {"id": "65"}};
//		batcheditForm.propunitsSelect.update("category=3");
		
		if (self.getValue("unitcategory") != "")
		{
//			batcheditForm.property = {"defaultUnit": {"id": batcheditForm.getValue("initial-unit")}};
//			batcheditForm.unitsSelect.update("category="+batcheditForm.getValue("initial-category"));
			
			// for the dynamic select of units it's important that the param is called "category"
			self.propunitsSelect.update("category="+self.getValue("unitcategory"), self.getValue("unit"));
			//setTimeout('batcheditForm.onPropertyChanged()', 1);
		}
	}
}

function ConditionBatchBrowser()
{
	this.controller = "properties";
	this.scope = ".conditionsscope";
	UnifiedBrowser.call(this);
	this.itemElement = "property-value";
	this.itemTemplate = "js/templates/conditionvalue-batchedit-edit.ejs";
	this.filters.scope = this.scope;
	this.url = "properties/listvaluesbatch.do";
	this.container = "ConditionsBrowser";
	this.minusCounter = -1;
	
	var self = this;
	
	this.doEditcondition = function()
	{
		var propWin = openTab("Select condition", webRoot+"properties/show.do?render-mode=popup&condition=true");
		propWin.moveBy(50, 50);
		
		propWin.callback = function(newProperty)
		{
			self.currentEntity.property = newProperty;
			if (newProperty["id"] == -1)
				alert("fetch Robert. there is a problem with conditions")
				
			self.setValue('cond-id', newProperty["id"]);
			if (self.getValue('old-cond-id') == "-1")
				self.setValue('old-cond-id', newProperty["id"]); // new condition added
			
			self.currentBlock.find('[name="condition-link"]').html(newProperty["name"]);
			self.currentBlock.find('[type="checkbox"]').attr("name","apply_cond_one_"+newProperty["id"]);
			propWin.closeTab();
			
			//setTimeout('conditionsBBrowser.currentBlock.get(0).unitsSelect.update("category='+newProperty.unitCategory["id"]+'")', 1);
			if (newProperty.qualitive == "true")
			{	
				// Update options dropbox
				setTimeout('conditionsBBrowser.currentBlock.get(0).optionsSelect.update("id='+newProperty.id+'")', 1);
			}
			else
			{	
				// Update units dropbox    newProperty.defaultUnit.name
//				setTimeout('conditionsBBrowser.currentBlock.get(0).unitsSelect.update("category='+newProperty.defaultUnit.id+'")', 1);
				setTimeout('conditionsBBrowser.currentBlock.get(0).unitsSelect.update("category='+newProperty.unitCategory["id"]+'&lightweight=1")', 1);
			}
			
			self.updateVisibility();
				
		}
	}
	
	this.onItemDrawn = function()
	{
		var dom = this.currentBlock.get(0);
		$(dom).find("[title]").tooltip();
		
		// Create dynamic select for units of condition
		var dynSelect = new UnitSelect('cond-unit', 'unit/list.do?lightweight=1', 'unit', this.currentBlock);
		dom.unitsSelect = dynSelect;
		
		
		// attach an event handler for each of unit selections
		// select the default unit for the selected condition
		dom.unitsSelect.listenEvent("update", function()
		{
			$(conditionsBBrowser.currentBlock.get(0).unitsSelect.select).val(conditionsBBrowser.currentEntity.property.defaultUnit.id);
		});
		
		// Create dynamic select for options for this condition ("qualitive conditions")
		var optSelect = new DynamicSelect('cond-option', 'properties/listoptions.do', 'option', this.currentBlock);
		dom.optionsSelect = optSelect;
		
		// Show/hide appropriate sections
		//dynSelect.fireEvent("update");
		this.updateVisibility();
	}
	
	this.updateVisibility = function()
	{
		this.currentBlock.find('span').addClass("invisible");
		this.currentBlock.find('span[name="type-' + this.currentEntity.property.type + '"]').removeClass("invisible");
	}
	
	this.doNew = function()
	{
		var entity = {"value":0.0, "property": {"id": this.minusCounter--, "name":"[...]", "unitCategory":{"unit": new Array()}}, "unit": {"name": "Click", "id": -1}};
		this.drawFromJSON(entity);
	}
	
	this.doDelete = function()
	{
		this.deleteRecord();
	}
	
	this.getRecordIdentifier = function(entity)
	{
		return entity.property["id"];
	}
}

include.plugins('view');
var batcheditForm = new BatcheditEditForm(); // batcheditform
var conditionsBBrowser = new ConditionBatchBrowser();

$(document).ready(function(){
	
	batcheditForm.initialize();
	conditionsBBrowser.initialize();
});
