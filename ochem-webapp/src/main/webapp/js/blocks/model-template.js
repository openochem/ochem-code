var form = new AjaxForm();
var self = form;
form.doSelectset = function(link)
{
	var targetLink = link;
	var baskWin = openTab("Select compound set", webRoot+"basket/show.do?render-mode=popup");
	baskWin.callback = function(basket)
	{
		form.onSetSelected(basket, targetLink);
		baskWin.closeTab();
	}
}


form.onSetSelected = function(basket, targetLink)
{
	var properties;
	var bindTo = $(targetLink).attr('bindto');
	
	if (bindTo.indexOf("validationsetid") != -1 && (form.propertyNames && form.propertyNames.length > 0))
	{
		// Check if the selected validation set contains necessary properties
		var props = array(basket.properties);
		var ok = false;
		for(var i = 0; i < props.length; i++)
		{
			if (form.propertyNames.indexOf(props[i].name) >= 0){
				ok = true;
				break;
			}
		}
		if (!ok)
		{
			window.alert("The selected validation set does not contain any of the predicted properties: " + form.propertyNames);
			return;
		}
	}
	
	self.setValue(bindTo, basket.id, basket.name);
	$("span#" + bindTo + "-more").removeClass("invisible");
	
	if (bindTo == "trainingsetid")
	{
		if (basket.mixtures > 0)
			$(".mixtures").removeClass("invisible")
		else
			$(".mixtures").addClass("invisible");
		
		form.propertyNames = new Array();
		// Enumerate the predicted properties and fetch their units
		$("#propertyBrowser").find("div:not('.invisible')").addClass("invisible");
		$(":input").find("[fromreset='true']").val('');
		
		if(basket.properties.length){
			$("#modelText").html("The model will predict these properties:");
			//$("select[name='template']").find("option[support='false']").hide();
			//$("select[name='template']").find("option[support='true']").show();
			//$("select[name='template']").find("option[support='true']").attr("selected", "selected");
			$("div[support='false']").hide();
			$("div[support='true']").show();
			selectDefaultTemplate();
			properties = basket.properties;
		} else {
			$("#modelText").html("The model will predict this property:");
			//$("select[name='template']").find("option").show();
			$("div[support]").show();
			selectDefaultTemplate();
			//$("option[support='false']").attr("selected", "selected");
			properties = array(basket.properties);
		}
		
		for(var i = 0; i < properties.length; i++){
			if(i < 41) // some magic number: changing it breaks code - some properties are not uploaded
				$("#property"+i).removeClass("invisible");
			else{
				//append property block
				$("#propertyBrowser")
					.append('<div id="property'+i+'"><a bindto="property'+i+'-id" storein="property'+i+'-title" href="properties/edit.do?id="'+properties[i].id+' tab="Property profile">'+properties[i].name +
							'</a> using unit: <select id="unit'+i+'" name="unit'+i+'"></select></div>');
				
				//append hidden values	
				$("#hidden-value")
					.append('<input type="hidden" name="property'+i+'-id" class="invisible" fromreset="true"/>' +
							'<input type="hidden" name="unitcategory'+i+'-id" class="invisible" fromreset="true"/>' +
							'<input type="hidden" name="unit'+i+'-id" class="invisible" fromreset="true"/>' +
							'<input type="hidden" name="property'+i+'-title" class="invisible" fromreset="true"/>');
			}
			
			self.propertyNames.push(properties[i].name);
			self.setValue("property"+i+"-id", properties[i].id, properties[i].name);
			self.propertyUpdated(properties[i],i);
		}
		
	}else
	{
		$(targetLink).parent("div").removeClass("invisible");
		$(targetLink).parent("div").find(".delete-link").removeClass("invisible");
	}
	
	$(document).trigger("set_selected", $(document));
	$(document).trigger("DOM_updated", $(document));
}

form.doCreate_from_another_model = function()
{
	var win = openTab("Select a template model", webRoot+"model/select.do?render-mode=popup&single=true");
	win.callback = function(model)
	{
		window.location.replace("modelconfigurator/createFromTemplate.do?model=" + model.id + "&trainingSet=" + form.getValue("trainingsetid"));
		win.closeTab();
	}
}

form.doCreate_from_template = function()
{
	var win = openTab("Select a configuration template", webRoot+"modeltemplate/show.do?render-mode=popup");
	win.callback = function(template)
	{
		window.location.replace("modelconfigurator/createFromTemplate.do?template=" + template.id + "&trainingSet=" + form.getValue("trainingsetid"));
		win.closeTab();
	}
}

form.doImporttemplate = function()
{
	window.alert("Select a file");
}

form.doMoreonset = function(link)
{
	var setLink = $("[bindto='" + $(link).attr("set") + "']");
	var setId = form.getValue($(link).attr("set"));
	openTab("Basket details - " + setLink.html(), webRoot + "basket/edit.do?render-mode=popup&id=" + setId);	
}

form.doPropertyclick = function(link)
{
	openTab("Property profile", "properties/edit.do?id=" + this.getValue($(link).attr("bindto")));
}

form.doSubmit = function()
{
	var id = self.getValue("property0-id");
	var ts = self.getValue("trainingsetid");
	var vs = self.getValue("validationsetid");
	if (ts == undefined || ts == "") 
		return window.alert("Please select a valid training set for the model");
	if (id == undefined || id == "") 
		return window.alert("Please select a valid property for the model");
		
	// Store in unit hidden field for correct back-button processing
		var i = 0;
	$("#propertyBrowser").find("div:not('.invisible')").each(function(){
		self.setValue("unit"+i+"-id", self.getValue("unit"+i));
	})	
	document.modelform.submit();
}

form.doDelete = function(link)
{
	var div = $(link).parent("div");
	var bindto = div.find("[bindto]").attr("bindto");
	//if (div.is(".first"))
	//{
		div.addClass("invisible");
		console.log(bindto);
		//div.find(".delete-link").addClass("invisible");
		div.find("[name='"+bindto+"']").removeAttr("value");
		div.find("[bindto='"+bindto+"']").html("[...]");
		$("#" + bindto + "-more").addClass("invisible");
	//}
	//else
	//{
	//	div.addClass("invisible");
	//}
}

form.propertyUpdated = function(newProperty, i)
{
	//if (newProperty.qualitive == "true")
	//	window.alert("Your property "+newProperty.name+" is qualitive. We do not currently support classification models.\nPlease select quantitive property for prediction");
	
		var selUnits = new DynamicSelect('unit'+i, 'unit/list.do', 'unit');
		selUnits.update("lightweight=1&category=" + newProperty.unitCategory.id);
		self.setValue("unitcategory"+i+"-id", newProperty.unitCategory.id);
		selUnits.listenEvent("update", function()
		{
			if (newProperty)
				if (QSPR.defaultPropertyUnit && QSPR.defaultPropertyUnit[newProperty.name]) // fix it.. how to find the unit by name
				{
					var unitId = $(selUnits.select).find('*:contains("'+QSPR.defaultPropertyUnit[newProperty.name]+'")').val();
					$(selUnits.select).val(unitId);
				}
				else
					$(selUnits.select).val(newProperty.defaultUnit.id);
		});
}

var selectDefaultTemplate = function()
{
	$(":radio").removeAttr("checked");
	if (!form.defaultTemplate)
		$(":radio:visible:first").attr("checked","checked");
	else
		$("[method="+form.defaultTemplate+"]").attr("checked","checked");
}

$(document).ready(function(){
	
	
	$("input[name=usenuminstances]").click(function(){
		if ($(this).is(':checked'))
			$("#usenuminstances").removeClass("invisible");
		else
			$("#usenuminstances").addClass("invisible");
		
	});
	
	if (getParams["test"])
	{
		$("input[name='test']").val("true");
	}

	if (getParams["upload"])
	{
		$("div[uploadable='']").remove();
		$(".not-uploadable").remove(); 
		$("input[name='upload']").val("true");
	}
	else{
		$(".uploadable").remove(); 
	}
	
	selectDefaultTemplate();
	//load property 
	var i=0;
	$("div[class='invisible']").each(function(){
		var storageElement = $("[name='"+$(this).attr('id')+"-id']");
		if (storageElement.val()){
			$(this).removeClass("invisible");
			
			var selUnits = new DynamicSelect('unit'+i, 'unit/list.do', 'unit');
			selUnits.update("lightweight=1&category=" + form.getValue("unitcategory"+i+"-id"));
			selUnits.listenEvent("update", function()
			{
				$(selUnits.select).val(form.getValue("unit"+i+"-id"));
			});
		}
		i++;
	});
	
	var onTemplateChanged;
	$("[name=template]").click(onTemplateChanged = function() {
		var consensusMethodSelected = $("[value=35]").is(":checked");
		$("#validation").setClass("invisible", consensusMethodSelected);
		$("#consensus-validation").setClass("invisible", !consensusMethodSelected);
	});
	onTemplateChanged();
	
	$("select[name='validation']").change(function()
	{
		$("div.validation").addClass("invisible");
		$("div#validation-" + $(this).val()).removeClass("invisible");	
	});
	$("select[name='validation']").change();
	
	//if (getParams['skip-configuration'])
	//{
		skipConf();
		setTimeout('skipConf()', 300); // due to a misteriously self-unchecking checkbox
//	}
});


function skipConf()
{
	var chkSkipConfig = $("input[name='skip-configuration']");
	
	chkSkipConfig.setChecked(getParams['skip-configuration']);
	chkSkipConfig.change(function(){
		$("#detailed-configuration").setClass("invisible", $(this).is(":checked"));
		$("#template-configuration").setClass("invisible", !$(this).is(":checked"));
	});
	chkSkipConfig.change();
}

