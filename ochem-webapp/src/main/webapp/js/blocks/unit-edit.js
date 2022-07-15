var ajaxForm = new EditForm();

ajaxForm.scope = ".EditForm";
ajaxForm.actionURL = "unit/action.do";
ajaxForm.itemElement = 'unit';
ajaxForm.restriction = "form";
ajaxForm.beforeEdit = function()
{	
	var propname = $("input[name='name']").attr("value");
	var desc = $("textarea[name='description']").attr("value");
	desc = desc.trim();

	if(propname.length < 2 || propname.length > 40)
	{
		window.alert("Unit name is too short or too long and can not be added");
		return false;
	}
	if(desc.length < 50)
	{
		window.alert("Unit description is too short and can not be added");
		return false;
	}
	return true;
}

$(document).ready(function(){
	$("select[name=category]").change(function(){
		$("#def-unit").html("Default unit for this category is <b>"+$(this).find("option:selected").attr("defaultUnit")+"</b>")
	});
	
	$("select[name=category]").change();
});

ajaxForm.onVerifySuccess = function(reply)
{ 
	var name = this.getValue("name");
	window.alert(reply["double-list"].value[0] + " "+ name+" = " + reply["double-list"].value[1] + " of default unit = " + reply["double-list"].value[2] + " " +name);
}
	
function checkUnitName(obj)
{
	if(obj.value.length <= 1)
	{
		$("#unitName").css("color","red").html("error: unit name length is "+obj.value.length+" and can not be saved (min. 2 characters)");
	}else
	if (obj.value.length > 1 && obj.value.length < 5)
	{
		$("#unitName").css("color","#FA58F4").html("warning: unit name length is "+obj.value.length+" and may not be informative");
	}
	else if(obj.value.length >= 5 && obj.value.length < 40)
	{
		$("#unitName").css("color","green").html("unit name is good and length is "+obj.value.length);
	}
	else if(obj.value.length >= 40)
	{
		$("#unitName").css("color","red").html("error: unit name length is "+obj.value.length+", its should be less than 40 characters and can not be saved (max. 40 characters)");
	}
}

function checkUnitDesc(obj)
{
	if(obj.value.length <= 50)
	{
		$("#unitDesc").css("color","red").html("error: unit description length is "+obj.value.length+" and can not be saved (min. 50 characters)");
	}else
	if(obj.value.length > 50 && obj.value.length < 200)
	{
		$("#unitDesc").css("color","#FA58F4").html("warning: unit description length is "+obj.value.length+" and may not be informative");
	}
	else if(obj.value.length > 200)
	{
		$("#unitDesc").css("color","green").html("unit description is good");
	}
}