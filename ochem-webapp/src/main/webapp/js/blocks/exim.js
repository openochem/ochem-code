var lastSubmittedForm = "";

function getIframeDataJson()
{
	var data = $("#iframe").contents();
	if ((data.text() != "null") && (data.text() != ""))
		return eval("("+data.text()+")");
	else
		return null;
}

function onFrameLoad(object)
{
	$("input[type=submit]").removeAttr("disabled");
	var json = getIframeDataJson();
	if (json == null)
		return;

	var params = array(json.param);
	var logline = array(json.log["log-line"]);
	
	for (i=0; i<logline.length; i++)
	{
		var line;
		if (logline[i].type == "error")
			line = "<p class='error'>Error: "+logline[i].content+"</p>";
		else
			line = logline[i].content;
		$("#log").append(line+"<br/>");
	}
	
	if (lastSubmittedForm == "exportform")
	{
		$("#status").html("Structure available for download by a <a href='"+params[0].content+"'>direct link</a>. The download should start automatically.");
		$("#iframe").get(0).src = params[0].content;
	} else
	{
		$("#status").html("Done!");
	}
}


function onFormSubmit(object)
{
	$("input[type=submit]").attr("disabled","disabled");
	$("#log").html("");

	lastSubmittedForm = object.currentTarget.name;
	
	if (lastSubmittedForm == "exportform")
		$("#status").html("<img src='img/roller.gif'/> Please wait while the structure ie being exported. ");
	else
		$("#status").html("<img src='img/roller.gif'/> Please wait while the structure ie being imported. ");
}

function onBasketSelect(object)
{
	link = object.currentTarget;
	var basketWin = openTab("Select molecule set", webRoot + "basket/show.do?render-mode=popup");
	basketWin.callback = function(basket)
	{
		$(link).parents("form:first").find("input[type='submit']:first").removeAttr("disabled");
		$(link).parents("form:first").find("input[type='hidden']:first").attr('value', basket["id"]);
		$(link).html(basket["name"]);
		basketWin.closeTab();
	}
}