var emptymolecule = "%0A%0A%0A%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200999%20V2000%0AM%20%20END%0A";

function p() {
	
	AjaxForm.call(this);
	this.actionURL = "physprop/action.do"; // define the controller method where the action will be handled
	var self = this;
	this.molEditor = new JSME("#mol-editor");
	this.molEditor.promise().then(function(){
		self.doSwitch();
		var moldata = self.getValue('moldata');
		
		if (moldata)
			self.molEditor.setMolecule(URLDecode(moldata));
	});
	
	this.doEditmolecule = function() {
		this.doSwitch();
		var moldata = self.getValue('moldata');
		
		if (moldata)
			this.molEditor.setMolecule(URLDecode(moldata));
	}
	
	this.getActionQuery = function(action)
	{
		if (action == "submit")
		{
			return "moldata=" + self.getValue('moldata') + "&modelIds=" + self.selectedIds.join(",");
		}
		if (action == "stop")
			return "";
	}
	
	this.doStore = function()
	{
		var div = $("<div/>");
		div.append("<img src='depiction.jsp?w=100&h=100&mol=" + self.predictedMol + "' align='left'/>");
		
		var table = $("<table/>");
		$("input[model-id]:checked").each(function(index){
			table.append("<tr><td>" + $(this).closest("tr").find(".value").html() + "</td><td></td></tr>");
		});
		div.append(table);
		div.append("<a class='delete' action='delete' title='Remove this prediction'>x</a>");
		
		$("#comparison").append(div);
		self.attachActions(div);
		
		if ($("#comparison > DIV").length >= 2)
			$("#comparison").removeClass("invisible");
	}
	
	this.doDelete = function(link) {
		$(link).closest('div').remove();
	}
	
	this.beforeSubmit = function() 
	{
		self.selectedIds = $("[model-id]:checked").map(function() {return $(this).attr("model-id")}).get();
		if (self.selectedIds.length == 0)
		{
			window.alert("Please, select at least one model to predict");
			return false;
		}
		
		if (self.getValue('moldata') == emptymolecule)
		{
			alert("Please select a valid molecule");
			return false;
		}
		
		return true;
	}
	
	this.onSubmitSuccess = function() 
	{
		self.predictedMol = self.getValue("moldata");
		$("#results input[model-id]").closest("tr").find(".value").html('-');
		$("#results input[model-id]:checked").closest("tr").find(".value").html('<img src="img/roller_transparent.gif"/>');
		//$("#results .value").html('<img src="img/roller_transparent.gif"/>');
		$(".running").removeClass("invisible");
		start();
	}
	
	this.doSwitch = function() {
		if ($('#dep').hasClass('invisible')) 
		{
			$('#mol-editor').addClass('invisible');
			$('#dep').removeClass('invisible');
			jQuery('#submitButton').removeClass('invisible');
			jQuery('#drawButton').removeClass('invisible');
		}
		else 
		{
			$('#dep').addClass('invisible');
			$('#mol-editor').removeClass('invisible');
			jQuery('#submitButton').addClass('invisible');
			jQuery('#drawButton').addClass('invisible');

		}
	}
	
	this.doSelect = function() {
		var molStructure = this.molEditor.getMolecule();
		molStructure = molStructure.replace(/JME[^\n]+/, "123");
		var encMol = URLEncode(molStructure);
		
		if (encMol == emptymolecule)
		{
			alert("Please select a valid molecule");
			return;
		}	
		
		jQuery('#depiction').attr('src', "depiction.jsp?mol=" + encMol + "&color=ffffff&w=400&h=300");
		
		self.setValue("moldata", encMol);
		
		this.doSwitch();
	}
	
	this.doCancel = function() {
		this.doSwitch();
	}
	
	this.initialize = function()
	{
		var moldata = self.getValue('moldata');
		$('#depiction').attr('src', "depiction.jsp?mol=" + moldata + "&color=ffffff&w=400&h=300");
		
		this.ajax.waitingDialog = new YAHOO.widget.Dialog("waitingDialog", { 
	    	width:"325px", 
			fixedcenter:true, 
			modal:true, 
			visible:false, 
			close: false,
			postmethod: "none"
	    });
		this.ajax.waitingDialog.render();
	}
}

var physprop = new p(); 

$(document).ready(function(){
	physprop.initialize();
});

var ajaxpp = new QSPR.Ajax();
var timer;
var lock = false;

checkStatus = function()
{
	if (lock)
		return;
	lock = true;
	ajaxpp.url = 'physprop/status.do?nodb=1';
	ajaxpp.send({
		data:'',
		success: function(xml)
		{
			
			var status = xml.message;
			
			if ( ! xml.message) //status == "Finished")
			{
				clearInterval(timer);	
				jQuery("#status").addClass('invisible');
				
				
				var dt = xml.datatable;
				var rows = array(dt.rows.row);
				var res = array(rows[0].string);
				console.log(res);
				var tdVals = $("#results .value");
				tdVals.html("-");

				// Restore the checkboxes to the state when the calculations were started
				$("input[model-id]").removeAttr("checked");
				physprop.selectedIds.map(function(id){
					$("input[model-id=" + id + "]").attr("checked", 1);
				});
				
				// Display the prediction values 
				$("input[model-id]:checked").each(function(index){
					$(this).closest("tr").find(".value").html(res[index]);
				});
				
				jQuery(".running").addClass("invisible");
				jQuery('#stopButton').addClass('invisible');
				jQuery(".foot").removeClass('invisible');
				$(".store").removeClass("invisible");
				physprop.doStore();
				
			}
			else 
			{
				if (status.message.substring(0, 5) == "Error")
				{
					$("#progress-img").attr('src', 'img/icons/error.jpg');		
					clearInterval(timer);	
				}
				$("#status").html(status.message.replace(/\_\$\$\_/g, "<br/>"));
			}
		},
		error: function(msg)
		{
			$("#status").html('QSPR server is not available');
		}
	}).done(function() {
		lock = false;
	});
}

function start()
{
	if (timer)
		clearInterval(timer);
	$("#status").removeClass('invisible');
	$("#status").html('Starting...');
	ajaxpp.url = 'physprop/start.do';
	ajaxpp.send({
		success: function()
		{
			$(".store").addClass("invisible");
			$("#status").html('Requesting status...');

			timer = setInterval("checkStatus()", 1000);
		}
	});
}


