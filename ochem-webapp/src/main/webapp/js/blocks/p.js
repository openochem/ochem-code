function p()
{		
	
	AjaxForm.call(this);
	this.controller = "p";
	this.actionURL = "p/action.do";
	var self = this;
	
		this.doEditmolecule = function()
		{
			this.doSwitch();
			
			var moldata = self.getValue('moldata');
			
			if(document.JME.isActive())
				if (moldata)
					document.JME.readMolFile(URLDecode(moldata));
			
			
//			var molWin = openTab("Edit molecule", webRoot+"molecule/show.do?id="+self.getValue('n-molecule'));
//			
//			molWin.callback_extended = function(mol)
//			{
//				self.setValue("n-molecule",mol.id,mol.id);
//				self.setValue("n-data", URLDecode(mol.encodedData), mol.id);
//				$('#depiction').attr('src', 'depiction.jsp?id='+mol.id);
//				//molWin.closeTab();
//			}
		}
	
		this.doSwitch = function() {
			if ($('#dep').hasClass('invisible')) 
			{
				$('#jme').addClass('invisible');
				$('#dep').removeClass('invisible');
				jQuery('#detectBlock').removeClass("invisible");
				jQuery('#submitButton').removeClass('invisible');
				jQuery('#drawButton').removeClass('invisible');
			}
			else 
			{
				$('#dep').addClass('invisible');
				$('#jme').removeClass('invisible');
				jQuery('#detectBlock').addClass("invisible");
				jQuery('#submitButton').addClass('invisible');
				jQuery('#drawButton').addClass('invisible');
				
				jQuery('#icBlock').addClass("invisible");
				jQuery('#modelBlock').addClass("invisible");
				jQuery('#resultBlock').addClass("invisible");

			}
		}
		
		this.getActionQuery = function(action)
		{
			if (action == "detect")
			{
				jQuery('#icBlock').addClass("invisible");
				jQuery('#modelBlock').addClass("invisible");
				jQuery('#resultBlock').addClass("invisible");
				
				return "moldata="+self.getValue('moldata');
			}
			if (action == "selectIC")
			{
				var selected = $('#tabelle input:radio');
				var actionQuery = "";
				for (i = 0; i < selected.length; i++)
					if (parseInt(selected[i].value) > -1)
						actionQuery += "&selectedICs=" + selected[i].value;
				return actionQuery;
			}
			if (action == "selectModel")
			{
				var models = $("input[name=modelradios]");
				var actionQuery = "";
				for (i = 0; i < models.length; i++)
					actionQuery += "&models=" + models[i].value;
				return actionQuery;
			}	
				
			if (action == "directpredict")
			{
				return "moldata="+self.getValue('moldata');
			}	
//			query = self.parentActionQuery();
//			window.alert(query);
//			return query;
		}
	
		this.onDetectSuccess = function(xml)
		{
			var icmols = array(xml.others.iCmol);
			var html="<table class=\"inner\">";
			for (i = 0; i < icmols.length; i++)
			{
				// at the moment there is only one molecule, update other images later
//				$('#depiction').attr('src', 'depictionPka.jsp?id='+ (parseInt(icmols[i].molId)));
				
				
				var groups = array(icmols[i].groups);
				for (j = 0; j < groups.length; j++)
				{	
					html += "<tr><td><h2>" + groups[j].group + "</h2>";
					
					var centers = array(groups[j].centers);
					for (k = 0; k < centers.length; k++)
					{	
						var centerIndex = parseInt(centers[k].index);
						html += "<input type=\"radio\" name=\"icradios\" send=\"1\" value=\"" + centerIndex + "\"> Center Nr.: " + (centerIndex + 1) + "<br/>";
					}
						
				
					html += "</td><td>";
					html += "<div>Suggested Models</div>";
					
					var models = array(groups[j].suggestedModel);
					for (k = 0; k < models.length; k++)
					{
						var modelname = models[k].name;
						html += "<input type=\"radio\" name=\"modelradios\" send=\"1\" value=\"" + models[k].name + "\">" + models[k].name + "<br/>";
					}
					
					html += "</td></tr>";
				}
			}			
			html+="</table>";

			$('#tabelle').html(html);
			
//			$('#selectICButton').removeClass("invisible");
			$('#icBlock').removeClass("invisible");
			
			$(".icselect").change(function () {
			    var str = "";
			    $(".icselect").each(function () {
			          str += $(this).text() + " ";
			          
			          if ( $("option:selected").length == 0 )
			          {
			        	  $(this).val($('option:first', this).val());
			          }
			          
			        });
			    $("#output").text(str);
			}).change();
			
			$('#selectModelButton').removeClass("invisible");
//			$('#modelBlock').removeClass("invisible");
			
		}
		
		this.onSelecticSuccess = function(xml)
		{
//			var html = 	"<table class=\"inner\">" +
//							"<tr><td><select id=\"mSelect\">" +
//								"<option value=\"0\">Tehan Phenols</option>" +	
//								"<option value=\"tehan11\">Tehan Phenols</option>" +
//							"</select></td></tr>" +
//						"<table>";
//			
//			$('#modelle').html(html);
			
			var models = array(xml.others.model);
			var html="<table class=\"inner\">";
			for (i = 0; i < models.length; i++)
			{
				html += "<tr><td><select name=\"models" + i + "\" send=\"1\">"; // +
							//"<option value=\"0\">--model--</option>";
				html += "<option value=\"" + models[i].id + "\">" + models[i].name + "</option>";
				html += "</select></td></tr>";
				
			}			
			html+="</table>";
			$('#modelle').html(html);
			
			
			$('#selectModelButton').removeClass("invisible");
			$('#modelBlock').removeClass("invisible");
		}
		
		this.onSelectmodelSuccess = function(xml)
		{
			var stop;
			alert("Model is submitted");
			$('#result').removeClass("invisible");
			$('#resultBlock').removeClass("invisible");
			
			start();
		}
		
		this.onDirectpredictSuccess = function(xml)
		{
//			var color = "green";
			var ics = new Array();
			var icmols = array(xml.others.iCmol);
			var html = "";
				
			for (i = 0; i < icmols.length; i++) // length is 1
			{	// two tables per icmol change if necessary

				html += "<h3>detected groups and ionisation sites</h3><table id=\"inner\">";
				html += "<tr><th>group</th><th>index of IC</th><th>prediction</th></tr>";
				
				var groups = array(icmols[i].groups);
				for (j = 0; j < groups.length; j++)
				{	
					var centers = array(groups[j].centers);
					for (k = 0; k < centers.length; k++)
					{
						html += "<tr>";
						var centerIndex = parseInt(centers[k].index) + 1;
						ics.push(centerIndex); // += centerIndex + ",";
						if (k == 0) { 
							html += ("<td>" + groups[j].group + "</td><td>" + centerIndex + "</td><td class=\"value\"><img src=\"img/roller_transparent.gif\"/></td>")
						} else {
							html += ("<td></td><td>" + centerIndex + "</td><td class=\"value\"><img src=\"img/roller_transparent.gif\"/></td>");
						}
						html += "</tr>";
					}
				}
				
				html += "</table><br/><br/>";
				html += "<h3>Used models</h3><table id=\"inner\">";
				
				for (j = 0; j < groups.length; j++)
				{	
					var models = array(groups[j].suggestedModels);
					if (models.length == 0)
					{
						html += "<tr>";
						html += ("<td>" + groups[j].group + "</td><td>No supporting model, yet.</td>") 
						html += "</tr>";
					} else
					{
						for (k = 0; k < models.length; k++)
						{
							html += "<tr>";
							var modelname = models[k].name;
	//						html += "<input type=\"radio\" name=\"modelradios\" send=\"1\" value=\"" + models[k].name + "\">" + models[k].name + "<br/>";
							if (k == 0) { 
								html += ("<td>" + groups[j].group + "</td><td>" + modelname + "</td>") 
							} else {
								html += ("<td></td><td>" + modelname + "</td>");
							}
							html += "</tr>";
						}
					}
				}
				html+="</table><br/>";
			}			
			
			$('#results').html(html);
			
			
			
			var moldata = self.getValue('moldata');
			var decorde = URLDecode(moldata);
			//decorde
			decorde += ">  <ics>\n" + ics.join() + "\n";
			var encMol = URLEncode(decorde);
			jQuery('#depiction').attr('src', "depictionPka.jsp?mol=" + encMol + "&color=ffffff&w=400&h=300");
			
			
			$('#results').removeClass("invisible");
			$('#resultBlock').removeClass("invisible");
			
			start();
		}
		
//		$(document).ready(function()
//		{
//			$("select[name|='ics']").change(function () {
//				alert("hey robbo");
//			}).trigger('change');
//		});
		
		
		
		this.doSelectjme = function() {
			var jmeMol2 = document.JME.molFile();
			jmeMol2 = jmeMol2.replace(/JME[^\n]+/, "123");
			var encMol = URLEncode(jmeMol2);
			
			jQuery('#depiction').attr('src', "depictionPka.jsp?mol=" + encMol + "&color=ffffff&w=400&h=300");
			
			self.setValue("moldata", encMol);
			
			this.doSwitch();
			
		}
		
		this.doCanceljme = function()
		{
			this.doSwitch();
		}
		
		this.initialize = function()
		{
			var moldata = self.getValue('moldata');
			$('#depiction').attr('src', "depictionPka.jsp?mol=" + moldata + "&color=ffffff&w=400&h=300");
		}
		
}		

include.plugins('view');
var pForm = new p(); // batcheditform

$(document).ready(function(){
	pForm.initialize();
});




var ajaxpp = new QSPR.Ajax();
var timer;

checkStatus = function()
{
	ajaxpp.url = 'p/status.do';
	ajaxpp.send({
		data:'',
		success: function(xml)
		{
			
			var status = xml.message;
			
			if ( ! xml.message) //status == "Finished")
			{
				clearInterval(timer);	
				jQuery("#status").addClass('invisible');
				
				var tdcount = 0;
				var tdVals = $("#results .value");
				
				var groups = array(xml.iCmol.groups);
				for (i = 0; i < groups.length; i++)
				{
					var ics = array(groups[i].centers);
					for (j = 0; j < ics.length; j++)
					{
						tdVals.eq(tdcount++).html(ics[j].prediction);
					}
				}
				
				
//				var dt = xml.datatable;
//				var rows = array(dt.rows.row);
//				var res = array(rows[0].string);
//				var tdVals = $("#results .value");
//				for (i = 0; i < res.length; i++) 
//				{
//					var td = res[i];
//					tdVals.eq(i).html(td);
//				}
				
				$("#result").html(xml.innerUrl);
				
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
	});
}

function start()
{
	$("#status").html('Starting...');
	ajaxpp.url = 'p/start.do';
	ajaxpp.send({
		success: function()
		{
			$("#status").html('Requesting status...');

			timer = setInterval("checkStatus()", 3000);
		}
	});
}