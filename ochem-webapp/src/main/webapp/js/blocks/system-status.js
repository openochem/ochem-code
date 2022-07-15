var sysstat = new SystemStatus();

function SystemStatus()
{
	var tests = new Array();
	var dates = new Array();
	var subTests = new Array();
	var table = $("table#results");
	
	this.addTestResult = function(id, testName, succeeded, time, duration, methodName)
	{
		if (methodName != '')
		{
			testName = methodName;
			testName = testName.replace('Test', '');
			testName = testName.replace('test', '');
		}
		
		var dtParts = time.split(" ");
		var date = dtParts[0];
		var hour = 1.0 * dtParts[1].split(":")[0];
		
		// Add a cell for a date
		if ($.inArray(date, dates) == -1)
		{
			var td;
			dates.push(date);
			$("tr.header").append($('<td></td>').attr('date', date).html(date));
			$("tr[res]").each(function(){
				$(this).append(td = $('<td></td>').attr('date', date));
				for (var i = 0; i < 24; i++)
					td.append($("<span></span>").attr("class", "gray").attr("hour", i));
			});
		}
		
		
		// the table row is identified by its class therefore no spaces, brackets, special chars are allowed
		var className = testName;
		className = className.replace(/\s/g, "_");
		className = className.replace(/\[/g, "_");
		className = className.replace(/\]/g, "");
		className = className.replace(/\(/g, "_");
		className = className.replace(/\)/g, "_");
		className = className.replace(/\+/g, "_");
		className = className.replace(/\./g, "_");
		className = className.replace(/\__/g, "_");
		className = className.replace(/\__/g, "_");
		
		// Add a row for a test
		if ($.inArray(testName, tests) == -1)
		{
			tests.push(testName);
			var td = $('<td class="first"></td>');
			td.html(testName);
			var tr = $('<tr res="1"></tr>');
			tr.attr("class", className);
			tr.append(td);
			for (var i in dates)
			{
				tr.append(td = $('<td></td>').attr('date', dates[i]));
				for (var i = 0; i < 24; i++)
					td.append($("<span></span>").attr("class", "gray").attr("hour", i));
			}
			table.append(tr);
		}
		
		var targetCell = table.find("tr." + className + " td[date=" + date + "]");
		var div = targetCell.find("span[hour=" + hour + "]");
		var tooltip = time;
		var detailedStatus = $("#details" + id).html();
		if (detailedStatus != null)
			tooltip += "\n" + detailedStatus;
		tooltip += " (" + (duration) + " sec.)"
		div.attr("title",  tooltip);
		div.attr("class", succeeded ? "green" : "red");
		//targetCell.append(div);
	}	
	
	this.highlightRecentlyFailedTests = function()
	{
		$("table#results").find("tr").each(function(){
			var results = $(this).find("span.green, span.red");
			$("#failures-summary").setClass("invisible", results.length == 0);
			if (results.length > 0)
			{
				$(this).find("td").eq(0).setClass("recently-failed", results.eq(results.length - 1).is(".red"));
				var failuresDiv = $("#failures-summary DIV");
				failuresDiv
			}
		});
	}
	
	this.getSeleniumTestArray = function(className, methodName)
	{
		if (className.length == 0 || methodName.length == 0)
			return;
		
		// typical examples are:
		// "qspr.tests.TestBatchUpload30" (class)
		// "qspr.tests.webDriver.TestControllerPhysprop" (class)
		// "testPagerItemAmount" (method)
		className = className.replace('qspr.tests.', '').replace('webDriver.', '').replace('Test', '');
		methodName = methodName.replace('test', '').replace('Test', '');
						
		// Add a row for a test class
		var len = tests.length;
		var index = $.inArray(className, tests);
		if (index == -1)
		{
			subTests[len] = new Array();
			subTests[len][0] = methodName;
			tests.push(className);
			
			var td = $('<td class="first"></td>');
			td.html(className);
			var tr = $('<tr res="1"></tr>');
			tr.attr("class", className);
			tr.append(td);
			table.append(tr);
		}
		else
		{
			if (subTests[index].indexOf(methodName) == -1)
			{
				var subLen = subTests[index].length;
				subTests[index][subLen] = methodName;
			}
		}
	}
	
	this.addSeleniumTest = function(id, className, methodName, succeeded, time, duration)
	{
		if (className.length == 0 || methodName.length == 0)
			return;
		
		className = className.replace('qspr.tests.', '').replace('webDriver.', '').replace('Test', '');
		methodName = methodName.replace('test', '').replace('Test', '');
		
		var dtParts = time.split(" ");
		var date = dtParts[0];
			
		var hour = 1.0 * dtParts[1].split(":")[0];
		var index = $.inArray(className, tests);
		var subTestIndex = subTests[index].indexOf(methodName);		
		
		// Add a cell for a date
		if ($.inArray(date, dates) == -1)
		{
			var td;
			dates.push(date);
			$("tr.header").append($('<td></td>').attr('date', date).html(date));
			
			// Fill the column with grey cells
			for(i in tests)
			{
				var targetTR = table.find("tr." + tests[i]);
				targetTR.append(td = $('<td></td>').attr('date', date));
				for (var j = 0; j < subTests[i].length; j++)
				{
					td.append($("<span></span>").attr("class", "gray").attr("subtest", j));
				}
				td.prepend($("<span></span>").attr("class", "dura").attr("title", "Total duration of all tests for this set of tests").text("0 sec"));
			}
		}
		
		// dye the grey cells with red or green			
		var targetCell = table.find("tr." + className + " td[date=" + date + "]");
		var span = targetCell.find("span[subtest=" + subTestIndex + "]");
		var tooltip = time + " [" + methodName + "] ";
		if ( ! succeeded)
			tooltip += $("#details" + id).html();
		tooltip += " (" + duration + " sec.)"
		span.attr("title",  tooltip);
		span.attr("class", succeeded ? "green" : "red");
		
		var dura = targetCell.find("span.dura");
		var d = parseInt(dura.text().split(" ")[0]);
		var dd = d + parseInt(duration);
		
		dura.text("" + dd + " sec");
		
		//targetCell.append(div);
	}

}

