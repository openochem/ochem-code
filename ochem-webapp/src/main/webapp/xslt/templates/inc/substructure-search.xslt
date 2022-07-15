<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template name="substructure-search">
	<script type="text/javascript" src="js/commons/mol-editors.js"></script>
	<style type="text/css">
		.sketcher-frame {
			position: absolute;
		}
	
		.hidden-frame
		{
			position: absolute;
			left: -1000px;
			top: -1000px;
		}
		
		.green-button, .green-button:visited, .green-button:active {
			-moz-border-radius: 15px;
			border-radius: 15px;
			background-color: #595;
			color: white;	
			padding: 3px 6px;
			font-size: 90%; text-transform: uppercase;
			margin-top: 20px;
		}
		.green-button:hover {text-decoration: none; background-color: #050;}
		
		#sketcherDialog TEXTAREA {width: 200px;}
	</style>
		<h2>Similarity/substructure search</h2>
			Draw a structure and search all the molecules containing it or similar to it<br/><br/>
			<a href="javascript:void()" onclick="sketcherDialog.drawStructure(); return false;" style="float: left; margin-right: 15px; border: 1px solid black;">
				<img src="img/click-to-draw.png" id="ss-depiction" title="Click to draw a molecular structure in a visual editor"/>
			</a>
			<input name="xemistry-smiles" type="hidden" filter="1"/>
			<div class="ss-details invisible" style="float: left; width: 130px; ">
				Search type:<br/>
				<input type="radio" name="xemistry-sstype" value="substructure" checked="1" filter="1"/>Substructure<br/>
				<input type="radio" name="xemistry-sstype" value="similarity" filter="1"/>Similarity<br/>
				<br/>
				<div id="xemistry-score">
				Similarity score at least <input type="text" name="xemistry-similarity-cutoff" filter="1" style="width: 50px; text-align: center;" value="85"/>%<br/>
				</div>
				
				<br/>
				<a href="#" onclick="sketcherDialog.clearFilter(); return false;" id="remove-ss-filter" class="invisible fb-button">remove the filter</a>
				
			</div>
			
			<div id="sketcherDialog"> 
			    <div class="hd">Draw a (sub)structure</div> 
			    <div class="bd">
			    <a href="#" onclick="sketcherDialog.submitMol(); return false;" class="fancy-button">submit the molecule below</a><br/><br/>
			    <table height="100%">
			    	<tr>
			    		<td width="550" valign="top">
			       			 <iframe src="../editor.html" id="sketch" class="sketcher-frame hidden-frame"></iframe>
			    		</td>
			    		<td valign="top"> Structure from text<br/><nobr>(SMILES, SD-file or chemical name):</nobr><br/>
			    		<textarea name="molecule-str"></textarea>
			    		<br/>
			    		<br/>
			    		<a href="#"  onclick="sketcherDialog.loadFromText(); return false;" class="fb-button">load from text</a>
			    		<span>
			    			<span id="progress" class="invisible"></span>
			    			<span class="lo-progress-bar invisible"><img src="img/roller_small.gif"/> Loading...</span>
			    		</span>
			    		</td>
			    	</tr>
			    </table>
			   
			    
			    </div> 
			</div>
			
			<script language="javascript">
			
			$("[name='xemistry-sstype']").change(function(){
				$("#xemistry-score").setClass("invisible", $(this).val() != "similarity");
			});
			$("[name='xemistry-sstype']").change();
			
			$("[name='xemistry-sstype']").change();
			
			
			
			var sketcherDialog = new YAHOO.widget.Dialog("sketcherDialog", { 
		    	width:"750px", 
		    	height: "600px",
				fixedcenter:true, 
				modal:true, 
				visible:false,
				postmethod: "none",
				buttons:[]
		    });
		    
		    $("[name='molecule-str']").keydown(function(e) {
				if(e.keyCode == 13){
				     sketcherDialog.loadFromText();  
			    }
			});
		    
		    sketcherDialog.editor = new JSME();
		    		    
		    var oldHide = sketcherDialog.hide;
			sketcherDialog.hide = function()
			{
				oldHide.call(sketcherDialog);
				$("#sketch").addClass("hidden-frame");
			}
		    
		    sketcherDialog.submitMol = function()
		    {
		    	var mol = "1" + sketcherDialog.editor.getMolecule();
		    	if (mol.length &lt;= 52)
		    	{
		    		window.alert("You did not draw any molecule");
		    	}
		    	else
		    	{
		    		sketcherDialog.setMolecule(mol, true);
		    	}
		    }
		    
		    sketcherDialog.setMolecule = function(mol, reload) {
		    	var tmp = sketcherDialog.filters.onFilterChange;
	    		sketcherDialog.filters.onFilterChange = function(){};
		    	sketcherDialog.filters.setValue("xemistry-smiles", mol);
		    	sketcherDialog.hide();
		    	$("#sketch").addClass("hidden-frame");
		    	if (mol.match(/M[0-9]+/))
		    		$("#ss-depiction").attr("src", "depiction.jsp?mp2=" + mol.substring(1));
		    	else if (mol.match(/DP[0-9]+/))
		    		$("#ss-depiction").attr("src", "depiction.jsp?id=" + mol.substring(2));
		    	else
		    		$("#ss-depiction").attr("src", "depiction.jsp?mol=" + URLEncode(mol));
		    	sketcherDialog.updateVisibility();
		    	sketcherDialog.filters.onFilterChange = tmp;
		    	if (reload)
		    		sketcherDialog.filters.onFilterChange();
		    }
		    
		    sketcherDialog.loadFromText = function() {
		    	var ajax = new QSPR.Ajax();
		    	ajax.send({
		    		url: "molecule/getFromText.do?out=json",
		    		data: "molecule-str=" + URLEncode($("[name='molecule-str']").val()),
		    		tracker: new SimpleTracker($("#progress")),
		    		success: function(resp) {
		    			var moleculeStr = URLDecode(resp.molecule.encodedData);
		    			console.log(moleculeStr);
		    			sketcherDialog.editor.setMolecule(moleculeStr);
		    		} 	
		    	});
		    }
		    
		    sketcherDialog.clearFilter = function()
		    {
		    	var tmp = sketcherDialog.filters.onFilterChange;
		    	sketcherDialog.filters.onFilterChange = function(){};
		    	
		    	sketcherDialog.filters.setValue("xemistry-smiles", "");
		    	$("#ss-depiction").attr("src", "img/click-to-draw.png");
		    	
		    	sketcherDialog.updateVisibility();
			    sketcherDialog.filters.onFilterChange = tmp;
			    sketcherDialog.filters.onFilterChange();
		    	
		    }
		    
		    sketcherDialog.drawStructure = function()
			{
				$("#sketch").removeClass("hidden-frame");
				sketcherDialog.show();
			}
			
			sketcherDialog.updateVisibility = function()
			{
				$(".ss-details").setClass("invisible", !sketcherDialog.filters.isSet("xemistry-smiles"));
				$("#remove-ss-filter").setClass("invisible",  !sketcherDialog.filters.isSet("xemistry-smiles"));
				$("[name='xemistry-sstype']").change();
			}
			
			function getSketcher() {
				if (typeof document.getElementById("sketch") != 'undefined') {
					var eframeWin = document.getElementById("sketch").contentWindow;
					if(typeof eframeWin != 'undefined') {
						return eframeWin.marvin.sketch;
					}
				}
				return null;
			}
		    
		    sketcherDialog.render();
		    	
			</script>
	</xsl:template>
</xsl:stylesheet>