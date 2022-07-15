<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../helper.xslt" />

	<xsl:template name="content">
		<style type="text/css">
			.actions A
			{
				margin-right: 10px;
			}
			
			#show-filters
			{
				position: fixed;
				left: 0px;
				top: 50%;
				background-color: #44F;
				color: white;
			}
			
		</style>
		<title>Property Predictor</title>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script>
			include.plugins('view');
			var compoundBrowser = new CompoundBrowser();
			var conditionsFilter = new ConditionsFilter();
			//var basketActionable = new BasketActionable();
			$(document).ready(function() {
			
			function hidefilters()
			{
				$("#applymodeldata").addClass("invisible");
				$(".show-button").removeClass("invisible");
			}
			
			this.doShowfilters = function()
			{
				$("#applymodeldata").removeClass("invisible");
				$(".show-button").addClass("invisible");
			}
			
			
				if (getParams["trash"])
					$("#deleteselected").attr("title", "Delete selected records permanently");
				compoundBrowser.options.highlight_row = "no";
				compoundBrowser.initialize();
				compoundBrowser.recordMenu = new YAHOO.widget.Menu("basicmenu");
				compoundBrowser.recordMenu.render();
				conditionsFilter.initialize(true); 
				//$('a[title]').tooltip({showURL: false, delay: 0});				
				//basketActionable.initialize();				
			});
			
			
			
			
			
		</script>

	<table width="100%" height="100%" cellspacing="0">
		<tr>
			<td class="itunes-up" colspan="2">
				<img src="img/icons/compounds.png" id="page-icon" />
				<h1>
					iPRIOR
				</h1>
			</td>
		</tr>
		<tr>
			<td class="browserScope itunes-left" id="applymodeldata">
				<div id="applymodeldataDiv">
					<a action="hidefilters" class="close-button" title="Hide filters panel">x</a>
					<h1>Molecular Data</h1>
					<div class="openable opened">
						<h1>Source</h1>
						<div class="openable-content" style="min-height: 74px;">
							Article/Source
							<a href="javascript:void(0)" action="editarticle" name="article-link">[select]</a>
							<br />
							<div id="ac-source" style="height: 25px;"></div>
							<table>
								<tr>
									<td>Page</td>
									<td>Table</td>
								</tr>
								<tr>
									<td>
										<input class="small right" name="page" type="text"
											prompt="All" value="" filter="1" />
									</td>
									<td>
										<input class="small right" name="table" type="text"
											prompt="All" value="" filter="1" />
									</td>
								</tr>
							</table>
						</div>
					</div>
					<div class="openable opened">
						<h1>Property</h1>
						<div class="openable-content">
							Activity/Property
							<a action="editproperty" name="property-link">[select]</a>
							<div id="ac-property" style="height: 25px;"></div>

						</div>
					</div>
					<div class="openable">
						<h1>Conditions</h1>
						<div class="openable-content conditionsscope">
							<a action="newfilter" restrict="conditionsfilter">[add]</a>
							<div id="ConditionsFilter"></div>
						</div>
					</div>

					<div class="openable opened">
						<h1>Molecule</h1>
						<div class="openable-content">
							Name / OCHEM ID
							<a help="OCHEM-id-hint">[?]</a>
							/ Inchi-Key
							<input type="text" name="name" prompt="" value="" filter="1" />
							<span class="invisible" id="OCHEM-id-hint">
								Examples of OCHEM IDs:
								<br />
								<br />
								M1234 for a specific molecule
								<br />
								M1234* for a molecule and all its isomers
								<br />
								R1234 for a particular record
								<br />
								R1234* for a record and all the similar records
								<br />
							</span>
							<p id="moleculeSearch">
								<a href="javascript:void()" onclick="compoundBrowser.fragSearch(-1); return false;">[search by fragment]</a>
							</p>
							<!-- <p id="cadasterSearch"><a href="javascript:void()" onclick="compoundBrowser.cadasterSearch(-1); 
								return false;">[cadaster substructure search]</a></p> -->
						</div>
						Molecular mass
						<br />
						between
						<input type="text" class="w50" name="mmlo" prompt="" value=""
							filter="1" />
						and
						<input type="text" class="w50" name="mmhi" prompt="" value=""
							filter="1" />
					</div>


					<div class="openable opened">
						<h1>Miscellaneous</h1>
						<div id="modification-time" class="invisible">
							Modification record id, minutes
							<input type="text" name="interval" prompt="" value=""
								filter="1" />
							<br />
						</div>
						<div class="openable-content">
							Current set
							<select id="basket-menu" name="basket-select" filter="1">
								<option value="-1">Loading...</option>
							</select>
							<input type="hidden" name="compare-baskets" value=""
								filter="1" />
							<input type="hidden" name="join-baskets" value="" filter="1" />
							<input type="hidden" name="unit" filter="1" />
							<br />

							<xsl:if test="/model/session/user">
								Records by introducers:
								<br />
								<select name="introducer" filter="1">
									<option value="">All users</option>
									<optgroup label="My records">
										<option value="my-records">Records I introduced or modified</option>
										<option value="{/model/session/user/@id}">Records I introduced</option>
										<option value="modified-by-me">Records I modified</option>
										<option value="modified-by-others">My records modified by others</option>
									</optgroup>
									<xsl:if test="/model/session/user/group">
										<optgroup label="{/model/session/user/group/name}">
											<option value="group">
												<xsl:value-of select="/model/session/user/group/name" />
												records
											</option>
											<xsl:for-each select="/model/session/user/group/member">
												<option value="{@id}">
													<xsl:value-of select="@login" />
												</option>
											</xsl:for-each>
										</optgroup>
									</xsl:if>
								</select>
								<br />
							</xsl:if>

							<input type="checkbox" name="experimental" value="true"
								filter="1" class="checkbox" />
							Original records
							<br />
							<input type="checkbox" name="primaryrecords" value="true"
								filter="1" class="checkbox" />
							Primary records
							<br />
							<input type="checkbox" name="validated" value="true"
								filter="1" class="checkbox" />
							Not validated
							<br />
							<input type="checkbox" name="error" value="true" filter="1"
								class="checkbox" />
							Error records
							<br />
							<input type="checkbox" name="inchierror" value="true"
								filter="1" class="checkbox" />
							Error inchies
							<br />
							<input type="checkbox" name="nameerror" value="true"
								filter="1" class="checkbox" />
							Mismatching names
							<br />
							<input type="checkbox" name="yesstereo" value="true"
								filter="1" class="checkbox padded" disabled="true" />
							Include stereochem.
							<br />
							<input type="checkbox" name="emptymol" value="true" filter="1"
								class="checkbox" />
							Empty molecules

							<div id="duplicate-filter" class="invisible">
								<br />
								<input type="checkbox" name="duplicates" value="true"
									filter="1" class="checkbox" />
								Duplicates
								<br />
								<input type="checkbox" name="nostereo" value="true"
									filter="1" class="checkbox padded" disabled="true" />
								No stereochemistry
							</div>
							<br />
							<br />
							Sort by:
							<br />
							<nobr>
								<select name="sortby" filter="1">
									<option value="">Creation time</option>
									<option value="time">Modification time</option>
									<option value="value">Numerical value</option>
									<option value="canonicalValue">Converted numerical value</option>
									<option value="molprop">Compound</option>
									<option value="mol.molWeight">Molecular weight</option>
									<option value="property">Property</option>
									<option value="artMolId">Identifier in article</option>
									<option value="artLineNum">Line in article</option>
								</select>
								<input type="checkbox" name="order" value="asc" filter="1" />
								Asc
							</nobr>
						</div>
					</div>
					<input type="hidden" name="predicate" value="" filter="1"
						class="invisible" />
					<input type="hidden" name="id" value="" filter="1" class="invisible"
						id="fixed-id" />
					<input type="hidden" name="molid" filter="1"
						value="{filter[@name='molid']/@value}" class="invisible" />
					<input type="hidden" name="value" filter="1"
						value="{filter[@name='value']/@value}" class="invisible" />
					<input type="hidden" name="canonicalValue" filter="1"
						value="{filter[@name='canonicalValue']/@value}" class="invisible" />
					<a href="javascript:void()" class="button-link"
						onclick="compoundBrowser.request(true); return false;">refresh</a>
					<a action="reset" class="button-link">reset</a>
				</div>
			</td>
			<td class="itunes-right">
				<div class="browserScope" style="position: relative;">
					<a action="showfilters" class="show-button invisible" title="Show filters panel">&lt;</a>
					<div class="upper-command-panel">
						<label class="no-trash">Basket</label>
						<a action="basketaction" type="addbasket"
							title="Add selected records matching current filters to a basket"
							class="no-trash">
							<img src="img/icons/add_to_set.gif" />
						</a>
						<a action="basketaction" type="removebasket"
							title="Remove selected records matching current filters from a basket"
							class="no-trash">
							<img src="img/icons/remove_from_set.gif" />
						</a>
						<a action="cleanbasketaction" type="cleanbasket" id="clean-basket"
							class="invisible" title="select duplicates within the basket">
							<img src="img/icons/duplicate.png" />
						</a>
						<label>Records</label>
						<a action="addselect" title="Select all records matching current filters"
							class="invisible">
							<img src="img/icons/select_all.gif" />
						</a>
						<a action="selectpage" title="Select records on currently visible page">
							<img src="img/icons/select_page.gif" />
						</a>
						<a action="removeselect" title="Unselect all records matching current filters">
							<img src="img/icons/unselect_all.gif" />
						</a>
						<a action="edit" title="Click to create new empty record" class="no-trash">
							<img src="img/icons/new.gif" />
						</a>
						<a href="batchupload/show.do?render-mode=full" target="_parent"
							title="Upload a set of properties/molecules from an excel or sdf file"
							class="no-trash">
							<img src="img/icons/batch_upload.gif" />
						</a>
						<a tab="Batch editing" href="batchedit/show.do?apply-basket=true"
							title="Batch edit selected molecules" class="no-trash">
							<img src="img/icons/edit.gif" />
						</a>
						<a action="deleteselected"
							title="Move selected records matching current filters to the trash"
							id="deleteselected">
							<img src="img/icons/delete.gif" />
						</a>
						<span class="trash invisible">
							<a action="restoreselected">
								<img src="img/icons/restore24.png" />
							</a>
						</span>
						<xsl:if test="session/user">
							<a href="epbrowser/show.do?trash=1&amp;render-mode=full" title="View your trash"
								target="_parent" class="no-trash" id="trash-button">
								<img src="img/icons/trash24.png" />
							</a>
						</xsl:if>
						<label>Tags</label>
						<a action="selecttag"
							title="Add a tag to selected compounds matching current filters">
							<img src="img/icons/tag.png" />
						</a>
						<a action="selecttag"
							title="Remove a tag from selected compounds matching current filters"
							removetag="1">
							<img src="img/icons/tag-orange.png" />
						</a>
					</div>
					<div class="pager-strip">
						<span>
							<b class="showed">none</b>
							of
							<b class="total">none</b>
						</span>
						<div id="pager" class="pgr">
						</div>
					</div>

					<div id="Browser">
					</div>
					<div class="pager-strip">
						<span>
							<b class="showed">none</b>
							of
							<b class="total">none</b>
						</span>
						<div id="pager" class="pgr">
						</div>
					</div>
				</div>
			</td>
		</tr>
	</table>
		
		
	</xsl:template>
</xsl:stylesheet>