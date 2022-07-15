<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="../../helper.xslt" />
	<xsl:include href="../model/inc/select-sets.xslt" />
	
	<xsl:template name="content">
	<style type="text/css">
		DIV.validation {margin-top: 5px; margin-left: 5px;}
		#template-configuration {padding: 10px; background-color: #FFE;}
	</style>
	<title>Model Builder</title>
		<table width="100%">
			<tr><td class="itunes-up silver">
				<p class="uploadable">
				Select model template, training set and descriptor-coefficient-sheet
				</p>				
				<h1><doc term="Creation+of+a+QSAR+model">Create a model</doc></h1>
				<p class="not-uploadable">
				Select the training and validation sets, the machine learning method and the validation protocol
				</p>
			</td></tr>
			<tr><td class="itunes-right">
				
				<form name="modelform" id="modelform" action="modelconfigurator/choosesubmit.do" method="post" enctype="multipart/form-data">
					<br/>
					<xsl:call-template name="select-sets"/>
					<xsl:if test="//ochem-model/method">
						<br/><input type="checkbox" name="skip-configuration" checked="true"/>&#160;Skip model configuration and use the predefined settings<br/>
						<div id="template-configuration">
						The model will be build used the predefined settings from a configuration template.
					</div><br/>
					</xsl:if>
					
					<div id="detailed-configuration">
					<b><doc term="Machine+learning+methods">Choose the learning method:</doc></b><br/>
					<div class="section">
					<i>Suggested modeling methods:</i>
					<xsl:for-each select="others/model-template">
						<xsl:if test="(@isSuggestedModel = 'true' and invisible != 'true') ">
							<div support="{@isSupportMultilearning}" uploadable="{@uploadable}"><input type="radio" name="template" id="{@id}" value="{@id}" support="{@isSupportMultilearning}" method="{@name}"/><label for="{@id}"><xsl:value-of select="displayedName"/></label></div>
						</xsl:if>
					</xsl:for-each><br/>
					<i>Methods under development:</i>
					<xsl:for-each select="others/model-template">
						<xsl:if test="(@isSuggestedModel = 'false' and invisible != 'true' and //user/ochem-labs = 'true')">
							<div support="{@isSupportMultilearning}" uploadable="{@uploadable}"><input type="radio" name="template" id="{@id}" value="{@id}" support="{@isSupportMultilearning}" method="{@name}"/><label for="{@id}"><xsl:value-of select="displayedName"/></label></div>
						</xsl:if>
					</xsl:for-each>
					</div>
					<br/>
					
					<div class="not-uploadable">
						<b>Model validation</b><br/>
						<div id="validation" class="section">
							Validation method: <select name="validation">
								<option value="no">No validation</option>
								<option value="cv" selected="true">N-Fold cross-validation</option>
								<option value="bagging">Bagging validation</option>
							</select>
							
							<div id="validation-cv" class="validation">
								Number of folds: <input type="text" name="cv-ensemble" class="small" value="5"/><br/>
								<input type="checkbox" name="cv-stratified"/> Stratified cross-validation (classification only)<a class="infolink" title="Will use under-sampling for imbalanced data: the number of samples in training set will be the same as for the smallest class."></a><br/>	
								<input type="checkbox" name="record-stratified"/> Treat each record as a new molecule<a class="infolink" title="The same molecule with different values can be in both training and validation sets simultaneously. This can contribute overfitted models. The option was added to mainly work with nano particles."></a><br/>	
							</div>
							<div id="validation-bagging" class="validation">
								Number of bagging models: <input type="text" name="bagging-ensemble" class="small" value="64"/><br/>
								<input type="checkbox" name="bagging-stratified"/> Stratified bagging (classification only)<br/>
								<input type="checkbox" name="usenuminstances"/> Use custom number of instances in training set for each bag<a class="infolink" title="Will use under-sampling for imbalanced data: the number of samples in training set will be the same as for the smallest class."></a><br/>	
								<div id="usenuminstances" class="invisible">Used fixed number of instances in bags: 
									<input type="text" name="bagging-instances" class="small" value="0"/>
								</div>
								<input type="checkbox" name="record-stratified"/> Treat each record as a new molecule<a class="infolink" title="The same molecule with different values can be in both training and validation sets simultaneously. This can contribute overfitted models. The option was added to mainly work with nano particles."></a><br/>	
							</div>
							<div class="warning validation" id="validation-no">
								Development of a model without validation may lead to incorrect estimation of the prediction accuracy of your model and incorrect assessment of its Applicability Domain.
							</div>
						</div>
						<div id="consensus-validation">
							Each member of the consensus model has its own validation method.
						</div>
						
						<br/>You can create a model from template: 
						<a href="modelconfigurator/importTemplate.do">import an XML model template</a> or  
						<a action="create_from_another_model">use another model as a template</a>
						
					</div>
					<br/>
					</div>
					
					<div class="formsubmit"> 
					<input type="button" action="submit" value="Next&gt;&gt;" name="next"/>
					</div>
					
				</form>
			</td></tr>
		</table>
		<script language="javascript" src="js/commons/actionable.js?ver=1.7.5"></script>
		<script language="javascript" src="js/commons/ajax-form.js"></script>
		<script language="javascript" src="js/commons/dynamic-select.js"></script>
		<script language="javascript" src="js/blocks/model-template.js?ver=1.8.9"></script>
		
		
		<script language="javascript">
			<xsl:if test="//ochem-model/method">
				form.defaultTemplate = "<xsl:value-of select="//ochem-model/method"/>";
			</xsl:if>
			
			<xsl:if test="//ochem-model/attachment/protocol/crossValidationConfiguration">
				$("select[name=validation]").val("cv").change();
				$("[name='cv-ensemble']").val(<xsl:value-of select="//ochem-model/attachment/protocol/crossValidationConfiguration/ensembleSize"/>);
				setCheckbox("cv-stratified", '<xsl:value-of select="//ochem-model/attachment/protocol/crossValidationConfiguration/validationType = 2"/>');
				setCheckbox("record-stratified", '<xsl:value-of select="//ochem-model/attachment/protocol/crossValidationConfiguration/mixtureValidation = RECORD"/>');
			</xsl:if>
			<xsl:if test="//ochem-model and not(//ochem-model/attachment/protocol/baggingConfiguration) and not(//ochem-model/attachment/protocol/crossValidationConfiguration)">
				$("select[name=validation]").val("no").change();
			</xsl:if>
			<xsl:if test="//ochem-model/attachment/protocol/baggingConfiguration">
				$("select[name=validation]").val("bagging").change();
				$("[name='bagging-ensemble']").val(<xsl:value-of select="//ochem-model/attachment/protocol/baggingConfiguration/ensembleSize"/>);
				setCheckbox("bagging-stratified", '<xsl:value-of select="//ochem-model/attachment/protocol/baggingConfiguration/validationType = 2"/>');
				$("[name='usenuminstances']").setChecked(<xsl:value-of select="//ochem-model/attachment/protocol/baggingConfiguration/numInstances &gt; 0"/>).change();
				setCheckbox("record-stratified", '<xsl:value-of select="//ochem-model/attachment/protocol/crossValidationConfiguration/validationType = RECORD"/>');
			</xsl:if>
			
		</script>
	</xsl:template>
</xsl:stylesheet>