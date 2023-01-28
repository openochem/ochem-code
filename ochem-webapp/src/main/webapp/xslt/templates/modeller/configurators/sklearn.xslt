<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:include href="configurators-common.xslt" />
	
	<xsl:template name="configuration-content">
		<title>Model builder - SKLearn</title>
		<h1>Select SKLearn method</h1>
		<table class="configuration">
			<tr>
				<td>Model:</td>
				<td>
				<select name="method">
					<option value="ALL_CV" selected="1">ALL-CV (automatic selection of the best cross-validated methods)</option>
					<option value="ALL">ALL - (automatic selection of the best methods)</option>
					<option value="ADA_BOOST">ADA_BOOST</option>
					<option value="ADA_BOOST_TREES">ADA_BOOST+TREES</option>
					<option value="ARD">ARD</option>
					<option value="BAGGING">BAGGING</option>
					<option value="BAGGING_TREES">BAGGING+TREES</option>
					<option value="BAYESIAN_RIDGE">BAYESIAN_RIDGE</option>
					<option value="CAT_BOOST">CatBoost</option>
					<option value="GBM">LightGBM</option>
					<option value="ELASTIC_NET_CV">ELASTIC_NET_CV</option>
					<option value="ELASTIC_NET">ELASTIC_NET</option>
					<option value="EXTRA_TREES">EXTRA_TREES</option>
					<option value="GAUSSIAN_PROCESS">GAUSSIAN_PROCESS</option>
					<option value="GRADIENT_BOOSTING">GRADIENT_BOOSTING</option>
					<option value="HUBER">HUBER</option>
					<option value="KERNEL_RIDGE">KERNEL_RIDGE</option>
					<option value="KERNEL_RIDGE_HYPER">KERNEL_RIDGE_HYPER</option>
					<option value="K_NEIGHBORS">K_NEIGHBORS</option>
					<option value="LASSO_LARS_CV">LASSO_LARS_CV</option>
					<option value="LASSO_LARS_IC">LASSO_LARS_IC</option>
					<option value="LINEAR">LINEAR</option>
					<option value="LOGISTIC_CV">LOGISTIC_CV</option>
					<option value="ORTHOGONAL_MATCHING_PURSUIT_CV">ORTHOGONAL_MATCHING_PURSUIT_CV</option>
					<option value="RANDOM_FOREST">RANDOM_FOREST</option>
					<option value="RANSAC_TREES">RANSAC+TREES</option>
					<option value="RIDGE_CV">RIDGE_CV</option>
					<option value="SVM">SVM</option>
				</select></td>
			</tr>
		</table>

		<script language="javascript">
				<xsl:if test="//method = 'SKLEARN'">
				setValue("method", '<xsl:value-of select="//attachment/configuration/modelConfiguration/method"/>');
				setValue("epochs", '<xsl:value-of select="//attachment/configuration/modelConfiguration/epochs"/>');
				setValue("ratio", '<xsl:value-of select="//attachment/configuration/modelConfiguration/ratio"/>');
				setValue("learning_rate", '<xsl:value-of select="//attachment/configuration/modelConfiguration/learning_rate"/>');
				setValue("dropout", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dropout"/>');
				setValue("n_hidden", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_hidden"/>');
				setValue("dense_layer_size", '<xsl:value-of select="//attachment/configuration/modelConfiguration/dense_layer_size"/>');
				setValue("graph_conv_layers", '<xsl:value-of select="//attachment/configuration/modelConfiguration/graph_conv_layers"/>');
				setValue("M", '<xsl:value-of select="//attachment/configuration/modelConfiguration/M"/>');
				setValue("T", '<xsl:value-of select="//attachment/configuration/modelConfiguration/T"/>');
				setValue("n_embedding", '<xsl:value-of select="//attachment/configuration/modelConfiguration/n_embedding"/>');
			</xsl:if>
		</script>
	</xsl:template>
</xsl:stylesheet>

	