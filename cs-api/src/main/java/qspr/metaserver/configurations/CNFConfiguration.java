/* Copyright (C) 2022 BIGCHEM GmbH <info@bigchem.de>
 *
 * Contact: info@bigchem.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License (AGPL)
 * as published by the Free Software Foundation; either version 3.0
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the Affero GNU General Public License for more details.
 *
 * You should have received a copy of the Affero GNU Lesser General Public License
 * along with this program; If not, see <https://www.gnu.org/licenses/>. 
 */

package qspr.metaserver.configurations;

import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.utils.NumericalValueStandardizer;

import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "cnf-configuration")
public class CNFConfiguration extends StandardNoDescriptorsConfiguration{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	// new params
	public Double dropout = 0.2;
	public TOKENIZER tokenizer = TOKENIZER.CHAR;
	public Integer FPDim = 32;
	public Integer nLayers = 8;
	public String FFNET_dim = "128,64";
	public Double rate = 0.001;
	public Integer batch = 32;
	public CNFTYPE type = CNFTYPE.V22;
	public NORMALIZE normalisation = NORMALIZE.LAYER;
	
	//public ERROR_FUNCTION errorFunction;
	//public Double nfp_alpha;
	//private Taxonomy tax;
	//public Boolean inversion[]; 
	//public String originalTaxonomy;

	//public Integer highway;
	//public Integer filterSize;

	public CNFConfiguration() {
		if(areDefault()) {
			nepochs = 1000;
			setAugmentations(10,10,false);
		}		
	}

	@Override
	public void setAugmentations(Integer training, Integer validation, boolean balance) {
		super.setAugmentations(training, validation, balance);
		if(training == -1) augmentation = -1;
	}

	public boolean isSupportRegression(){
		return true;
	}

	/*
	public ERROR_FUNCTION getErrorFunction(){
		return errorFunction == null? ERROR_FUNCTION.RMSE : errorFunction;
	}

	public enum ERROR_FUNCTION{
		RMSE, ENTROPY,  TAXONOMY, MULTILABELNN
	}
	 */

	@Override
	public boolean isSupportAugmentation(){
		return true;
	}

	public enum TOKENIZER{
		ATOM, CHAR, OLD
	};


	@Override
	public String getDefaultName() {
		return QSPRConstants.CNF;
	}

	public TOKENIZER getTokeniser() {
		return tokenizer == null? TOKENIZER.CHAR : tokenizer;
	}

	public String getTokeniserChar() {
		String t = getTokeniser().toString();
		return getTokeniser() != TOKENIZER.CHAR?" "+t.substring(0, 1):"";
	}


	@Override
	public String getInformativeName() {
		String s = super.getInformativeName() + (type == null ? "" : getType()) + " ";

		return s + 
				getTokeniserChar()
		//+ (highway != null ?" h="+highway :"")
		+ (getDropout() < 0 ?" D":"") + ((asPeptides != null && asPeptides)?" pept":"")
		;
	}


	public double getDropout(){
		return dropout == null? 0.3:dropout;
	}

	@Override
	public String toString(){

		return  " FPDim:" + FPDim + " nLayers:" + nLayers + " FFNET_dim:" + FFNET_dim + 
				(" nnf_type=" +  getType()) +  
				//(isUsesFilter() ?(" filterSize=" + filterSize):"") + 
				//(isUsesAlpha()?" alpha=" + getAlpha():"") +
				//				(originalTaxonomy == null ? "": (" " + originalTaxonomy + 
				//						(tax != null && tax.skeptTaxonomy != null ?" skipped: " + tax.skeptTaxonomy : ""))) +

				" batch=" + batch + " learning rate=" + rate + 

				//" loss: " + getErrorFunction() +

				//(highway != null?" highway=" + highway:"") +
				" " + getTokeniser() +
				(" dropout=" + getDropout()) + " " +super.toString();
	}

	public enum NORMALIZE{
		NONE, LAYER, BATCH
	}; 

	public enum CNFTYPE{
		V0, V11, V12, V21, V22
	};

	public boolean isUsesAlpha(){
		return type == CNFTYPE.V21 || type == CNFTYPE.V22;
	}

	public boolean isUsesFilter(){
		return type == CNFTYPE.V0 || type == CNFTYPE.V11 || type == CNFTYPE.V21;
	}

	public String getNormalize(){
		return normalisation.toString().toLowerCase();
	}

	public String getType(){
		return "v"+NumericalValueStandardizer.getSignificantDigits(getTypeInteger()/10.);
	}

	public Integer getTypeInteger(){
		return Integer.valueOf(type.toString().substring(1).toLowerCase());
	}

	/*
	@Override
	public void setTaxonomy(Taxonomy tax) {
		this.tax = tax;		
	}

	@Override
	public Taxonomy getTaxonomy() {
		return tax;
	}

	@Override
	public String getOriginalTaxonomyString() {
		return originalTaxonomy;
	}

	@Override
	public Boolean[] getInversions() {
		return inversion;
	}

	@Override
	public void setInversions(Boolean[] map) {
		inversion = map;
	}

	@Override
	public void setOriginalTaxonomy(String taxa) throws IOException {
		originalTaxonomy = SupportsTaxonomy.Taxonomy.validateOriginalSyntax(taxa);
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

	 */
	@Override
	public boolean isSupportConditions() {
		return true;
	}

	/*
	public Double getAlpha() {
		return nfp_alpha == null?0.5:nfp_alpha;
	}
	 */


}
