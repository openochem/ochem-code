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

package qspr.services;

public class ModelSummary {

	private Double rmse;
	private Double rmseConfidence;
	private Double r2;
	private Double r2Confidence;
	private Double q2;
	private Double q2Confidence;
	private Double mae;
	private Double maeConfidence;

	private Double accuracy;
	private Double accuracyConfidence;
	private Double accuracyBalanced;
	private Double accuracyBalancedConfidence;
	private Double mcc;
	private Double mccConfidence;
	private Integer tp;
	private Integer tn;
	private Integer fp;
	private Integer fn;
	private Integer n;

	private String propertyName;
	private String validationSetName;
	
	private String modelName;

	public Double getRmse() {
		return rmse;
	}

	public void setRmse(Double rmse) {
		this.rmse = rmse;
	}

	public Double getRmseConfidence() {
		return rmseConfidence;
	}

	public void setRmseConfidence(Double rmseConfidence) {
		this.rmseConfidence = rmseConfidence;
	}

	public Double getR2() {
		return r2;
	}

	public void setR2(Double r2) {
		this.r2 = r2;
	}

	public Double getR2Confidence() {
		return r2Confidence;
	}

	public void setR2Confidence(Double r2Confidence) {
		this.r2Confidence = r2Confidence;
	}

	public Double getQ2() {
		return q2;
	}

	public void setQ2(Double q2) {
		this.q2 = q2;
	}

	public Double getQ2Confidence() {
		return q2Confidence;
	}

	public void setQ2Confidence(Double q2Confidence) {
		this.q2Confidence = q2Confidence;
	}

	public Double getMae() {
		return mae;
	}

	public void setMae(Double mae) {
		this.mae = mae;
	}

	public Double getMaeConfidence() {
		return maeConfidence;
	}

	public void setMaeConfidence(Double maeConfidence) {
		this.maeConfidence = maeConfidence;
	}

	public Double getAccuracy() {
		return accuracy;
	}

	public void setAccuracy(Double accuracy) {
		this.accuracy = accuracy;
	}

	public Double getAccuracyConfidence() {
		return accuracyConfidence;
	}

	public void setAccuracyConfidence(Double accuracyConfidence) {
		this.accuracyConfidence = accuracyConfidence;
	}

	public Double getMcc() {
		return mcc;
	}

	public void setMcc(Double mcc) {
		this.mcc = mcc;
	}

	public Double getMccConfidence() {
		return mccConfidence;
	}

	public void setMccConfidence(Double mccConfidence) {
		this.mccConfidence = mccConfidence;
	}

	public Double getAccuracyBalanced() {
		return accuracyBalanced;
	}

	public void setAccuracyBalanced(Double accuracyBalanced) {
		this.accuracyBalanced = accuracyBalanced;
	}

	public Double getAccuracyBalancedConfidence() {
		return accuracyBalancedConfidence;
	}

	public void setAccuracyBalancedConfidence(Double accuracyBalancedConfidence) {
		this.accuracyBalancedConfidence = accuracyBalancedConfidence;
	}

	public String getPropertyName() {
		return propertyName;
	}

	public void setPropertyName(String propertyName) {
		this.propertyName = propertyName;
	}

	public String getValidationSetName() {
		return validationSetName;
	}

	public void setValidationSetName(String validationSetName) {
		this.validationSetName = validationSetName;
	}

	public Integer getTP() {
		return tp;
	}

	public void setTP(Integer tp) {
		this.tp = tp;
	}

	public Integer getTN() {
		return tn;
	}

	public void setTN(Integer tn) {
		this.tn = tn;
	}

	public Integer getFP() {
		return fp;
	}

	public void setFP(Integer fp) {
		this.fp = fp;
	}

	public void setFN(Integer fn) {
		this.fn = fn;
	}

	public Integer getFN() {
		return fn;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("== Model Summary ==");
		sb.append("\n L Property Name   : ");
		sb.append(getPropertyName());
		sb.append("\n L R Squared       : ");
		sb.append(getR2());
		sb.append("\n L Q Squared       : ");
		sb.append(getQ2());
		sb.append("\n L Validation Set  : ");
		sb.append(getValidationSetName());
		sb.append("\n L RMSE            : ");
		sb.append(getRmse());
		sb.append("\n L RMSE Confidence : ");
		sb.append(getRmseConfidence());
		sb.append("\n L MAE             : ");
		sb.append(getMae());
		sb.append("\n L MAE Confidence  : ");
		sb.append(getMaeConfidence());
		sb.append("\n L Accuracy (Bal)  : ");
		sb.append(getAccuracyBalanced());
		sb.append("\n L Accuracy        : ");
		sb.append(getAccuracy());
		sb.append("\n L Acc. Bal. Conf. : ");
		sb.append(getAccuracyBalancedConfidence());
		sb.append("\n L N : ");
		sb.append(getN());
		return sb.toString();
	}

	public Integer getN() {
		return n;
	}

	public void setN(Integer n) {
		this.n = n;
	}

	public String getModelName() {
		return modelName;
	}

	public void setModelName(String modelName) {
		this.modelName = modelName;
	}

}
