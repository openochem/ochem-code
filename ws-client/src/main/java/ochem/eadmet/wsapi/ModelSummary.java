/**
 * ModelSummary.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public class ModelSummary  implements java.io.Serializable {
    private java.lang.Integer FN;

    private java.lang.Integer FP;

    private java.lang.Integer TN;

    private java.lang.Integer TP;

    private java.lang.Double accuracy;

    private java.lang.Double accuracyBalanced;

    private java.lang.Double accuracyBalancedConfidence;

    private java.lang.Double accuracyConfidence;

    private java.lang.Double mae;

    private java.lang.Double maeConfidence;

    private java.lang.Double mcc;

    private java.lang.Double mccConfidence;

    private java.lang.String modelName;

    private java.lang.Integer n;

    private java.lang.String propertyName;

    private java.lang.Double q2;

    private java.lang.Double q2Confidence;

    private java.lang.Double r2;

    private java.lang.Double r2Confidence;

    private java.lang.Double rmse;

    private java.lang.Double rmseConfidence;

    private java.lang.String validationSetName;

    public ModelSummary() {
    }

    public ModelSummary(
           java.lang.Integer FN,
           java.lang.Integer FP,
           java.lang.Integer TN,
           java.lang.Integer TP,
           java.lang.Double accuracy,
           java.lang.Double accuracyBalanced,
           java.lang.Double accuracyBalancedConfidence,
           java.lang.Double accuracyConfidence,
           java.lang.Double mae,
           java.lang.Double maeConfidence,
           java.lang.Double mcc,
           java.lang.Double mccConfidence,
           java.lang.String modelName,
           java.lang.Integer n,
           java.lang.String propertyName,
           java.lang.Double q2,
           java.lang.Double q2Confidence,
           java.lang.Double r2,
           java.lang.Double r2Confidence,
           java.lang.Double rmse,
           java.lang.Double rmseConfidence,
           java.lang.String validationSetName) {
           this.FN = FN;
           this.FP = FP;
           this.TN = TN;
           this.TP = TP;
           this.accuracy = accuracy;
           this.accuracyBalanced = accuracyBalanced;
           this.accuracyBalancedConfidence = accuracyBalancedConfidence;
           this.accuracyConfidence = accuracyConfidence;
           this.mae = mae;
           this.maeConfidence = maeConfidence;
           this.mcc = mcc;
           this.mccConfidence = mccConfidence;
           this.modelName = modelName;
           this.n = n;
           this.propertyName = propertyName;
           this.q2 = q2;
           this.q2Confidence = q2Confidence;
           this.r2 = r2;
           this.r2Confidence = r2Confidence;
           this.rmse = rmse;
           this.rmseConfidence = rmseConfidence;
           this.validationSetName = validationSetName;
    }


    /**
     * Gets the FN value for this ModelSummary.
     * 
     * @return FN
     */
    public java.lang.Integer getFN() {
        return FN;
    }


    /**
     * Sets the FN value for this ModelSummary.
     * 
     * @param FN
     */
    public void setFN(java.lang.Integer FN) {
        this.FN = FN;
    }


    /**
     * Gets the FP value for this ModelSummary.
     * 
     * @return FP
     */
    public java.lang.Integer getFP() {
        return FP;
    }


    /**
     * Sets the FP value for this ModelSummary.
     * 
     * @param FP
     */
    public void setFP(java.lang.Integer FP) {
        this.FP = FP;
    }


    /**
     * Gets the TN value for this ModelSummary.
     * 
     * @return TN
     */
    public java.lang.Integer getTN() {
        return TN;
    }


    /**
     * Sets the TN value for this ModelSummary.
     * 
     * @param TN
     */
    public void setTN(java.lang.Integer TN) {
        this.TN = TN;
    }


    /**
     * Gets the TP value for this ModelSummary.
     * 
     * @return TP
     */
    public java.lang.Integer getTP() {
        return TP;
    }


    /**
     * Sets the TP value for this ModelSummary.
     * 
     * @param TP
     */
    public void setTP(java.lang.Integer TP) {
        this.TP = TP;
    }


    /**
     * Gets the accuracy value for this ModelSummary.
     * 
     * @return accuracy
     */
    public java.lang.Double getAccuracy() {
        return accuracy;
    }


    /**
     * Sets the accuracy value for this ModelSummary.
     * 
     * @param accuracy
     */
    public void setAccuracy(java.lang.Double accuracy) {
        this.accuracy = accuracy;
    }


    /**
     * Gets the accuracyBalanced value for this ModelSummary.
     * 
     * @return accuracyBalanced
     */
    public java.lang.Double getAccuracyBalanced() {
        return accuracyBalanced;
    }


    /**
     * Sets the accuracyBalanced value for this ModelSummary.
     * 
     * @param accuracyBalanced
     */
    public void setAccuracyBalanced(java.lang.Double accuracyBalanced) {
        this.accuracyBalanced = accuracyBalanced;
    }


    /**
     * Gets the accuracyBalancedConfidence value for this ModelSummary.
     * 
     * @return accuracyBalancedConfidence
     */
    public java.lang.Double getAccuracyBalancedConfidence() {
        return accuracyBalancedConfidence;
    }


    /**
     * Sets the accuracyBalancedConfidence value for this ModelSummary.
     * 
     * @param accuracyBalancedConfidence
     */
    public void setAccuracyBalancedConfidence(java.lang.Double accuracyBalancedConfidence) {
        this.accuracyBalancedConfidence = accuracyBalancedConfidence;
    }


    /**
     * Gets the accuracyConfidence value for this ModelSummary.
     * 
     * @return accuracyConfidence
     */
    public java.lang.Double getAccuracyConfidence() {
        return accuracyConfidence;
    }


    /**
     * Sets the accuracyConfidence value for this ModelSummary.
     * 
     * @param accuracyConfidence
     */
    public void setAccuracyConfidence(java.lang.Double accuracyConfidence) {
        this.accuracyConfidence = accuracyConfidence;
    }


    /**
     * Gets the mae value for this ModelSummary.
     * 
     * @return mae
     */
    public java.lang.Double getMae() {
        return mae;
    }


    /**
     * Sets the mae value for this ModelSummary.
     * 
     * @param mae
     */
    public void setMae(java.lang.Double mae) {
        this.mae = mae;
    }


    /**
     * Gets the maeConfidence value for this ModelSummary.
     * 
     * @return maeConfidence
     */
    public java.lang.Double getMaeConfidence() {
        return maeConfidence;
    }


    /**
     * Sets the maeConfidence value for this ModelSummary.
     * 
     * @param maeConfidence
     */
    public void setMaeConfidence(java.lang.Double maeConfidence) {
        this.maeConfidence = maeConfidence;
    }


    /**
     * Gets the mcc value for this ModelSummary.
     * 
     * @return mcc
     */
    public java.lang.Double getMcc() {
        return mcc;
    }


    /**
     * Sets the mcc value for this ModelSummary.
     * 
     * @param mcc
     */
    public void setMcc(java.lang.Double mcc) {
        this.mcc = mcc;
    }


    /**
     * Gets the mccConfidence value for this ModelSummary.
     * 
     * @return mccConfidence
     */
    public java.lang.Double getMccConfidence() {
        return mccConfidence;
    }


    /**
     * Sets the mccConfidence value for this ModelSummary.
     * 
     * @param mccConfidence
     */
    public void setMccConfidence(java.lang.Double mccConfidence) {
        this.mccConfidence = mccConfidence;
    }


    /**
     * Gets the modelName value for this ModelSummary.
     * 
     * @return modelName
     */
    public java.lang.String getModelName() {
        return modelName;
    }


    /**
     * Sets the modelName value for this ModelSummary.
     * 
     * @param modelName
     */
    public void setModelName(java.lang.String modelName) {
        this.modelName = modelName;
    }


    /**
     * Gets the n value for this ModelSummary.
     * 
     * @return n
     */
    public java.lang.Integer getN() {
        return n;
    }


    /**
     * Sets the n value for this ModelSummary.
     * 
     * @param n
     */
    public void setN(java.lang.Integer n) {
        this.n = n;
    }


    /**
     * Gets the propertyName value for this ModelSummary.
     * 
     * @return propertyName
     */
    public java.lang.String getPropertyName() {
        return propertyName;
    }


    /**
     * Sets the propertyName value for this ModelSummary.
     * 
     * @param propertyName
     */
    public void setPropertyName(java.lang.String propertyName) {
        this.propertyName = propertyName;
    }


    /**
     * Gets the q2 value for this ModelSummary.
     * 
     * @return q2
     */
    public java.lang.Double getQ2() {
        return q2;
    }


    /**
     * Sets the q2 value for this ModelSummary.
     * 
     * @param q2
     */
    public void setQ2(java.lang.Double q2) {
        this.q2 = q2;
    }


    /**
     * Gets the q2Confidence value for this ModelSummary.
     * 
     * @return q2Confidence
     */
    public java.lang.Double getQ2Confidence() {
        return q2Confidence;
    }


    /**
     * Sets the q2Confidence value for this ModelSummary.
     * 
     * @param q2Confidence
     */
    public void setQ2Confidence(java.lang.Double q2Confidence) {
        this.q2Confidence = q2Confidence;
    }


    /**
     * Gets the r2 value for this ModelSummary.
     * 
     * @return r2
     */
    public java.lang.Double getR2() {
        return r2;
    }


    /**
     * Sets the r2 value for this ModelSummary.
     * 
     * @param r2
     */
    public void setR2(java.lang.Double r2) {
        this.r2 = r2;
    }


    /**
     * Gets the r2Confidence value for this ModelSummary.
     * 
     * @return r2Confidence
     */
    public java.lang.Double getR2Confidence() {
        return r2Confidence;
    }


    /**
     * Sets the r2Confidence value for this ModelSummary.
     * 
     * @param r2Confidence
     */
    public void setR2Confidence(java.lang.Double r2Confidence) {
        this.r2Confidence = r2Confidence;
    }


    /**
     * Gets the rmse value for this ModelSummary.
     * 
     * @return rmse
     */
    public java.lang.Double getRmse() {
        return rmse;
    }


    /**
     * Sets the rmse value for this ModelSummary.
     * 
     * @param rmse
     */
    public void setRmse(java.lang.Double rmse) {
        this.rmse = rmse;
    }


    /**
     * Gets the rmseConfidence value for this ModelSummary.
     * 
     * @return rmseConfidence
     */
    public java.lang.Double getRmseConfidence() {
        return rmseConfidence;
    }


    /**
     * Sets the rmseConfidence value for this ModelSummary.
     * 
     * @param rmseConfidence
     */
    public void setRmseConfidence(java.lang.Double rmseConfidence) {
        this.rmseConfidence = rmseConfidence;
    }


    /**
     * Gets the validationSetName value for this ModelSummary.
     * 
     * @return validationSetName
     */
    public java.lang.String getValidationSetName() {
        return validationSetName;
    }


    /**
     * Sets the validationSetName value for this ModelSummary.
     * 
     * @param validationSetName
     */
    public void setValidationSetName(java.lang.String validationSetName) {
        this.validationSetName = validationSetName;
    }

    private java.lang.Object __equalsCalc = null;
    public synchronized boolean equals(java.lang.Object obj) {
        if (!(obj instanceof ModelSummary)) return false;
        ModelSummary other = (ModelSummary) obj;
        if (obj == null) return false;
        if (this == obj) return true;
        if (__equalsCalc != null) {
            return (__equalsCalc == obj);
        }
        __equalsCalc = obj;
        boolean _equals;
        _equals = true && 
            ((this.FN==null && other.getFN()==null) || 
             (this.FN!=null &&
              this.FN.equals(other.getFN()))) &&
            ((this.FP==null && other.getFP()==null) || 
             (this.FP!=null &&
              this.FP.equals(other.getFP()))) &&
            ((this.TN==null && other.getTN()==null) || 
             (this.TN!=null &&
              this.TN.equals(other.getTN()))) &&
            ((this.TP==null && other.getTP()==null) || 
             (this.TP!=null &&
              this.TP.equals(other.getTP()))) &&
            ((this.accuracy==null && other.getAccuracy()==null) || 
             (this.accuracy!=null &&
              this.accuracy.equals(other.getAccuracy()))) &&
            ((this.accuracyBalanced==null && other.getAccuracyBalanced()==null) || 
             (this.accuracyBalanced!=null &&
              this.accuracyBalanced.equals(other.getAccuracyBalanced()))) &&
            ((this.accuracyBalancedConfidence==null && other.getAccuracyBalancedConfidence()==null) || 
             (this.accuracyBalancedConfidence!=null &&
              this.accuracyBalancedConfidence.equals(other.getAccuracyBalancedConfidence()))) &&
            ((this.accuracyConfidence==null && other.getAccuracyConfidence()==null) || 
             (this.accuracyConfidence!=null &&
              this.accuracyConfidence.equals(other.getAccuracyConfidence()))) &&
            ((this.mae==null && other.getMae()==null) || 
             (this.mae!=null &&
              this.mae.equals(other.getMae()))) &&
            ((this.maeConfidence==null && other.getMaeConfidence()==null) || 
             (this.maeConfidence!=null &&
              this.maeConfidence.equals(other.getMaeConfidence()))) &&
            ((this.mcc==null && other.getMcc()==null) || 
             (this.mcc!=null &&
              this.mcc.equals(other.getMcc()))) &&
            ((this.mccConfidence==null && other.getMccConfidence()==null) || 
             (this.mccConfidence!=null &&
              this.mccConfidence.equals(other.getMccConfidence()))) &&
            ((this.modelName==null && other.getModelName()==null) || 
             (this.modelName!=null &&
              this.modelName.equals(other.getModelName()))) &&
            ((this.n==null && other.getN()==null) || 
             (this.n!=null &&
              this.n.equals(other.getN()))) &&
            ((this.propertyName==null && other.getPropertyName()==null) || 
             (this.propertyName!=null &&
              this.propertyName.equals(other.getPropertyName()))) &&
            ((this.q2==null && other.getQ2()==null) || 
             (this.q2!=null &&
              this.q2.equals(other.getQ2()))) &&
            ((this.q2Confidence==null && other.getQ2Confidence()==null) || 
             (this.q2Confidence!=null &&
              this.q2Confidence.equals(other.getQ2Confidence()))) &&
            ((this.r2==null && other.getR2()==null) || 
             (this.r2!=null &&
              this.r2.equals(other.getR2()))) &&
            ((this.r2Confidence==null && other.getR2Confidence()==null) || 
             (this.r2Confidence!=null &&
              this.r2Confidence.equals(other.getR2Confidence()))) &&
            ((this.rmse==null && other.getRmse()==null) || 
             (this.rmse!=null &&
              this.rmse.equals(other.getRmse()))) &&
            ((this.rmseConfidence==null && other.getRmseConfidence()==null) || 
             (this.rmseConfidence!=null &&
              this.rmseConfidence.equals(other.getRmseConfidence()))) &&
            ((this.validationSetName==null && other.getValidationSetName()==null) || 
             (this.validationSetName!=null &&
              this.validationSetName.equals(other.getValidationSetName())));
        __equalsCalc = null;
        return _equals;
    }

    private boolean __hashCodeCalc = false;
    public synchronized int hashCode() {
        if (__hashCodeCalc) {
            return 0;
        }
        __hashCodeCalc = true;
        int _hashCode = 1;
        if (getFN() != null) {
            _hashCode += getFN().hashCode();
        }
        if (getFP() != null) {
            _hashCode += getFP().hashCode();
        }
        if (getTN() != null) {
            _hashCode += getTN().hashCode();
        }
        if (getTP() != null) {
            _hashCode += getTP().hashCode();
        }
        if (getAccuracy() != null) {
            _hashCode += getAccuracy().hashCode();
        }
        if (getAccuracyBalanced() != null) {
            _hashCode += getAccuracyBalanced().hashCode();
        }
        if (getAccuracyBalancedConfidence() != null) {
            _hashCode += getAccuracyBalancedConfidence().hashCode();
        }
        if (getAccuracyConfidence() != null) {
            _hashCode += getAccuracyConfidence().hashCode();
        }
        if (getMae() != null) {
            _hashCode += getMae().hashCode();
        }
        if (getMaeConfidence() != null) {
            _hashCode += getMaeConfidence().hashCode();
        }
        if (getMcc() != null) {
            _hashCode += getMcc().hashCode();
        }
        if (getMccConfidence() != null) {
            _hashCode += getMccConfidence().hashCode();
        }
        if (getModelName() != null) {
            _hashCode += getModelName().hashCode();
        }
        if (getN() != null) {
            _hashCode += getN().hashCode();
        }
        if (getPropertyName() != null) {
            _hashCode += getPropertyName().hashCode();
        }
        if (getQ2() != null) {
            _hashCode += getQ2().hashCode();
        }
        if (getQ2Confidence() != null) {
            _hashCode += getQ2Confidence().hashCode();
        }
        if (getR2() != null) {
            _hashCode += getR2().hashCode();
        }
        if (getR2Confidence() != null) {
            _hashCode += getR2Confidence().hashCode();
        }
        if (getRmse() != null) {
            _hashCode += getRmse().hashCode();
        }
        if (getRmseConfidence() != null) {
            _hashCode += getRmseConfidence().hashCode();
        }
        if (getValidationSetName() != null) {
            _hashCode += getValidationSetName().hashCode();
        }
        __hashCodeCalc = false;
        return _hashCode;
    }

    // Type metadata
    private static org.apache.axis.description.TypeDesc typeDesc =
        new org.apache.axis.description.TypeDesc(ModelSummary.class, true);

    static {
        typeDesc.setXmlType(new javax.xml.namespace.QName("http://services.qspr/xsd", "ModelSummary"));
        org.apache.axis.description.ElementDesc elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("FN");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "FN"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "int"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("FP");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "FP"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "int"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("TN");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "TN"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "int"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("TP");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "TP"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "int"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("accuracy");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "accuracy"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("accuracyBalanced");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "accuracyBalanced"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("accuracyBalancedConfidence");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "accuracyBalancedConfidence"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("accuracyConfidence");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "accuracyConfidence"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("mae");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "mae"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("maeConfidence");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "maeConfidence"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("mcc");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "mcc"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("mccConfidence");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "mccConfidence"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("modelName");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "modelName"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("n");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "n"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "int"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("propertyName");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "propertyName"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("q2");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "q2"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("q2Confidence");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "q2Confidence"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("r2");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "r2"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("r2Confidence");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "r2Confidence"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("rmse");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "rmse"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("rmseConfidence");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "rmseConfidence"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("validationSetName");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "validationSetName"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
    }

    /**
     * Return type metadata object
     */
    public static org.apache.axis.description.TypeDesc getTypeDesc() {
        return typeDesc;
    }

    /**
     * Get Custom Serializer
     */
    public static org.apache.axis.encoding.Serializer getSerializer(
           java.lang.String mechType, 
           java.lang.Class _javaType,  
           javax.xml.namespace.QName _xmlType) {
        return 
          new  org.apache.axis.encoding.ser.BeanSerializer(
            _javaType, _xmlType, typeDesc);
    }

    /**
     * Get Custom Deserializer
     */
    public static org.apache.axis.encoding.Deserializer getDeserializer(
           java.lang.String mechType, 
           java.lang.Class _javaType,  
           javax.xml.namespace.QName _xmlType) {
        return 
          new  org.apache.axis.encoding.ser.BeanDeserializer(
            _javaType, _xmlType, typeDesc);
    }

}
