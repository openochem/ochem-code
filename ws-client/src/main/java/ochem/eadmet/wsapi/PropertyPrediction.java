/**
 * PropertyPrediction.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public class PropertyPrediction  implements java.io.Serializable {
    private java.lang.Double accuracy;

    private java.lang.Double dm;

    private java.lang.String error;

    private java.lang.Boolean insideAD;

    private java.lang.String moleculeId;

    private java.lang.String predictedValueString;

    private java.lang.String property;

    private java.lang.Double realValue;

    private java.lang.String realValueString;

    private java.lang.String unit;

    private java.lang.Double value;

    public PropertyPrediction() {
    }

    public PropertyPrediction(
           java.lang.Double accuracy,
           java.lang.Double dm,
           java.lang.String error,
           java.lang.Boolean insideAD,
           java.lang.String moleculeId,
           java.lang.String predictedValueString,
           java.lang.String property,
           java.lang.Double realValue,
           java.lang.String realValueString,
           java.lang.String unit,
           java.lang.Double value) {
           this.accuracy = accuracy;
           this.dm = dm;
           this.error = error;
           this.insideAD = insideAD;
           this.moleculeId = moleculeId;
           this.predictedValueString = predictedValueString;
           this.property = property;
           this.realValue = realValue;
           this.realValueString = realValueString;
           this.unit = unit;
           this.value = value;
    }


    /**
     * Gets the accuracy value for this PropertyPrediction.
     * 
     * @return accuracy
     */
    public java.lang.Double getAccuracy() {
        return accuracy;
    }


    /**
     * Sets the accuracy value for this PropertyPrediction.
     * 
     * @param accuracy
     */
    public void setAccuracy(java.lang.Double accuracy) {
        this.accuracy = accuracy;
    }


    /**
     * Gets the dm value for this PropertyPrediction.
     * 
     * @return dm
     */
    public java.lang.Double getDm() {
        return dm;
    }


    /**
     * Sets the dm value for this PropertyPrediction.
     * 
     * @param dm
     */
    public void setDm(java.lang.Double dm) {
        this.dm = dm;
    }


    /**
     * Gets the error value for this PropertyPrediction.
     * 
     * @return error
     */
    public java.lang.String getError() {
        return error;
    }


    /**
     * Sets the error value for this PropertyPrediction.
     * 
     * @param error
     */
    public void setError(java.lang.String error) {
        this.error = error;
    }


    /**
     * Gets the insideAD value for this PropertyPrediction.
     * 
     * @return insideAD
     */
    public java.lang.Boolean getInsideAD() {
        return insideAD;
    }


    /**
     * Sets the insideAD value for this PropertyPrediction.
     * 
     * @param insideAD
     */
    public void setInsideAD(java.lang.Boolean insideAD) {
        this.insideAD = insideAD;
    }


    /**
     * Gets the moleculeId value for this PropertyPrediction.
     * 
     * @return moleculeId
     */
    public java.lang.String getMoleculeId() {
        return moleculeId;
    }


    /**
     * Sets the moleculeId value for this PropertyPrediction.
     * 
     * @param moleculeId
     */
    public void setMoleculeId(java.lang.String moleculeId) {
        this.moleculeId = moleculeId;
    }


    /**
     * Gets the predictedValueString value for this PropertyPrediction.
     * 
     * @return predictedValueString
     */
    public java.lang.String getPredictedValueString() {
        return predictedValueString;
    }


    /**
     * Sets the predictedValueString value for this PropertyPrediction.
     * 
     * @param predictedValueString
     */
    public void setPredictedValueString(java.lang.String predictedValueString) {
        this.predictedValueString = predictedValueString;
    }


    /**
     * Gets the property value for this PropertyPrediction.
     * 
     * @return property
     */
    public java.lang.String getProperty() {
        return property;
    }


    /**
     * Sets the property value for this PropertyPrediction.
     * 
     * @param property
     */
    public void setProperty(java.lang.String property) {
        this.property = property;
    }


    /**
     * Gets the realValue value for this PropertyPrediction.
     * 
     * @return realValue
     */
    public java.lang.Double getRealValue() {
        return realValue;
    }


    /**
     * Sets the realValue value for this PropertyPrediction.
     * 
     * @param realValue
     */
    public void setRealValue(java.lang.Double realValue) {
        this.realValue = realValue;
    }


    /**
     * Gets the realValueString value for this PropertyPrediction.
     * 
     * @return realValueString
     */
    public java.lang.String getRealValueString() {
        return realValueString;
    }


    /**
     * Sets the realValueString value for this PropertyPrediction.
     * 
     * @param realValueString
     */
    public void setRealValueString(java.lang.String realValueString) {
        this.realValueString = realValueString;
    }


    /**
     * Gets the unit value for this PropertyPrediction.
     * 
     * @return unit
     */
    public java.lang.String getUnit() {
        return unit;
    }


    /**
     * Sets the unit value for this PropertyPrediction.
     * 
     * @param unit
     */
    public void setUnit(java.lang.String unit) {
        this.unit = unit;
    }


    /**
     * Gets the value value for this PropertyPrediction.
     * 
     * @return value
     */
    public java.lang.Double getValue() {
        return value;
    }


    /**
     * Sets the value value for this PropertyPrediction.
     * 
     * @param value
     */
    public void setValue(java.lang.Double value) {
        this.value = value;
    }

    private java.lang.Object __equalsCalc = null;
    public synchronized boolean equals(java.lang.Object obj) {
        if (!(obj instanceof PropertyPrediction)) return false;
        PropertyPrediction other = (PropertyPrediction) obj;
        if (obj == null) return false;
        if (this == obj) return true;
        if (__equalsCalc != null) {
            return (__equalsCalc == obj);
        }
        __equalsCalc = obj;
        boolean _equals;
        _equals = true && 
            ((this.accuracy==null && other.getAccuracy()==null) || 
             (this.accuracy!=null &&
              this.accuracy.equals(other.getAccuracy()))) &&
            ((this.dm==null && other.getDm()==null) || 
             (this.dm!=null &&
              this.dm.equals(other.getDm()))) &&
            ((this.error==null && other.getError()==null) || 
             (this.error!=null &&
              this.error.equals(other.getError()))) &&
            ((this.insideAD==null && other.getInsideAD()==null) || 
             (this.insideAD!=null &&
              this.insideAD.equals(other.getInsideAD()))) &&
            ((this.moleculeId==null && other.getMoleculeId()==null) || 
             (this.moleculeId!=null &&
              this.moleculeId.equals(other.getMoleculeId()))) &&
            ((this.predictedValueString==null && other.getPredictedValueString()==null) || 
             (this.predictedValueString!=null &&
              this.predictedValueString.equals(other.getPredictedValueString()))) &&
            ((this.property==null && other.getProperty()==null) || 
             (this.property!=null &&
              this.property.equals(other.getProperty()))) &&
            ((this.realValue==null && other.getRealValue()==null) || 
             (this.realValue!=null &&
              this.realValue.equals(other.getRealValue()))) &&
            ((this.realValueString==null && other.getRealValueString()==null) || 
             (this.realValueString!=null &&
              this.realValueString.equals(other.getRealValueString()))) &&
            ((this.unit==null && other.getUnit()==null) || 
             (this.unit!=null &&
              this.unit.equals(other.getUnit()))) &&
            ((this.value==null && other.getValue()==null) || 
             (this.value!=null &&
              this.value.equals(other.getValue())));
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
        if (getAccuracy() != null) {
            _hashCode += getAccuracy().hashCode();
        }
        if (getDm() != null) {
            _hashCode += getDm().hashCode();
        }
        if (getError() != null) {
            _hashCode += getError().hashCode();
        }
        if (getInsideAD() != null) {
            _hashCode += getInsideAD().hashCode();
        }
        if (getMoleculeId() != null) {
            _hashCode += getMoleculeId().hashCode();
        }
        if (getPredictedValueString() != null) {
            _hashCode += getPredictedValueString().hashCode();
        }
        if (getProperty() != null) {
            _hashCode += getProperty().hashCode();
        }
        if (getRealValue() != null) {
            _hashCode += getRealValue().hashCode();
        }
        if (getRealValueString() != null) {
            _hashCode += getRealValueString().hashCode();
        }
        if (getUnit() != null) {
            _hashCode += getUnit().hashCode();
        }
        if (getValue() != null) {
            _hashCode += getValue().hashCode();
        }
        __hashCodeCalc = false;
        return _hashCode;
    }

    // Type metadata
    private static org.apache.axis.description.TypeDesc typeDesc =
        new org.apache.axis.description.TypeDesc(PropertyPrediction.class, true);

    static {
        typeDesc.setXmlType(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "PropertyPrediction"));
        org.apache.axis.description.ElementDesc elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("accuracy");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "accuracy"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(false);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("dm");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "dm"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(false);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("error");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "error"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("insideAD");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "insideAD"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "boolean"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("moleculeId");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "moleculeId"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("predictedValueString");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "predictedValueString"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("property");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "property"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("realValue");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "realValue"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("realValueString");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "realValueString"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("unit");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "unit"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("value");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "value"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "double"));
        elemField.setMinOccurs(0);
        elemField.setNillable(false);
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
