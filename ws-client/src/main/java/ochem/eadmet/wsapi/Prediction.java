/**
 * Prediction.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public class Prediction  implements java.io.Serializable {
    private java.lang.Long depictionID;

    private java.lang.String error;

    private java.lang.String inChIKey;

    private java.lang.Integer moleculeID;

    private ochem.eadmet.wsapi.PropertyPrediction[] predictions;

    private java.lang.String smiles;

    public Prediction() {
    }

    public Prediction(
           java.lang.Long depictionID,
           java.lang.String error,
           java.lang.String inChIKey,
           java.lang.Integer moleculeID,
           ochem.eadmet.wsapi.PropertyPrediction[] predictions,
           java.lang.String smiles) {
           this.depictionID = depictionID;
           this.error = error;
           this.inChIKey = inChIKey;
           this.moleculeID = moleculeID;
           this.predictions = predictions;
           this.smiles = smiles;
    }


    /**
     * Gets the depictionID value for this Prediction.
     * 
     * @return depictionID
     */
    public java.lang.Long getDepictionID() {
        return depictionID;
    }


    /**
     * Sets the depictionID value for this Prediction.
     * 
     * @param depictionID
     */
    public void setDepictionID(java.lang.Long depictionID) {
        this.depictionID = depictionID;
    }


    /**
     * Gets the error value for this Prediction.
     * 
     * @return error
     */
    public java.lang.String getError() {
        return error;
    }


    /**
     * Sets the error value for this Prediction.
     * 
     * @param error
     */
    public void setError(java.lang.String error) {
        this.error = error;
    }


    /**
     * Gets the inChIKey value for this Prediction.
     * 
     * @return inChIKey
     */
    public java.lang.String getInChIKey() {
        return inChIKey;
    }


    /**
     * Sets the inChIKey value for this Prediction.
     * 
     * @param inChIKey
     */
    public void setInChIKey(java.lang.String inChIKey) {
        this.inChIKey = inChIKey;
    }


    /**
     * Gets the moleculeID value for this Prediction.
     * 
     * @return moleculeID
     */
    public java.lang.Integer getMoleculeID() {
        return moleculeID;
    }


    /**
     * Sets the moleculeID value for this Prediction.
     * 
     * @param moleculeID
     */
    public void setMoleculeID(java.lang.Integer moleculeID) {
        this.moleculeID = moleculeID;
    }


    /**
     * Gets the predictions value for this Prediction.
     * 
     * @return predictions
     */
    public ochem.eadmet.wsapi.PropertyPrediction[] getPredictions() {
        return predictions;
    }


    /**
     * Sets the predictions value for this Prediction.
     * 
     * @param predictions
     */
    public void setPredictions(ochem.eadmet.wsapi.PropertyPrediction[] predictions) {
        this.predictions = predictions;
    }

    public ochem.eadmet.wsapi.PropertyPrediction getPredictions(int i) {
        return this.predictions[i];
    }

    public void setPredictions(int i, ochem.eadmet.wsapi.PropertyPrediction _value) {
        this.predictions[i] = _value;
    }


    /**
     * Gets the smiles value for this Prediction.
     * 
     * @return smiles
     */
    public java.lang.String getSmiles() {
        return smiles;
    }


    /**
     * Sets the smiles value for this Prediction.
     * 
     * @param smiles
     */
    public void setSmiles(java.lang.String smiles) {
        this.smiles = smiles;
    }

    private java.lang.Object __equalsCalc = null;
    public synchronized boolean equals(java.lang.Object obj) {
        if (!(obj instanceof Prediction)) return false;
        Prediction other = (Prediction) obj;
        if (obj == null) return false;
        if (this == obj) return true;
        if (__equalsCalc != null) {
            return (__equalsCalc == obj);
        }
        __equalsCalc = obj;
        boolean _equals;
        _equals = true && 
            ((this.depictionID==null && other.getDepictionID()==null) || 
             (this.depictionID!=null &&
              this.depictionID.equals(other.getDepictionID()))) &&
            ((this.error==null && other.getError()==null) || 
             (this.error!=null &&
              this.error.equals(other.getError()))) &&
            ((this.inChIKey==null && other.getInChIKey()==null) || 
             (this.inChIKey!=null &&
              this.inChIKey.equals(other.getInChIKey()))) &&
            ((this.moleculeID==null && other.getMoleculeID()==null) || 
             (this.moleculeID!=null &&
              this.moleculeID.equals(other.getMoleculeID()))) &&
            ((this.predictions==null && other.getPredictions()==null) || 
             (this.predictions!=null &&
              java.util.Arrays.equals(this.predictions, other.getPredictions()))) &&
            ((this.smiles==null && other.getSmiles()==null) || 
             (this.smiles!=null &&
              this.smiles.equals(other.getSmiles())));
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
        if (getDepictionID() != null) {
            _hashCode += getDepictionID().hashCode();
        }
        if (getError() != null) {
            _hashCode += getError().hashCode();
        }
        if (getInChIKey() != null) {
            _hashCode += getInChIKey().hashCode();
        }
        if (getMoleculeID() != null) {
            _hashCode += getMoleculeID().hashCode();
        }
        if (getPredictions() != null) {
            for (int i=0;
                 i<java.lang.reflect.Array.getLength(getPredictions());
                 i++) {
                java.lang.Object obj = java.lang.reflect.Array.get(getPredictions(), i);
                if (obj != null &&
                    !obj.getClass().isArray()) {
                    _hashCode += obj.hashCode();
                }
            }
        }
        if (getSmiles() != null) {
            _hashCode += getSmiles().hashCode();
        }
        __hashCodeCalc = false;
        return _hashCode;
    }

    // Type metadata
    private static org.apache.axis.description.TypeDesc typeDesc =
        new org.apache.axis.description.TypeDesc(Prediction.class, true);

    static {
        typeDesc.setXmlType(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "Prediction"));
        org.apache.axis.description.ElementDesc elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("depictionID");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "depictionID"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "long"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("error");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "error"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("inChIKey");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "inChIKey"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("moleculeID");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "moleculeID"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "int"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("predictions");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "predictions"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "PropertyPrediction"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        elemField.setMaxOccursUnbounded(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("smiles");
        elemField.setXmlName(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "smiles"));
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
