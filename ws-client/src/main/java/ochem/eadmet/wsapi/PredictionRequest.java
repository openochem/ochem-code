/**
 * PredictionRequest.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public class PredictionRequest  implements java.io.Serializable {
    private java.lang.Long basketId;

    private java.lang.Long modelId;

    private java.lang.String predictionScenario;

    private java.lang.String[] sdfs;

    private java.lang.String sessionGUID;

    public PredictionRequest() {
    }

    public PredictionRequest(
           java.lang.Long basketId,
           java.lang.Long modelId,
           java.lang.String predictionScenario,
           java.lang.String[] sdfs,
           java.lang.String sessionGUID) {
           this.basketId = basketId;
           this.modelId = modelId;
           this.predictionScenario = predictionScenario;
           this.sdfs = sdfs;
           this.sessionGUID = sessionGUID;
    }


    /**
     * Gets the basketId value for this PredictionRequest.
     * 
     * @return basketId
     */
    public java.lang.Long getBasketId() {
        return basketId;
    }


    /**
     * Sets the basketId value for this PredictionRequest.
     * 
     * @param basketId
     */
    public void setBasketId(java.lang.Long basketId) {
        this.basketId = basketId;
    }


    /**
     * Gets the modelId value for this PredictionRequest.
     * 
     * @return modelId
     */
    public java.lang.Long getModelId() {
        return modelId;
    }


    /**
     * Sets the modelId value for this PredictionRequest.
     * 
     * @param modelId
     */
    public void setModelId(java.lang.Long modelId) {
        this.modelId = modelId;
    }


    /**
     * Gets the predictionScenario value for this PredictionRequest.
     * 
     * @return predictionScenario
     */
    public java.lang.String getPredictionScenario() {
        return predictionScenario;
    }


    /**
     * Sets the predictionScenario value for this PredictionRequest.
     * 
     * @param predictionScenario
     */
    public void setPredictionScenario(java.lang.String predictionScenario) {
        this.predictionScenario = predictionScenario;
    }


    /**
     * Gets the sdfs value for this PredictionRequest.
     * 
     * @return sdfs
     */
    public java.lang.String[] getSdfs() {
        return sdfs;
    }


    /**
     * Sets the sdfs value for this PredictionRequest.
     * 
     * @param sdfs
     */
    public void setSdfs(java.lang.String[] sdfs) {
        this.sdfs = sdfs;
    }

    public java.lang.String getSdfs(int i) {
        return this.sdfs[i];
    }

    public void setSdfs(int i, java.lang.String _value) {
        this.sdfs[i] = _value;
    }


    /**
     * Gets the sessionGUID value for this PredictionRequest.
     * 
     * @return sessionGUID
     */
    public java.lang.String getSessionGUID() {
        return sessionGUID;
    }


    /**
     * Sets the sessionGUID value for this PredictionRequest.
     * 
     * @param sessionGUID
     */
    public void setSessionGUID(java.lang.String sessionGUID) {
        this.sessionGUID = sessionGUID;
    }

    private java.lang.Object __equalsCalc = null;
    public synchronized boolean equals(java.lang.Object obj) {
        if (!(obj instanceof PredictionRequest)) return false;
        PredictionRequest other = (PredictionRequest) obj;
        if (obj == null) return false;
        if (this == obj) return true;
        if (__equalsCalc != null) {
            return (__equalsCalc == obj);
        }
        __equalsCalc = obj;
        boolean _equals;
        _equals = true && 
            ((this.basketId==null && other.getBasketId()==null) || 
             (this.basketId!=null &&
              this.basketId.equals(other.getBasketId()))) &&
            ((this.modelId==null && other.getModelId()==null) || 
             (this.modelId!=null &&
              this.modelId.equals(other.getModelId()))) &&
            ((this.predictionScenario==null && other.getPredictionScenario()==null) || 
             (this.predictionScenario!=null &&
              this.predictionScenario.equals(other.getPredictionScenario()))) &&
            ((this.sdfs==null && other.getSdfs()==null) || 
             (this.sdfs!=null &&
              java.util.Arrays.equals(this.sdfs, other.getSdfs()))) &&
            ((this.sessionGUID==null && other.getSessionGUID()==null) || 
             (this.sessionGUID!=null &&
              this.sessionGUID.equals(other.getSessionGUID())));
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
        if (getBasketId() != null) {
            _hashCode += getBasketId().hashCode();
        }
        if (getModelId() != null) {
            _hashCode += getModelId().hashCode();
        }
        if (getPredictionScenario() != null) {
            _hashCode += getPredictionScenario().hashCode();
        }
        if (getSdfs() != null) {
            for (int i=0;
                 i<java.lang.reflect.Array.getLength(getSdfs());
                 i++) {
                java.lang.Object obj = java.lang.reflect.Array.get(getSdfs(), i);
                if (obj != null &&
                    !obj.getClass().isArray()) {
                    _hashCode += obj.hashCode();
                }
            }
        }
        if (getSessionGUID() != null) {
            _hashCode += getSessionGUID().hashCode();
        }
        __hashCodeCalc = false;
        return _hashCode;
    }

    // Type metadata
    private static org.apache.axis.description.TypeDesc typeDesc =
        new org.apache.axis.description.TypeDesc(PredictionRequest.class, true);

    static {
        typeDesc.setXmlType(new javax.xml.namespace.QName("http://services.qspr/xsd", "PredictionRequest"));
        org.apache.axis.description.ElementDesc elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("basketId");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "basketId"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "long"));
        elemField.setMinOccurs(0);
        elemField.setNillable(false);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("modelId");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "modelId"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "long"));
        elemField.setMinOccurs(0);
        elemField.setNillable(false);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("predictionScenario");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "predictionScenario"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("sdfs");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "sdfs"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        elemField.setMaxOccursUnbounded(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("sessionGUID");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "sessionGUID"));
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
