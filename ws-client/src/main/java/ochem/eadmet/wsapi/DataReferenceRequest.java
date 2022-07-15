/**
 * DataReferenceRequest.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public class DataReferenceRequest  implements java.io.Serializable {
    private java.lang.String dataReference;

    private java.lang.Integer datasize;

    private java.lang.Long modelId;

    private java.lang.String predictionScenario;

    private java.lang.String sessionGUID;

    public DataReferenceRequest() {
    }

    public DataReferenceRequest(
           java.lang.String dataReference,
           java.lang.Integer datasize,
           java.lang.Long modelId,
           java.lang.String predictionScenario,
           java.lang.String sessionGUID) {
           this.dataReference = dataReference;
           this.datasize = datasize;
           this.modelId = modelId;
           this.predictionScenario = predictionScenario;
           this.sessionGUID = sessionGUID;
    }


    /**
     * Gets the dataReference value for this DataReferenceRequest.
     * 
     * @return dataReference
     */
    public java.lang.String getDataReference() {
        return dataReference;
    }


    /**
     * Sets the dataReference value for this DataReferenceRequest.
     * 
     * @param dataReference
     */
    public void setDataReference(java.lang.String dataReference) {
        this.dataReference = dataReference;
    }


    /**
     * Gets the datasize value for this DataReferenceRequest.
     * 
     * @return datasize
     */
    public java.lang.Integer getDatasize() {
        return datasize;
    }


    /**
     * Sets the datasize value for this DataReferenceRequest.
     * 
     * @param datasize
     */
    public void setDatasize(java.lang.Integer datasize) {
        this.datasize = datasize;
    }


    /**
     * Gets the modelId value for this DataReferenceRequest.
     * 
     * @return modelId
     */
    public java.lang.Long getModelId() {
        return modelId;
    }


    /**
     * Sets the modelId value for this DataReferenceRequest.
     * 
     * @param modelId
     */
    public void setModelId(java.lang.Long modelId) {
        this.modelId = modelId;
    }


    /**
     * Gets the predictionScenario value for this DataReferenceRequest.
     * 
     * @return predictionScenario
     */
    public java.lang.String getPredictionScenario() {
        return predictionScenario;
    }


    /**
     * Sets the predictionScenario value for this DataReferenceRequest.
     * 
     * @param predictionScenario
     */
    public void setPredictionScenario(java.lang.String predictionScenario) {
        this.predictionScenario = predictionScenario;
    }


    /**
     * Gets the sessionGUID value for this DataReferenceRequest.
     * 
     * @return sessionGUID
     */
    public java.lang.String getSessionGUID() {
        return sessionGUID;
    }


    /**
     * Sets the sessionGUID value for this DataReferenceRequest.
     * 
     * @param sessionGUID
     */
    public void setSessionGUID(java.lang.String sessionGUID) {
        this.sessionGUID = sessionGUID;
    }

    private java.lang.Object __equalsCalc = null;
    public synchronized boolean equals(java.lang.Object obj) {
        if (!(obj instanceof DataReferenceRequest)) return false;
        DataReferenceRequest other = (DataReferenceRequest) obj;
        if (obj == null) return false;
        if (this == obj) return true;
        if (__equalsCalc != null) {
            return (__equalsCalc == obj);
        }
        __equalsCalc = obj;
        boolean _equals;
        _equals = true && 
            ((this.dataReference==null && other.getDataReference()==null) || 
             (this.dataReference!=null &&
              this.dataReference.equals(other.getDataReference()))) &&
            ((this.datasize==null && other.getDatasize()==null) || 
             (this.datasize!=null &&
              this.datasize.equals(other.getDatasize()))) &&
            ((this.modelId==null && other.getModelId()==null) || 
             (this.modelId!=null &&
              this.modelId.equals(other.getModelId()))) &&
            ((this.predictionScenario==null && other.getPredictionScenario()==null) || 
             (this.predictionScenario!=null &&
              this.predictionScenario.equals(other.getPredictionScenario()))) &&
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
        if (getDataReference() != null) {
            _hashCode += getDataReference().hashCode();
        }
        if (getDatasize() != null) {
            _hashCode += getDatasize().hashCode();
        }
        if (getModelId() != null) {
            _hashCode += getModelId().hashCode();
        }
        if (getPredictionScenario() != null) {
            _hashCode += getPredictionScenario().hashCode();
        }
        if (getSessionGUID() != null) {
            _hashCode += getSessionGUID().hashCode();
        }
        __hashCodeCalc = false;
        return _hashCode;
    }

    // Type metadata
    private static org.apache.axis.description.TypeDesc typeDesc =
        new org.apache.axis.description.TypeDesc(DataReferenceRequest.class, true);

    static {
        typeDesc.setXmlType(new javax.xml.namespace.QName("http://services.qspr/xsd", "DataReferenceRequest"));
        org.apache.axis.description.ElementDesc elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("dataReference");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "dataReference"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("datasize");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "datasize"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "int"));
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
