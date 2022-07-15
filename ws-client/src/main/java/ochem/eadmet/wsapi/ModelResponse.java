/**
 * ModelResponse.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public class ModelResponse  implements java.io.Serializable {
    private java.lang.String detailedStatus;

    private java.lang.Long metaserverTaskId;

    private java.lang.String modelDescriptionUrl;

    private ochem.eadmet.wsapi.Prediction[] predictions;

    private java.lang.String status;

    private java.lang.Long taskId;

    public ModelResponse() {
    }

    public ModelResponse(
           java.lang.String detailedStatus,
           java.lang.Long metaserverTaskId,
           java.lang.String modelDescriptionUrl,
           ochem.eadmet.wsapi.Prediction[] predictions,
           java.lang.String status,
           java.lang.Long taskId) {
           this.detailedStatus = detailedStatus;
           this.metaserverTaskId = metaserverTaskId;
           this.modelDescriptionUrl = modelDescriptionUrl;
           this.predictions = predictions;
           this.status = status;
           this.taskId = taskId;
    }


    /**
     * Gets the detailedStatus value for this ModelResponse.
     * 
     * @return detailedStatus
     */
    public java.lang.String getDetailedStatus() {
        return detailedStatus;
    }


    /**
     * Sets the detailedStatus value for this ModelResponse.
     * 
     * @param detailedStatus
     */
    public void setDetailedStatus(java.lang.String detailedStatus) {
        this.detailedStatus = detailedStatus;
    }


    /**
     * Gets the metaserverTaskId value for this ModelResponse.
     * 
     * @return metaserverTaskId
     */
    public java.lang.Long getMetaserverTaskId() {
        return metaserverTaskId;
    }


    /**
     * Sets the metaserverTaskId value for this ModelResponse.
     * 
     * @param metaserverTaskId
     */
    public void setMetaserverTaskId(java.lang.Long metaserverTaskId) {
        this.metaserverTaskId = metaserverTaskId;
    }


    /**
     * Gets the modelDescriptionUrl value for this ModelResponse.
     * 
     * @return modelDescriptionUrl
     */
    public java.lang.String getModelDescriptionUrl() {
        return modelDescriptionUrl;
    }


    /**
     * Sets the modelDescriptionUrl value for this ModelResponse.
     * 
     * @param modelDescriptionUrl
     */
    public void setModelDescriptionUrl(java.lang.String modelDescriptionUrl) {
        this.modelDescriptionUrl = modelDescriptionUrl;
    }


    /**
     * Gets the predictions value for this ModelResponse.
     * 
     * @return predictions
     */
    public ochem.eadmet.wsapi.Prediction[] getPredictions() {
        return predictions;
    }


    /**
     * Sets the predictions value for this ModelResponse.
     * 
     * @param predictions
     */
    public void setPredictions(ochem.eadmet.wsapi.Prediction[] predictions) {
        this.predictions = predictions;
    }

    public ochem.eadmet.wsapi.Prediction getPredictions(int i) {
        return this.predictions[i];
    }

    public void setPredictions(int i, ochem.eadmet.wsapi.Prediction _value) {
        this.predictions[i] = _value;
    }


    /**
     * Gets the status value for this ModelResponse.
     * 
     * @return status
     */
    public java.lang.String getStatus() {
        return status;
    }


    /**
     * Sets the status value for this ModelResponse.
     * 
     * @param status
     */
    public void setStatus(java.lang.String status) {
        this.status = status;
    }


    /**
     * Gets the taskId value for this ModelResponse.
     * 
     * @return taskId
     */
    public java.lang.Long getTaskId() {
        return taskId;
    }


    /**
     * Sets the taskId value for this ModelResponse.
     * 
     * @param taskId
     */
    public void setTaskId(java.lang.Long taskId) {
        this.taskId = taskId;
    }

    private java.lang.Object __equalsCalc = null;
    public synchronized boolean equals(java.lang.Object obj) {
        if (!(obj instanceof ModelResponse)) return false;
        ModelResponse other = (ModelResponse) obj;
        if (obj == null) return false;
        if (this == obj) return true;
        if (__equalsCalc != null) {
            return (__equalsCalc == obj);
        }
        __equalsCalc = obj;
        boolean _equals;
        _equals = true && 
            ((this.detailedStatus==null && other.getDetailedStatus()==null) || 
             (this.detailedStatus!=null &&
              this.detailedStatus.equals(other.getDetailedStatus()))) &&
            ((this.metaserverTaskId==null && other.getMetaserverTaskId()==null) || 
             (this.metaserverTaskId!=null &&
              this.metaserverTaskId.equals(other.getMetaserverTaskId()))) &&
            ((this.modelDescriptionUrl==null && other.getModelDescriptionUrl()==null) || 
             (this.modelDescriptionUrl!=null &&
              this.modelDescriptionUrl.equals(other.getModelDescriptionUrl()))) &&
            ((this.predictions==null && other.getPredictions()==null) || 
             (this.predictions!=null &&
              java.util.Arrays.equals(this.predictions, other.getPredictions()))) &&
            ((this.status==null && other.getStatus()==null) || 
             (this.status!=null &&
              this.status.equals(other.getStatus()))) &&
            ((this.taskId==null && other.getTaskId()==null) || 
             (this.taskId!=null &&
              this.taskId.equals(other.getTaskId())));
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
        if (getDetailedStatus() != null) {
            _hashCode += getDetailedStatus().hashCode();
        }
        if (getMetaserverTaskId() != null) {
            _hashCode += getMetaserverTaskId().hashCode();
        }
        if (getModelDescriptionUrl() != null) {
            _hashCode += getModelDescriptionUrl().hashCode();
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
        if (getStatus() != null) {
            _hashCode += getStatus().hashCode();
        }
        if (getTaskId() != null) {
            _hashCode += getTaskId().hashCode();
        }
        __hashCodeCalc = false;
        return _hashCode;
    }

    // Type metadata
    private static org.apache.axis.description.TypeDesc typeDesc =
        new org.apache.axis.description.TypeDesc(ModelResponse.class, true);

    static {
        typeDesc.setXmlType(new javax.xml.namespace.QName("http://services.qspr/xsd", "ModelResponse"));
        org.apache.axis.description.ElementDesc elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("detailedStatus");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "detailedStatus"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("metaserverTaskId");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "metaserverTaskId"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "long"));
        elemField.setMinOccurs(0);
        elemField.setNillable(false);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("modelDescriptionUrl");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "modelDescriptionUrl"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("predictions");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "predictions"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://applier.modelling.qspr/xsd", "Prediction"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        elemField.setMaxOccursUnbounded(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("status");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "status"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
        elemField.setMinOccurs(0);
        elemField.setNillable(true);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("taskId");
        elemField.setXmlName(new javax.xml.namespace.QName("http://services.qspr/xsd", "taskId"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "long"));
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
