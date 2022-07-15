/**
 * IdentitySearchOptions.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

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

package qspr.util.ncbi.pubchem;

public class IdentitySearchOptions  implements java.io.Serializable {
    private qspr.util.ncbi.pubchem.IdentityType eIdentity;

    private java.lang.String toWebEnv;

    public IdentitySearchOptions() {
    }

    public IdentitySearchOptions(
           qspr.util.ncbi.pubchem.IdentityType eIdentity,
           java.lang.String toWebEnv) {
           this.eIdentity = eIdentity;
           this.toWebEnv = toWebEnv;
    }


    /**
     * Gets the eIdentity value for this IdentitySearchOptions.
     * 
     * @return eIdentity
     */
    public qspr.util.ncbi.pubchem.IdentityType getEIdentity() {
        return eIdentity;
    }


    /**
     * Sets the eIdentity value for this IdentitySearchOptions.
     * 
     * @param eIdentity
     */
    public void setEIdentity(qspr.util.ncbi.pubchem.IdentityType eIdentity) {
        this.eIdentity = eIdentity;
    }


    /**
     * Gets the toWebEnv value for this IdentitySearchOptions.
     * 
     * @return toWebEnv
     */
    public java.lang.String getToWebEnv() {
        return toWebEnv;
    }


    /**
     * Sets the toWebEnv value for this IdentitySearchOptions.
     * 
     * @param toWebEnv
     */
    public void setToWebEnv(java.lang.String toWebEnv) {
        this.toWebEnv = toWebEnv;
    }

    private java.lang.Object __equalsCalc = null;
    public synchronized boolean equals(java.lang.Object obj) {
        if (!(obj instanceof IdentitySearchOptions)) return false;
        IdentitySearchOptions other = (IdentitySearchOptions) obj;
        if (obj == null) return false;
        if (this == obj) return true;
        if (__equalsCalc != null) {
            return (__equalsCalc == obj);
        }
        __equalsCalc = obj;
        boolean _equals;
        _equals = true && 
            ((this.eIdentity==null && other.getEIdentity()==null) || 
             (this.eIdentity!=null &&
              this.eIdentity.equals(other.getEIdentity()))) &&
            ((this.toWebEnv==null && other.getToWebEnv()==null) || 
             (this.toWebEnv!=null &&
              this.toWebEnv.equals(other.getToWebEnv())));
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
        if (getEIdentity() != null) {
            _hashCode += getEIdentity().hashCode();
        }
        if (getToWebEnv() != null) {
            _hashCode += getToWebEnv().hashCode();
        }
        __hashCodeCalc = false;
        return _hashCode;
    }

    // Type metadata
    private static org.apache.axis.description.TypeDesc typeDesc =
        new org.apache.axis.description.TypeDesc(IdentitySearchOptions.class, true);

    static {
        typeDesc.setXmlType(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/", "IdentitySearchOptions"));
        org.apache.axis.description.ElementDesc elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("EIdentity");
        elemField.setXmlName(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/", "eIdentity"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/", "IdentityType"));
        elemField.setNillable(false);
        typeDesc.addFieldDesc(elemField);
        elemField = new org.apache.axis.description.ElementDesc();
        elemField.setFieldName("toWebEnv");
        elemField.setXmlName(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/", "ToWebEnv"));
        elemField.setXmlType(new javax.xml.namespace.QName("http://www.w3.org/2001/XMLSchema", "string"));
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
