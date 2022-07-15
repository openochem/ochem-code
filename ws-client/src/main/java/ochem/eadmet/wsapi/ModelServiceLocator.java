/**
 * ModelServiceLocator.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public class ModelServiceLocator extends org.apache.axis.client.Service implements ochem.eadmet.wsapi.ModelService {

    public ModelServiceLocator() {
    }


    public ModelServiceLocator(org.apache.axis.EngineConfiguration config) {
        super(config);
    }

    public ModelServiceLocator(java.lang.String wsdlLoc, javax.xml.namespace.QName sName) throws javax.xml.rpc.ServiceException {
        super(wsdlLoc, sName);
    }

    // Use to get a proxy class for ModelServiceHttpSoap11Endpoint
    private java.lang.String ModelServiceHttpSoap11Endpoint_address = "http://localhost:8080/services/ModelService.ModelServiceHttpSoap11Endpoint/";

    public java.lang.String getModelServiceHttpSoap11EndpointAddress() {
        return ModelServiceHttpSoap11Endpoint_address;
    }

    // The WSDD service name defaults to the port name.
    private java.lang.String ModelServiceHttpSoap11EndpointWSDDServiceName = "ModelServiceHttpSoap11Endpoint";

    public java.lang.String getModelServiceHttpSoap11EndpointWSDDServiceName() {
        return ModelServiceHttpSoap11EndpointWSDDServiceName;
    }

    public void setModelServiceHttpSoap11EndpointWSDDServiceName(java.lang.String name) {
        ModelServiceHttpSoap11EndpointWSDDServiceName = name;
    }

    public ochem.eadmet.wsapi.ModelServicePortType getModelServiceHttpSoap11Endpoint() throws javax.xml.rpc.ServiceException {
       java.net.URL endpoint;
        try {
            endpoint = new java.net.URL(ModelServiceHttpSoap11Endpoint_address);
        }
        catch (java.net.MalformedURLException e) {
            throw new javax.xml.rpc.ServiceException(e);
        }
        return getModelServiceHttpSoap11Endpoint(endpoint);
    }

    public ochem.eadmet.wsapi.ModelServicePortType getModelServiceHttpSoap11Endpoint(java.net.URL portAddress) throws javax.xml.rpc.ServiceException {
        try {
            ochem.eadmet.wsapi.ModelServiceSoap11BindingStub _stub = new ochem.eadmet.wsapi.ModelServiceSoap11BindingStub(portAddress, this);
            _stub.setPortName(getModelServiceHttpSoap11EndpointWSDDServiceName());
            return _stub;
        }
        catch (org.apache.axis.AxisFault e) {
            return null;
        }
    }

    public void setModelServiceHttpSoap11EndpointEndpointAddress(java.lang.String address) {
        ModelServiceHttpSoap11Endpoint_address = address;
    }


    // Use to get a proxy class for ModelServiceHttpSoap12Endpoint
    private java.lang.String ModelServiceHttpSoap12Endpoint_address = "http://localhost:8080/services/ModelService.ModelServiceHttpSoap12Endpoint/";

    public java.lang.String getModelServiceHttpSoap12EndpointAddress() {
        return ModelServiceHttpSoap12Endpoint_address;
    }

    // The WSDD service name defaults to the port name.
    private java.lang.String ModelServiceHttpSoap12EndpointWSDDServiceName = "ModelServiceHttpSoap12Endpoint";

    public java.lang.String getModelServiceHttpSoap12EndpointWSDDServiceName() {
        return ModelServiceHttpSoap12EndpointWSDDServiceName;
    }

    public void setModelServiceHttpSoap12EndpointWSDDServiceName(java.lang.String name) {
        ModelServiceHttpSoap12EndpointWSDDServiceName = name;
    }

    public ochem.eadmet.wsapi.ModelServicePortType getModelServiceHttpSoap12Endpoint() throws javax.xml.rpc.ServiceException {
       java.net.URL endpoint;
        try {
            endpoint = new java.net.URL(ModelServiceHttpSoap12Endpoint_address);
        }
        catch (java.net.MalformedURLException e) {
            throw new javax.xml.rpc.ServiceException(e);
        }
        return getModelServiceHttpSoap12Endpoint(endpoint);
    }

    public ochem.eadmet.wsapi.ModelServicePortType getModelServiceHttpSoap12Endpoint(java.net.URL portAddress) throws javax.xml.rpc.ServiceException {
        try {
            ochem.eadmet.wsapi.ModelServiceSoap12BindingStub _stub = new ochem.eadmet.wsapi.ModelServiceSoap12BindingStub(portAddress, this);
            _stub.setPortName(getModelServiceHttpSoap12EndpointWSDDServiceName());
            return _stub;
        }
        catch (org.apache.axis.AxisFault e) {
            return null;
        }
    }

    public void setModelServiceHttpSoap12EndpointEndpointAddress(java.lang.String address) {
        ModelServiceHttpSoap12Endpoint_address = address;
    }

    /**
     * For the given interface, get the stub implementation.
     * If this service has no port for the given interface,
     * then ServiceException is thrown.
     * This service has multiple ports for a given interface;
     * the proxy implementation returned may be indeterminate.
     */
    public java.rmi.Remote getPort(Class serviceEndpointInterface) throws javax.xml.rpc.ServiceException {
        try {
            if (ochem.eadmet.wsapi.ModelServicePortType.class.isAssignableFrom(serviceEndpointInterface)) {
                ochem.eadmet.wsapi.ModelServiceSoap11BindingStub _stub = new ochem.eadmet.wsapi.ModelServiceSoap11BindingStub(new java.net.URL(ModelServiceHttpSoap11Endpoint_address), this);
                _stub.setPortName(getModelServiceHttpSoap11EndpointWSDDServiceName());
                return _stub;
            }
            if (ochem.eadmet.wsapi.ModelServicePortType.class.isAssignableFrom(serviceEndpointInterface)) {
                ochem.eadmet.wsapi.ModelServiceSoap12BindingStub _stub = new ochem.eadmet.wsapi.ModelServiceSoap12BindingStub(new java.net.URL(ModelServiceHttpSoap12Endpoint_address), this);
                _stub.setPortName(getModelServiceHttpSoap12EndpointWSDDServiceName());
                return _stub;
            }
        }
        catch (java.lang.Throwable t) {
            throw new javax.xml.rpc.ServiceException(t);
        }
        throw new javax.xml.rpc.ServiceException("There is no stub implementation for the interface:  " + (serviceEndpointInterface == null ? "null" : serviceEndpointInterface.getName()));
    }

    /**
     * For the given interface, get the stub implementation.
     * If this service has no port for the given interface,
     * then ServiceException is thrown.
     */
    public java.rmi.Remote getPort(javax.xml.namespace.QName portName, Class serviceEndpointInterface) throws javax.xml.rpc.ServiceException {
        if (portName == null) {
            return getPort(serviceEndpointInterface);
        }
        java.lang.String inputPortName = portName.getLocalPart();
        if ("ModelServiceHttpSoap11Endpoint".equals(inputPortName)) {
            return getModelServiceHttpSoap11Endpoint();
        }
        else if ("ModelServiceHttpSoap12Endpoint".equals(inputPortName)) {
            return getModelServiceHttpSoap12Endpoint();
        }
        else  {
            java.rmi.Remote _stub = getPort(serviceEndpointInterface);
            ((org.apache.axis.client.Stub) _stub).setPortName(portName);
            return _stub;
        }
    }

    public javax.xml.namespace.QName getServiceName() {
        return new javax.xml.namespace.QName("http://services.qspr", "ModelService");
    }

    private java.util.HashSet ports = null;

    public java.util.Iterator getPorts() {
        if (ports == null) {
            ports = new java.util.HashSet();
            ports.add(new javax.xml.namespace.QName("http://services.qspr", "ModelServiceHttpSoap11Endpoint"));
            ports.add(new javax.xml.namespace.QName("http://services.qspr", "ModelServiceHttpSoap12Endpoint"));
        }
        return ports.iterator();
    }

    /**
    * Set the endpoint address for the specified port name.
    */
    public void setEndpointAddress(java.lang.String portName, java.lang.String address) throws javax.xml.rpc.ServiceException {
        
if ("ModelServiceHttpSoap11Endpoint".equals(portName)) {
            setModelServiceHttpSoap11EndpointEndpointAddress(address);
        }
        else 
if ("ModelServiceHttpSoap12Endpoint".equals(portName)) {
            setModelServiceHttpSoap12EndpointEndpointAddress(address);
        }
        else 
{ // Unknown Port Name
            throw new javax.xml.rpc.ServiceException(" Cannot set Endpoint Address for Unknown Port" + portName);
        }
    }

    /**
    * Set the endpoint address for the specified port name.
    */
    public void setEndpointAddress(javax.xml.namespace.QName portName, java.lang.String address) throws javax.xml.rpc.ServiceException {
        setEndpointAddress(portName.getLocalPart(), address);
    }

}
