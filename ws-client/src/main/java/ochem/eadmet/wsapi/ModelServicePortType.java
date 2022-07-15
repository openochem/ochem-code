/**
 * ModelServicePortType.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis 1.4 Apr 22, 2006 (06:55:48 PDT) WSDL2Java emitter.
 */

package ochem.eadmet.wsapi;

public interface ModelServicePortType extends java.rmi.Remote {
    public ochem.eadmet.wsapi.ModelResponse postModelWithSession(java.lang.String sessionGUID, java.lang.Long modelId, java.lang.String[] sdfs) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.PropertyPrediction[] getModelPredictions(java.lang.String sessionGUID, java.lang.Long publicModelId) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.ModelResponse fetchModel(java.lang.String sessionGUID, java.lang.Long taskId) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.ModelResponse postModelRequest(ochem.eadmet.wsapi.PredictionRequest request) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.ModelResponse postDataReferenceRequest(ochem.eadmet.wsapi.DataReferenceRequest request) throws java.rmi.RemoteException;
    public java.lang.String login(java.lang.String username, java.lang.String password) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.ModelResponse applyModelSingleSDF(java.lang.String sessionGUID, java.lang.Long modelId, java.lang.String sdf) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.ModelResponse deleteTask(java.lang.String sessionGUID, java.lang.Long taskId) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.ModelResponse getPredictions(java.lang.String sessionGUID, java.lang.Long taskId, java.lang.String[] measurementUnit) throws java.rmi.RemoteException;
    public ochem.eadmet.wsapi.ModelSummary[] getModelSummary(java.lang.String sessionGUID, java.lang.Long publicModelId) throws java.rmi.RemoteException;
}
