package qspr.metaserver.cs;

import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsOSMORDREDConfiguration;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class OSMORDREDDescriptorsServer extends DescriptorsAbstractExecutableServer{

	@Override
	protected DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration configuration, int start, int batchSize)
			throws Exception {
		DescriptorsOSMORDREDConfiguration conf = (DescriptorsOSMORDREDConfiguration) configuration;

		System.out.print("running "+ conf);
		saveMolecules(dtMolecules, datain, QSPRConstants.SDFNOAROM_WITHH, start, batchSize);
		String[] commandsLin = {OSType.isMac()?"python" : "/opt/conda/envs/osmordred/bin/python3", getExeFile()};
		runPython(commandsLin, dataout, CONDA.OSMORDRED, batchSize>10?batchSize:10);
		return readStandardOutputResults(getResults(),dataout,false);
	}

	@Override
	int getBatchSize(){
		return 1000;
	}


	public OSMORDREDDescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.OSMORDRED;
		startPosition = 1;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		DELIMITER =",";
	}

}
