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

package qspr.workflow.utils;

import qspr.metaserver.protocol.Task;

public class QSPRConstants
{

	/** MOLECULES */
	public static final String SDF_COLUMN = "SDF"; // distinguishing column in the data table and format
	public static final String SDF3D_COLUMN = SDF_COLUMN; /// distinguish path for 3D

	//  conversion formats as defined in ChemAxon
	public static final String PDB = "PDB";
	public static final String SDF = "SDF";
	public static final String SDFH = "SDF:H"; //add all  hydrogens
	public static final String SDFNOAROM_NOH = "SDF:-a-H"; //remove all  hydrogens, KEKULE representation
	public static final String SDFNOAROM_WITHH = "SDF:-a+H"; // added implicit hydrogens
	public static final String SDFAROM_BASIC_WITHH = "SDF:a_bas+H"; // aromatized basic with hydrogens basic
	public static final String SDFAROM_BASIC_NOH = "SDF:a_bas-H"; // aromatized basic without hydrogens basic
	public static final String SDFAROM_GENERAL_WITHH = "SDF:a_loose+H"; // aromatized general with hydrogens basic

	public static final String SMILES_FORMAT = "SMILES";
	public static final String SMILES_ATTACHMENT = "SMILES";
	public static final String SMILESH = "SMILES:-H";
	public static final String SMILESNOSTEREO = "SMILES:0-H";
	public static final String SMILESNOAROM = "SMILES:-a-H";
	public static final String SMILESNOAROMNOSTEREO = "SMILES:0-a-H";
	public static final String SMILESNOAROM_WITHH = "SMILES:-a+H"; // added implicit hydrogens
	public static final String SMILESUniqueNoHAromatic = "SMILES:u-H"; // Unique aromatic SMILES with suppressed implicit hydrogens
	public static final String SMILESSRC = "SRC"; 

	public static final String MOL2 = "MOL2";
	public static final String PEPTIDES = "peptide:1";
	public static final String MOPIN = "mopin";
	public static final String CXSMARTS = "cxsmarts";
	public static final String INCHIKEYS = "InChIKey";
	public static final String ASIS = "RAW"; // not changed - as is

	/** MODEL */
	// limitations
	public static final Integer MODEL_REPOST_SIZE = 32768; // to fit in 5*128*3 = 2 GB servers

	// fields
	public static final String DM = "DM"; 
	public static final String DM_PREFIX = "DM:";
	public static final String CLASS = "CLASS"; 
	public static final String VALUE = "VALUE"; 
	public static final String CLASSLAG = ":CLASS-LAG"; 
	public static final String PROBSTD = ":PROB-STD"; 
	public static final String STDEV = ":ASNN-STDEV";
	public static final String STDEVAUG = ":AUG-STDEV";
	public static final String CORREL = ":ASNN-CORREL"; // why DM names start with ":"? its logical
	public static final String BAGGING_STD = "BAGGING-STD";
	public static final String CONSENSUS_STD = "CONSENSUS-STD";
	public static final String WHOLE = "whole";
	public static final String EXCLUDED = "excluded";
	public static final String ADDED_MOLECULES = "ADDED_MOLECULES";
	public static final String MIXTURE_CONDITION = "MIXTURE";
	public static final String MIXTURE_ATTACHMENT = MIXTURE_CONDITION; // different names just for logic
	public static final String SOLVENT_CONDITION = "solvent";
	public static final String SOLVENT_ATTACHMENT = "solvent_attachment"; // different names just for logic
	public static final String REACTION_CONDITION = "REACTION";

	// Attachment fields
	public static final String IS_CONDITION_COLUMN = "isCondition"; //  maybe not used?!
	public static final String IS_QUALITATIVE_COLUMN = "isQualitative";
	public static final String PREDICTION_RESULT_COLUMN = "Result";
	public static final String PREDICATE_ATTACHMENT = "Predicate";
	public static final String RECORD_ID_ATTACHMENT = "EPID";
	public static final String MOLECULE_ID_STEREOCHEM = "MP2"; // provides mapping to molecules with this stereochemisty, Mapping2 id
	public static final String MOLECULE_ID_NO_STEREOCHEM = "ID"; // provides mapping to molecules without stereochemistry, Mapping1 id
	public static final String MOLECULE_ID_OCHEM = "ID_OCHEM"; // original ID in the database, the most full reference
	public static final String PROPERTY_ID_OCHEM = "ID_PROPERTY"; // original ID in the database, the most full reference
	public static final String OPTION_ID_OCHEM = "ID_OPTION"; // original ID in the database, the most full reference
	public static final String EXTERNAL_ID = "ORIGINAL_TAG";
	public static final String VALIDATION = "validation";
	public static final String TRAINING = "training";
	public static final String IMPLICIT_ATTACHMENT = "IMPLICIT";

	//public static final String mixtureConditionMolId = "Mix_C2_MID";
	//public static final String mixtureConditionMolName = "Mix_C2_Name";
	//public static final String mixtureConditionMolarFraction = "Mix_MolFr1"; 

	// atomic descriptor constants
	public static final String IC_INDEX = "IonizableCenterIndex";
	public static final String MOPAC_DERIVED_CHARGES="MOPAC_DERIVED_CHARGES";


	/** OCHEM */
	//  global management parameters
	public static final Integer XEMISTRY_INDEXING_TASK = 1000;
	public static final Integer ECXEMISTRY_INDEXING_TASK = 50000;
	public static final Integer CACHING_SIZE_TASK = 100000;
	public static final int SUPERUSER_RECORDS_LIMIT = 2000000; 
	public static final int VALIDATEDUSER_RECORDS_LIMIT = 100000; 
	public static final int USER_RECORDS_LIMIT = 30000;
	public static final int INHOUSE_XEMISTRY = 1000000000;
	public static final int MIN_PROPERTY_COUNT = 50;

	/** Mongodb */
	//  database names
	public static final String DESCRIPTORSCACHE = "descriptorscache";

	/** TASK */
	// General server constants
	public static final String EXE="EXE";
	public static final String ERROR = "ERROR";

	/** Various constants */
	public static final String MODIFIED_POINTS_ID = "modified-on-the-fly";
	public static final String TASK_TEMPORAL_FAILURE="TemporaryFailure";
	public static final String INDIVIDUAL_PREDICTIONS = "IndividualPredictions"; // predictions of all networks in Bagging; to determine nearets neighbours in prediction spaces
	public static final Integer SEED = 10666; // default seed number
	public static final String OMP_NUM_THREADS = "OMP_NUM_THREADS"; // number of threads for parallel calculations
	public static final Integer FILE_BUFFER = 65536;

	public static final String NO_VALIDATION = "no-validation";

	//statuses
	public static final String ERROR_STATUS = Task.ERROR; // not nice, to be changes for better logic one day
	public static final String WEBSERVICE_DATABASE = "webservice";

	public static final String PDF = "pdf";
	public static final String EXCEL = "xls";
	public static final String EXCEL_NEW = "xlsx";
	public static final String ANONYMOUS = "anonymous";
	public static final String CACHED = "cached";

	// Workflow types
	public static final String Workflow = Task.Workflow;
	public static final String WEBSERVICE_TASK = "Applied using WebService";

	// Model template types
	public static final String UPLOADED_MODEL = "Uploaded model";

	//MOLECULE identificators
	public static final String MOLECULE   = "MOLECULE";
	public static final String MOLECULEID = "MOLECULEID";
	public static final String EXTERNALID = "EXTERNALID";

	// various constants
	public static final String EMPTY_MOL ="1\n2\n3\n  1  0  0  0  0  0            999 V2000\n    0.0000   -0.0000   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\nM  END";
	public static final String EMPTY_MD5 ="98dd04fe45ae60877d0c88efbf61a7fd";	
	public static final String CONDITIONS = "conditions";
	public static final String DEFAULTDATABASEDRIVER = "org.mariadb.jdbc.Driver";
	public static final Long   UNPUBLISHED_JOURNAL = 1l;
	public static final int    MODEL_BONUS = 15;
	public static final int    APPLIER_BONUS = 100;
	public static final int    DESCRIPTORS_BONUS = 100;
	public static final int    UPDATE_RECORDS_COUNT_TIME = 300; // time to update sizes of searches
	public static final String NOSTORED = "No stored descriptors found for this molecule";
	public static final int    MAXIMAL_PUBLIC_ID = 1000000;

	// Publishing
	public static final long     PUBLISHER_ID = 1; // PUBLISHER id of a user which formally publish models, data, stored descriptors
	public static final String   PUBLISHER = "published"; // The login of this special user

	// Methods
	public static final String COMPARE_SCAFFOLDS = "CompareScaffolds";
	public static final String CONSENSUS = "Consensus";
	public static final String BAGGING = "Bagging";
	public static final String CROSSVALIDATION = "CrossValidation";
	public static final String CV = "cv";

	//public static final String DOPTIMAL = "DOptimal";
	public static final String SELECTION = "Selection";
	public static final String ExternalDescriptors = "ExternalDescriptors";
	public static final String LEVERAGE = "Leverage";
	public static final String MolStandartizer = "MolStandartizer";
	public static final String DESCRIPTORS = "Descriptors";
	public static final String DESCRIPTOR = "desc";
	public static final String MMPFrag = "MMPFrag";

	public static final String J48 = "WEKA-J48";
	public static final String RF = "WEKA-RF";

	public static final String KNN = "KNN";
	public static final String KPLS = "KPLS";
	public static final String DEEP = "DEEP";
	public static final String LIBSVM = "LibSVM";
	public static final String FSMLR = "FSMLR";
	public static final String MLRA = "MLRA";
	public static final String CHEMPROP = "ChemProp";
	public static final String CNF = "CNF"; //TensorFlow
	public static final String GNN = "GNN"; //TensorFlow
	public static final String KGCNN = "KGCNN"; //KGCNN TensorFlow servers
	public static final String ATTFP = "ATTFP"; //TensorFlow
	public static final String IMGQSAR = "IMGQSAR"; //TensorFlow
	public static final String TRANSNN = "TRANSNN"; //TensorFlow
	public static final String PYTORCH = "PYTORCH"; 
	public static final String DIMENET = "DIMENET"; 
	public static final String HAMNET = "HAMNET"; 
	public static final String PLS = "PLS";
	public static final String RFR = "RFR";
	public static final String XGBOOST = "XGBOOST";
	public static final String SKLEARN = "SKLEARN";
	public static final String EAGCNG = "EAGCNG";

	// multilearning approaches
	public static final String ANN = "ANN";
	public static final String DNN = "DNN";
	public static final String ASNN = "ASNN";
	public static final String MACAU = "MACAU";
	public static final String LSSVMG = "LSSVMG";
	public static final String DEEPCHEM = "DEEPCHEM";
	public static final String DLCA = "DLCA";

	// structure optimisation servers
	public static final String CORINA = "Corina";
	public static final String BALLOON = "BALLOON";
	public static final String OBABEL = "OBABEL";
	public static final String OBGEN = "OBGEN";	

	public static final String EMPTY_MOLECULE_ERROR = "Empty molecule was provided.";
	public static final String SKIP = "skip";
	public static final String USE = "use";
	public static final String USED = ":used";
	public static final String TEST_USER_PREFIX = Task.TEST_USER_PREFIX;
	public static final Integer MAX_AUGMENTATION = 100;
	public static final Integer NO_REPLICAS = 0;
	public static final int MAX_ATOMS = 600;
	public static final String AUGMENTATIONS = "AUGMENTATIONS";
	public static final long MAXMODELSIZE = 1024 * 1024 * 1024; // 1GB
	public static final String INFOEMAIL = "info@ochem.eu";
	public static final int MAXBAGGING = 128;
	public static final int MINBAGGING = 32;
	public static final int MAXCV = 16;
	public static final int MINCV = 2;
	public static final String FEATURE_NET_PREFIX = "feature_net";
	public static final String UNSAVED = " unsaved";
	public static final String RESTRICT = "Because the OCHEM web site was heavily overloaded due to the missuse of this feature, it is currently accessible only to the validated users.";
	public static final String RDF_SMILES = "RDFSMILE";
	public static final String RDF_PROP = ">  <"+RDF_SMILES+">";
	public static final String UNPUBLISHED_ARTICLE = "Unpublished data of ";
	public static final Integer JAVA_MIN_MEMORY = 2048;
	public static final String IMMEDIATE_FAILURE = "FAILURE:"; // to be also in MetServer
	public static final int SMARTTIMEOUT = 900; // in ms
	//public static final String DYNAMICBOND = "dynbond";
	public static final int MAX_MODELS = 40;
	public static final String RESET = "NNNNN>>NNN.NN"; // to reset cache in reaction upload
	public static final int ERRORS_CACHE = 3;

	// libraries for conversion of molecules - N.B.! to be fixed later -- code is non reliable here
	public static final String TOOLS ="tools"; // should match same directory in ServerRunner 
	public static final String SOURCE = "/etc/source/";

	public static final String TETKO_ANACONDA = "/anaconda3/";
	// libraries for conversion of molecules - N.B.! to be fixed later -- code is non reliable above

	public static final String FINISHED = "Finished"; // report that task is finished to browser
	public static final String NOSESSION = "Unknown session ID:";

	public static final String MISSED_DESCRIPTOR = "MISSED_DESCRIPTOR";
	public static final String ERROR_SMILES = "An invalid or empty molecule";
	public static final String SUCCESS = "success";
	public static final String FAILED = "failed";
	public static final String ERROR_ONE_ATOM = "one atom only";
	public static final String CHEMENGINE = "ChemEngine";
	public static final String DUMMY = "Dummy";

	public static final String ENV = "OCHEM_ENV";

	public static final double MAXVALUE = 99999; // if condition values are large, the user should use another Unit scale

	public static final int VARCHAR = 255; // to avoid problem in database due to the length of session

	public static final int WEBSERVICE_SIZE_STORED_MOLECULE = 500; // maximum number of molecules that are processed and stored in OCHEM

	public static final String EXTENDED_USER = "qspr.entities.ExtendedUser";
	public static final String DEFAULT_USER = "qspr.entities.User";

	public static final String PYTHON36 = "/opt/conda/envs/map4/bin/python3.6";
	public static final String RDKITPYTHON ="export RDBASE=/opt/conda/envs/map4/share/RDKit/; export PYTHONPATH=/opt/conda/envs/map4/include/rdkit; " + PYTHON36;

}
