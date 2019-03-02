
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
//import ambit2.
import utils.FileUtilities;

import org.openscience.cdk.smiles.SmilesGenerator;

import wishartlab.cfmid_plus.fragmentation.FPLists;
import wishartlab.cfmid_plus.fragmentation.FragmentationPattern;
import wishartlab.cfmid_plus.fragmentation.Fragmenter;
import wishartlab.cfmid_plus.molecules.StructuralClass.ClassName;
import wishartlab.cfmid_plus.molecules.StructuralClass;
import wishartlab.cfmid_plus.molecules.StructureExplorer;

/**
 * Hello world!
 *
 */
public class RuleBasedFrag 
{
	
	public RuleBasedFrag(){
		
	}
	static String VERSION = "1.0.7";
	
	private static Options generateOptions(){
		
		final Option smiInputOption = Option.builder("ismi")
				.required(false)
				.hasArg(true)
				.argName("SMILES Input")
				.longOpt("ismiles")
				.desc("The input, which can be a SMILES string")
				.build();
		
		final Option sdfInputOption = Option.builder("isdf")
				.required(false)
				.hasArg(true)
				.argName("Sdf Input")
				.longOpt("sdfinput")
				.desc("Accept the input as a SD file. You must enter an output filename")
				.build();
		
		final Option adductTypes = Option.builder("a")
				.required(false)
				.hasArg(true)
				.argName("Adduct types")
				.longOpt("adducts")
				.desc("A semicolon-separated list of adduct types to consider for fragmentation. ")
				.build();	
		
		final Option outputDestinationOption = Option.builder("o")
				.required(false)
				.hasArg(true)
				.argName("Spectra Destination")
				.longOpt("output")
				.desc("Destination for the created spectra files. If -ismi is select and the output is specified, this must be a file. If -sdf is selected and the output is specified, this must be an existing folder.")
				.build();		
	
		final Option helpOption = Option.builder("h")
				.required(false)
				.hasArg(false)
				.argName("help")
				.longOpt("help")
				.desc("Prints the usage.")
				.build();	
		
		final Options options = new Options();
		options.addOption(smiInputOption);
		options.addOption(sdfInputOption);		
		options.addOption(adductTypes);
		options.addOption(outputDestinationOption);
		options.addOption(helpOption);
		
		return options;
	}
	
	public static CommandLine generateCommandLine(
			final Options options, final String[] commandLineArguments) throws ParseException{
			final CommandLineParser cmdLineParser = new DefaultParser();
			CommandLine commandLine = null;
			
			String header ="This is msrb-fragmenter-"+VERSION + ". It is a tool that uses a rule-based fragmentation algorithm to predict ESI-MS/MS spectra at 10eV, 20eV, and 40 eV."
					+ " The library of fragmentation rules currently covers 24 classes of lipids and seven adduct types.";
			String footer ="";			
			HelpFormatter formatter = new HelpFormatter();
			
			try{
				commandLine = cmdLineParser.parse(options, commandLineArguments);
			}
			catch (MissingOptionException missingOptionException){
				
				if( Arrays.asList(commandLineArguments).contains("-h") || Arrays.asList(commandLineArguments).contains("--help")){
					formatter.printHelp("\njava -jar msrb-fragmenter-"+VERSION + ".jar --help", header, options, footer, true);
				}
				else {
					System.out.println(missingOptionException.getLocalizedMessage());
				}			
			}
			catch (ParseException parseException){
				System.out.println("Could not parse the command line arguments "
						+ Arrays.toString(commandLineArguments) + "\nfor the following reaons:" 
						+ parseException);
			}			
			
			return commandLine;
		}	
	
	
    public static void main( String[] args ) throws Exception
    {
		Options options = generateOptions();
		CommandLine commandLine = generateCommandLine(options, args);

//		System.out.println(commandLine.getOptionValue("b"));
//		System.out.println(commandLine.getOptionValue("f").length());
//		System.out.println(commandLine.getOptionValue("m"));
		
		String iFormat = null;
		IAtomContainer singleInput = null;
		String inputFileName = null;   	
		String adducts = null;
    	String outputName	= null;
    	


//    	if(!(Arrays.asList(args).contains("-ismi") || Arrays.asList(args)
//				.contains("--ismiles") || Arrays.asList(args).contains("-isdf") 
//				|| Arrays.asList(args).contains("--sdfinput"))){
//    		throw new MissingOptionException("You must specify a format");
 //   	}
//    	
    	
    	if(Arrays.asList(args).contains("-ismi") || Arrays.asList(args)
    			.contains("--ismiles")){
    		iFormat = "smiles";
    	}
    	
    	else if(Arrays.asList(args).contains("-isdf") || Arrays.asList(args)
				.contains("--sdfinput")){
			iFormat = "sdf";
		}
    	
    	if(Arrays.asList(args).contains("-a")|| Arrays.asList(args)
				.contains("--adducts")){
    		adducts = commandLine.getOptionValue("a");
	    	if( adducts == null){
	    		throw new MissingOptionException("Missing argument for option 'a'. Using this option requires that the user specifies a list of adducts. "
	    				+ "If the '-a' option is not used, then only selected adducts will be considered, based on the query compound's chemical class.");
	    	}
    	}
    	
    	if(Arrays.asList(args).contains("-o")|| Arrays.asList(args)
				.contains("--output")){
    		outputName = commandLine.getOptionValue("o");
    	}
    	
		if(commandLine !=null){
			
			IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
			SmilesParser	smiParser		= new SmilesParser(builder);
	        IChemObjectBuilder bldr = DefaultChemObjectBuilder.getInstance();
	        Fragmenter fr = new Fragmenter();	
	        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();


	        if(adducts == null || adducts.contains("all;") || adducts.contains(";all")){
				if(commandLine.getOptionValue("ismi") != null){
					String molSmiles = commandLine.getOptionValue("ismi");
					singleInput = smiParser.parseSmiles(molSmiles.replace("[O-]", "O"));
					
					if(outputName == null){
						InChIGenerator gen = factory.getInChIGenerator(singleInput);
						outputName = gen.getInchiKey()+".log";
					}
					
					try{
						
						fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName);
					}
		            catch(NullPointerException e){
		            	System.err.println("Could not compute spectra for " + molSmiles);
		            	System.err.println(e.getMessage());
		            }

				}
				else if(commandLine.getOptionValue("isdf") != null){
					
					inputFileName = commandLine.getOptionValue("isdf");
					if(inputFileName == null){
						throw new MissingOptionException("You must be specify an input file name (Molfile or SDF). For more information, type java -jar biotransformer-1.0.8 --help.");
					} else{				
						File file = new File(inputFileName);
						if((!file.exists()) & file.isDirectory()){
							throw new IllegalArgumentException("Invalid argument: Please make sure to enter a valid existing directory name if you select the -isdf or -isdfInput option.");
						}else{
							
							if(outputName != null){
								File directory_ = new File(outputName);
								if(!(file.exists() & directory_.exists())){
									throw new IllegalArgumentException("Invalid argument: Please make sure to enter a valid existing directory name if you select the -isdf or -isdfInput option.");
								}
								
							}else{
								outputName = new File("").getAbsolutePath();
							}
							
							IAtomContainerSet atContainers = FileUtilities.parseSdfAndAddTitles(inputFileName, factory);
							int counter = 0;
							
							for(IAtomContainer atc : atContainers.atomContainers()){
								counter++;
								try{
									fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName+"/"+ atc.getProperty(CDKConstants.TITLE) + ".log");
								}
					            catch(NullPointerException e){
					            	System.err.println("Could not compute spectra for molecule no. " + counter + "(" + atc.getProperty(CDKConstants.TITLE) + ")");
					            	System.err.println(e.getMessage());
					            }								
							}
	
						}
					}
					
					
				}
					
				}
				else{					
					ArrayList<String> adduct_list = new ArrayList<String>();
					if(adducts != null){
						adduct_list = new ArrayList<String>(Arrays.asList(adducts.split(";")));
					}
					
					if(commandLine.getOptionValue("ismi") != null){
						String molSmiles = commandLine.getOptionValue("ismi");
						singleInput = smiParser.parseSmiles(molSmiles.replace("[O-]", "O"));
						
						if(outputName == null){
							InChIGenerator gen = factory.getInChIGenerator(singleInput);
							outputName = gen.getInchiKey()+".log";
						}
						//System.out.println("NIL ATOM CONTAINER: " + (singleInput.isEmpty()));
						//SmilesGenerator sg = new SmilesGenerator().unique();
						//System.out.println(sg.create(singleInput));
//						try{
							fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName, adduct_list);
//						}
//			            catch(NullPointerException e){
//			            	System.err.println("Could not compute spectra for " + molSmiles);
//			            	System.err.println(e.getMessage());
//			            }

					}
					else if(commandLine.getOptionValue("isdf") != null){
						inputFileName = commandLine.getOptionValue("isdf");
						if(inputFileName == null){
							throw new MissingOptionException("You must be specify an input file name (Molfile or SDF). For more information, type java -jar biotransformer-1.0.8 --help.");
						} else{				
							File file = new File(inputFileName);
							if((!file.exists()) & file.isDirectory()){
								throw new IllegalArgumentException("Invalid argument: Please make sure to enter a valid existing directory name if you select the -isdf or -isdfInput option.");
							}else{
								
								if(outputName != null){
									File directory_ = new File(outputName);
									if(!(file.exists() & directory_.exists())){
										throw new IllegalArgumentException("Invalid argument: Please make sure to enter a valid existing directory name if you select the -isdf or -isdfInput option.");
									}
									
								}else{
									outputName = new File("").getAbsolutePath();
								}
								
								IAtomContainerSet atContainers = FileUtilities.parseSdfAndAddTitles(inputFileName, factory);
								int counter = 0;
								
								for(IAtomContainer atc : atContainers.atomContainers()){
									counter++;
									try{
										fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName+"/"+ atc.getProperty(CDKConstants.TITLE) + ".log");
										fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName+"/"+ atc.getProperty(CDKConstants.TITLE) + ".log", adduct_list);
									
									}
						            catch(NullPointerException e){
						            	System.err.println("Could not compute spectra for molecule no. " + counter + "(" + atc.getProperty(CDKConstants.TITLE) + ")");
						            	System.err.println(e.getMessage());
						            }								
								}
		
							}
						}

					}					
				}
			}
	
		
		
		
		
		
//    	String molSmiles 	= args[0];
//    	String adductType 	= args[1];
//    	String outputName	= args[1];
//    			
//  	
//	
//  	SmilesParser sp =  new SmilesParser(SilentChemObjectBuilder.getInstance());
//        SmilesGenerator sg = new SmilesGenerator().unique();
//       StructureExplorer se = new StructureExplorer();
//        IAtomContainer molecule = sp.parseSmiles(molSmiles.replace("[O-]", "O"));
//        IChemObjectBuilder bldr = DefaultChemObjectBuilder.getInstance();
//        Fragmenter fr = new Fragmenter();
//        
//		IAtomContainer standardized_mol = se.standardizeMolecule(molecule);
//		StructuralClass.ClassName type = se.findClassName(standardized_mol);
//		System.out.println("The type of this molecule is " + String.valueOf(type));
//       
//        if(FPLists.classSpecificFragmentationPatterns.containsKey(type)){
//            try {
//            	fr.saveSingleCfmidLikeMSPeakList(standardized_mol, bldr, type, outputName, false);
//            }
//            catch(NullPointerException e){
//            	System.err.println("Could not compute spectra for " + args[0]);
//            	System.err.println(e.getMessage());
//            }        	
//        } 
//        else{
//        	System.err.println("Could not compute spectra for " + args[0]);
//        	System.err.println("The compound does not belong to any of the covered lipid classes.");
//       }
        

        
        
        
////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O[C@H]3[C@@H]([C@H]([C@H]([C@H](O3)CO)O)O)NC(=O)C)O[C@@]4(C[C@@H]([C@H](C(O4)[C@@H]([C@@H](CO)O)O)NC(=O)C)O)C(=O)O)O)O)O)[C@@H](/C=C/CCCCCCCCCCCCC)O");
////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@H]1C([C@H]([C@@H](C(O1)CO)O[C@H]2C([C@H]([C@H](C(O2)CO)O)O[C@@]3(CC([C@H](C(O3)[C@@H]([C@@H](CO)O)O)NC(=O)C)O)C(=O)O)O)O)O)[C@@H](/C=C/CCCCCCCCCCCCC)O");
////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1OC(CO)[C@H](O)[C@H](O[C@@]2(CC(O)[C@@H](NC(C)=O)C(O2)[C@H](O)[C@H](O)CO)C(O)=O)C1O)[C@H](O)C=CCCCCCCCCCCCCC");
//
////        IAtomContainer molecule = sp.parseSmiles("OCC(O)COP(O)(=O)OCC(CO(C(=O)C))O(C(=O)CC)");
//        	IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@H](CO)O)OC(=O)CCCCCCCCC/C=C\\CCCCCC");
//        	IAtomContainer molecule = sp.parseSmiles("CCCCC(=O)OC[C@H](COP(O)(=O)OC[C@@H](O)CO)OC(=O)CCCC");
//        	IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCC(O)COP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC");
//        IAtomContainer molecule = sp.parseSmiles("[H]O[C@]([H])(C([H])([H])OC([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])C([H])([H])OP(=O)(O[H])OC([H])([H])C([H])([H])N([H])[H]");


        
//        System.out.println("Molecule standardized: " + sg.create(se.standardizeMolecule(molecule)) + "\n\n");


//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);		
//		System.out.println("MOLECULE WITH EXPLICIT HYDROGENS: " + sg.create(molecule));      
      

//        fr.saveSingleCfmidLikeMSPeakList(molecule, bldr, adductType, outputName);
//        System.out.println("Bla bla");
    
    
    
        
    
    
    
    }
	
	
	
//    public static void main( String[] args ) throws Exception
//    {
//    
//    	
//        SmilesParser sp =  new SmilesParser(SilentChemObjectBuilder.getInstance());
//        SmilesGenerator sg = new SmilesGenerator().unique();
//        StructureExplorer se = new StructureExplorer();
//        
//////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O[C@H]3[C@@H]([C@H]([C@H]([C@H](O3)CO)O)O)NC(=O)C)O[C@@]4(C[C@@H]([C@H](C(O4)[C@@H]([C@@H](CO)O)O)NC(=O)C)O)C(=O)O)O)O)O)[C@@H](/C=C/CCCCCCCCCCCCC)O");
//////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@H]1C([C@H]([C@@H](C(O1)CO)O[C@H]2C([C@H]([C@H](C(O2)CO)O)O[C@@]3(CC([C@H](C(O3)[C@@H]([C@@H](CO)O)O)NC(=O)C)O)C(=O)O)O)O)O)[C@@H](/C=C/CCCCCCCCCCCCC)O");
//////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1OC(CO)[C@H](O)[C@H](O[C@@]2(CC(O)[C@@H](NC(C)=O)C(O2)[C@H](O)[C@H](O)CO)C(O)=O)C1O)[C@H](O)C=CCCCCCCCCCCCCC");
////
//////        IAtomContainer molecule = sp.parseSmiles("OCC(O)COP(O)(=O)OCC(CO(C(=O)C))O(C(=O)CC)");
////        	IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@H](CO)O)OC(=O)CCCCCCCCC/C=C\\CCCCCC");
////        	IAtomContainer molecule = sp.parseSmiles("CCCCC(=O)OC[C@H](COP(O)(=O)OC[C@@H](O)CO)OC(=O)CCCC");
////        	IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCC(O)COP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC");
//        	IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC(O)COP([O-])(=O)OC[C@@H](COC(=O)CCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC");
//        	IAtomContainer molecule2 = sp.parseSmiles("CCCCCCCC\\C=C/CCCCCCCC(=O)O[C@H]1CC[C@]2(C)C3CC[C@]4(C)[C@H](CC[C@H]4C3CCC2=C1)[C@H](C)CCCC(C)C");
//        	
//        	IAtomContainer molecule3 = sp.parseSmiles("[O-]C(CCCCCCCCCCCCCCC)=O");
//        	
//        	System.out.println("Molecule3 standardized: " + sg.create(se.standardizeMolecule(molecule3)) + "\n\n");
//        	
//        	
////         System.out.println("[#6:7]-[#8;X2:1][P;X4:2]([#8;X1-:4])(=[O;X1:5])[#8;X2:3]-[#6:6]".replaceAll("\\:[0-9]+\\]", "\\]"));
////         IAtomContainer molecule_c = se.standardizeMolecule(molecule);
////         AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule_c);
////         CDKHydrogenAdder.getInstance(molecule_c.getBuilder()).addImplicitHydrogens(molecule_c);
////         String molecule_cs = sg.create(molecule_c);
////         System.out.println("Standardized molecule: "+ molecule_cs);
//        	
//        // CCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCC(O)COP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC
//        //        
////        String molecule_c = sg.create(StructureExplorer.cleanMolecule(molecule,new SMIRKSManager(SilentChemObjectBuilder.getInstance())));
////        System.out.println("Cleaned molecule: "+ molecule_c);
////        System.out.println(AtomContainerManipulator.getNaturalExactMass(molecule));
////        System.out.println("Type of this molecule: " + StructureExplorer.fin dClassName(molecule));
//        Fragmenter fr = new Fragmenter();
//
//    	
//
//        
//////        IAtomContainerSet fragments = fr.fragmentMolecule(molecule);
//////        
//////        for (int i = 0; i < fragments.getAtomContainerCount(); i++){
//////        	System.out.println(sg.create(fragments.getAtomContainer(i)));
//////        }
////        LinkedHashMap<String,ArrayList<String>> peaksResults = fr.generateMSPeakList(molecule2);
////        ArrayList<String> peaks  =  peaksResults;
////        fr.generateGraphFromPeakListAndSaveToDir(peaksResults);   
////        System.out.println("WRITING  PEAKS");
////        for(Map.Entry<String, ArrayList<String>> p : peaksResults.entrySet() ){
////        	for (int k =0; k < p.getValue().size(); k++){
////        		System.out.println(p.getValue().get(k));
////        	}
////        }
//        
//        System.out.println(FragmentationPattern.patterns.get(ClassName._1_ALKYL_GLYCEROPHOSPHOETHANOLAMINES).get("[M+H]+")[0][1]);
//        
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
//		
//		System.out.println("MOLECULE WITH EXPLICIT HYDROGENS: " + sg.create(molecule));      
//        
//        IChemObjectBuilder bldr = DefaultChemObjectBuilder.getInstance();
////        ArrayList<LinkedHashMap<String,ArrayList<String>>> allPeaks = fr.generatePeakListFromSDF("data/lipids.sdf", bldr);
////        fr.generateGraphFromPeakListAndSaveToDir(allPeaks.get(0));
//              
////         fr.generateSpectraFromPeakListAndSaveToDir(allPeaks);
//
//        
////      IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
////      fr.savePeakListFromSDF("data/ymdb_test.sdf",bldr,"data/ymdb_spectra.msp"); 
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/test3.sdf",bldr,"[M+H]+");
////		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/DG_IDS/DG_IDS.sdf",bldr);
////		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/TG_IDS/TG_IDS.sdf",bldr);
////		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PA_IDS/PA_IDS.sdf",bldr);
////		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/Lyso-PC_IDS/Lyso-PC_IDS.sdf",bldr);
////		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/Lyso-PE_IDS/Lyso-PE_IDS.sdf",bldr);
////		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/MG_IDS/MG_IDS.sdf",bldr); 
////	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PI_IDS/PI_IDS.sdf",bldr);
////	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PG_IDS/PG_IDS.sdf",bldr);
////	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PE_IDS/PE_IDS.sdf",bldr);
////	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PS_IDS/PS_IDS.sdf",bldr);   
////	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PC_IDS/PC_IDS.sdf",bldr);
////	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/CL_IDS/CL_IDS.sdf",bldr);
//      
//      
//      
//       
//
//
//      
//        
////        System.out.println(f.size());
////        System.out.println(f.get(f.size()-1));
////        
////        
////        System.out.println("\nH: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("H")));
////        System.out.println("Li: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("[Li]")));
////        System.out.println("H2O: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("[H]O[H]")));
////        System.out.println("O: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("O")));
////        System.out.println("C: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("C")));
////        System.out.println(sp.parseSmiles("C"));
////        System.out.println(AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("CC(=O)N[C@@H]1[C@@H](O)C[C@@H](O[C@H]1[C@H](O)[C@H](O)CO)C(O)=O")));
////        System.out.println(AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("OCC(O)CO")));
////        
////        //196.41083318690977
//////        10.0794075382579
////        
////        
////        System.out.println("C6H10O5: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("OC1CC(O)C(O)C(O)C1O")));
////        
//
//        
////        fr.generateMSSpecImage();
//    }
}

