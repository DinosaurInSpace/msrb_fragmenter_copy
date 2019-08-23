
//import java.io.BufferedReader;
//import java.io.BufferedWriter;
import java.io.File;
//import java.io.FileReader;
//import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
//import java.util.LinkedHashMap;
//import java.util.Map;

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
//import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
//import org.openscience.cdk.tools.CDKHydrogenAdder;
//import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
//
//import ambit2.smarts.SMIRKSManager;
//import ambit2.
import utils.FileUtilities;

//import org.openscience.cdk.smiles.SmilesGenerator;
//
//import wishartlab.cfmid_plus.fragmentation.FPLists;
//import wishartlab.cfmid_plus.fragmentation.FragmentationPattern;
import wishartlab.cfmid_plus.fragmentation.Fragmenter;
//import wishartlab.cfmid_plus.molecules.StructuralClass.ClassName;
//import wishartlab.cfmid_plus.molecules.StructuralClass;
import wishartlab.cfmid_plus.molecules.StructureExplorer;

/**
 * Hello world!
 *
 */
public class RuleBasedFrag 
{
	
	public RuleBasedFrag(){
		
	}
	static String VERSION = "1.0.8";
	
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
		
		final Option classifyOption = Option.builder("c")
				.required(false)
				.hasArg(false)
				.argName("Classify")
				.longOpt("classify")
				.desc("Find the class of a compound (Limited to 20 chemical categories so far.")
				.build();		

		final Option nopredictOption = Option.builder("n")
				.required(false)
				.hasArg(false)
				.argName("No Predict")
				.longOpt("nopredict")
				.desc("Do not predict spectra. This will work only if the classify option is selected.")
				.build();	
		
//		final Option versionOption = Option.builder("v")
//				.required(false)
//				.hasArg(false)
//				.argName("version")
//				.longOpt("version")
//				.desc("Prints the version.")
//				.build();
		
		final Options options = new Options();
		options.addOption(smiInputOption);
		options.addOption(sdfInputOption);		
		options.addOption(adductTypes);
		options.addOption(outputDestinationOption);
		options.addOption(helpOption);
		options.addOption(classifyOption);
		options.addOption(nopredictOption);
		
		return options;
	}
	
	public static CommandLine generateCommandLine(
			final Options options, final String[] commandLineArguments) throws ParseException{
			final CommandLineParser cmdLineParser = new DefaultParser();
			CommandLine commandLine = null;
			
			String header ="This is version "+ VERSION +" of msrb-fragmenter. " + ". It is a tool that uses a rule-based fragmentation algorithm to predict ESI-MS/MS spectra at 10eV, 20eV, and 40 eV."
					+ " The library of fragmentation rules currently covers 21 classes of lipids and seven adduct types.";
			String footer ="For more information/help, contact cfmid@wishartlab.com.";
			
			HelpFormatter formatter = new HelpFormatter();
			
			try{
				commandLine = cmdLineParser.parse(options, commandLineArguments);
				if( Arrays.asList(commandLineArguments).contains("-h") || Arrays.asList(commandLineArguments).contains("--help")){
					formatter.printHelp("\njava -jar msrb-fragmenter.jar --help", header, options, footer, true);
				}
			}
			catch (MissingOptionException missingOptionException){
				
				if( Arrays.asList(commandLineArguments).contains("-h") || Arrays.asList(commandLineArguments).contains("--help")){
					formatter.printHelp("\njava -jar msrb-fragmenter.jar --help", header, options, footer, true);
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
    	boolean nopredict = false;
    	boolean classify = false;
    	
    	


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
    	
    	if((Arrays.asList(args).contains("-n") || Arrays.asList(args).contains("--nopredict")) &&
    			!(Arrays.asList(args).contains("-c") || Arrays.asList(args).contains("--classify")
    			)) {
    		throw new MissingOptionException("No task selected. You must select at least one of the two tasks: MS-spectra prediction or classify. "
    										+ "To only classify, select both the -n and -c options.\n");
    	}
    	
    	if(Arrays.asList(args).contains("-o")|| Arrays.asList(args)
				.contains("--output")){
    		outputName = commandLine.getOptionValue("o");
    	}
    	
    	if(Arrays.asList(args).contains("-n") || Arrays.asList(args).contains("--nopredict")){
    		nopredict = true;
    	}
    	
    	if(Arrays.asList(args).contains("-c") || Arrays.asList(args).contains("--classify")){
    		classify=true;
    	}
    	
    	
		if(commandLine !=null){
			
			IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
			SmilesParser	smiParser		= new SmilesParser(builder);
	        IChemObjectBuilder bldr = DefaultChemObjectBuilder.getInstance();
	        Fragmenter fr = new Fragmenter();	
	        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
	        
	        
	        if(classify){
	        	
				if(commandLine.getOptionValue("ismi") != null){
					String molSmiles = commandLine.getOptionValue("ismi");
					singleInput = smiParser.parseSmiles(molSmiles.replace("[O-]", "O"));
					
					
					try{						
						System.out.println("Chemical Class = " + StructureExplorer.findClassName(singleInput));
					}
		            catch(NullPointerException e){
		            	System.err.println("Could not classify the compound wth smiles " + molSmiles);
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
								outputName = new File("data").getAbsolutePath();
							}
							
							IAtomContainerSet atContainers = FileUtilities.parseSdfAndAddTitles(inputFileName, factory);
							int counter = 0;
							
							for(IAtomContainer atc : atContainers.atomContainers()){
								counter++;
//								String title = atc.getProperty(CDKConstants.TITLE);
								try{
									System.out.println("MOLECULE " + counter + " - Chemical Class: " + StructureExplorer.findClassName(atc));
//									if(title != null){
//										System.out.println(title + " - Chemical Class : " + StructureExplorer.findClassName(atc));
//									}
//									else{
//										System.out.println(counter + " - Chemical Class : " + StructureExplorer.findClassName(atc));
//									}
									
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
	        
	        if(nopredict == false){
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
//								System.out.println(inputFileName);
								if(outputName != null){
									File directory_ = new File(outputName);
									if(!(file.exists() & directory_.exists())){
										throw new IllegalArgumentException("Invalid argument: Please make sure to enter a valid existing directory name if you select the -isdf or -isdfInput option.");
									}
									
								}else{
									outputName = new File("data").getAbsolutePath();
								}
								
								IAtomContainerSet atContainers = FileUtilities.parseSdfAndAddTitles(inputFileName, factory);
								int counter = 0;
								
								for(IAtomContainer atc : atContainers.atomContainers()){
//									System.out.println(atc.isEmpty());
//									System.out.println(atc == null);
//									System.out.println(atc);
									counter++;
									System.out.println("\nMOLECULE " + counter);
									System.out.println("TITLE: " + atc.getProperty(CDKConstants.TITLE));
									try{
										String title = String.valueOf(atc.getProperty(CDKConstants.TITLE)).replace("/", "_");
										fr.saveSingleCfmidLikeMSPeakList(atc, bldr, outputName+"/"+ title + ".log");
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
							try{
								fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName, adduct_list);
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
										outputName = new File("data").getAbsolutePath();
									}
									
									IAtomContainerSet atContainers = FileUtilities.parseSdfAndAddTitles(inputFileName, factory);
									int counter = 0;
									
									System.out.println(atContainers.getAtomContainerCount());
									for(IAtomContainer atc : atContainers.atomContainers()){
										counter++;
										System.out.println("MOLECULE " + counter);
										System.out.println("TITLE: " + atc.getProperty(CDKConstants.TITLE));
										try{
											String title = String.valueOf(atc.getProperty(CDKConstants.TITLE)).replace("/", "_");
//											fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName+"/"+ atc.getProperty(CDKConstants.TITLE) + ".log");
											fr.saveSingleCfmidLikeMSPeakList(singleInput, bldr, outputName+"/"+ title + ".log", adduct_list);
//										
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
	
		}
    }

}

