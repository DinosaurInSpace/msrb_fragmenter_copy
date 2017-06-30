/**
 * 
 */
package wishartlab.cfmid_plus.fragmentation;

/**
 * @author Yannick Djoumbou Feunang
 *
 */

import wishartlab.cfmid_plus.molecules.StructureExplorer;
import wishartlab.cfmid_plus.molecules.StructuralClass;
import wishartlab.cfmid_plus.molecules.StructuralClass.ClassName;
import wishartlab.cfmid_plus.fragmentation.FragmentationPattern;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.query.SMARTSException;



//import org.rosuda.REngine.*;
//import org.rosuda.REngine.Rserve.*;

import org.rosuda.REngine.REXP;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.RList;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;


import org.rosuda.REngine.Rserve.*;
//import RserveEngine;

import org.math.R.Rsession;



/**
 * This class describes the Fragmenter, which takes a molecule and returns the possible fragments according
 * to more elaborate fragmentation patterns
 * 
 * References:
 * R1: Murphy RC et al (2001) Chem. Rev. 2001, 101, 479−526
 * 
 * 
 * @author yandj
 *
 */

/**
 */



public class Fragmenter {

	
	IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
	private SMIRKSManager smrkMan = new SMIRKSManager(bldr);


	private AtomContainerManipulator acm = new AtomContainerManipulator();
	private SmilesParser sParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
	private SmilesGenerator sGen = new SmilesGenerator().isomeric();
	private String host;
	private int port;
	StructureExplorer sExplorer = new StructureExplorer();

	public void Fragmenter(){
		this.host = "127.0.0.1";
		this.port = 6311;
		this.smrkMan.setFlagFilterEquivalentMappings(true);
	}
	
	public void Fragmenter(String host, int port){
		this.host = host;
		this.port = port;
		this.smrkMan.setFlagFilterEquivalentMappings(true);
	}
	
	
	public static void main(String[] args) throws Exception{
		// DG(18:0/20:4) - "[H]C([H])(O)C([H])(OC(=O)CCC\C=C/C\C=C/C\C=C/C\C=C/CCCCC)C([H])([H])OC(=O)CCCCCCCCCCCCCCCCC"
		// MG(18:2/0:0/0:0) - "CCCCCC=CCC=CCCCCCCCC(=O)OCC(O)CO"
		// MG(16:0/0:0/0:0) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO"
		// TG(16:0/16:0/16:0) - "CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
		// TG(20:0/20:1/20:0) - "[H]C([H])(OC(=O)CCCCCCCCCCCCCCCCCCC)C([H])(OC(=O)CCCCCCCCC\C=C/CCCCCCCC)C([H])([H])OC(=O)CCCCCCCCCCCCCCCCCCC"
		// TG(22:6/22:6/22:6) - "CCC=CCC=CCC=CCC=CCC=CCC=CCCC(=O)OCC(COC(=O)CCC=CCC=CCC=CCC=CCC=CCC=CCC)OC(=O)CCC=CCC=CCC=CCC=CCC=CCC=CCC"		
		// PA(16:0/18:1(9Z)) - "OP(O)(=O)OCC(CO(C(=O)CCCCCCCCCCCCCCC))O(C(=O)CCCCCCCC=CCCCCCCCC)"
		// PE(16:0/18:1(9Z)) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OCCN)OC(=O)CCCCCCCC=CCCCCCCCC"
		// PE(17:0/14:1(9Z)) - "NCCOP(O)(=O)OC[C@@H](CO(C(=O)CCCCCCCCCCCCCCCC))O(C(=O)CCCCCCCC=CCCCC)"
		// PE(16:0/20:4(5Z,8Z,11Z,14Z)) - "CCCCCC=CCC=CCC=CCC=CCCCC(=O)O[C@@H](COP(O)(=O)OCCN)COC(=O)CCCCCCCCCCCCCCC"
		// PE(17:0/20:4(5Z,8Z,11Z,14Z)) - "CCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCCN)OC(=O)CCCC=CCC=CCC=CCC=CCCCCC"
		// PC(16:0/18:1(9Z)) - "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC\C=C/CCCCCCCC"
		// PC(16:0/0:0) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCC[N+](C)(C)C"
		// PC(18:1(9Z):0) - "CCCCCCCCC=CCCCCCCCC(=O)OCC(O)COP(O)(=O)OCC[N+](C)(C)C"
		// PA(18:1(9Z):0) - "CCCCCCCC/C=C\CCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O"
		// PS(16:0/18:1(9Z)) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OC[C@@H](N)C(O)=O)OC(=O)CCCCCCCC=CCCCCCCCC"
		// Cer(18:1(9Z)/16:0) - "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCC)=O)CO"
		// Cer(18:1(9Z)/OH-18:0) - "CCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\C=C\CCCCCCCCCCCCC"
		// SM(18:1(9Z):16:0) - "CCCCCCCCCCCCCCCC(=O)N[C@H](COP(O)(=O)OCC[N+](C)(C)C)[C@@H](O)C=CCCCCCCCCCCCCC"
		// CL(1'-[18:1(9Z)/18:1(9Z)],3'-[18:1(9Z)/18:1(9Z)]) - "[H]C(O)(COP(O)(=O)OC[C@@]([H])(COC(=O)CCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCC=CCCCCCCCC)COP(O)(=O)OC[C@@]([H])(COC(=O)CCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCC=CCCCCCCCC"
		// Carnitine 16 - "CCCCCCCCCCCCCCCC(=O)OC(CC(O)=O)C[N+](C)(C)C"
		
    	String molSmiles = "CCCCCCCCCCCCCCCC(=O)OC(CC(O)=O)C[N+](C)(C)C";    	
    	String adductType 	= "[M+]";
    	String outputName	= "data/Carnitine_16_plus.log";
    	
    	Fragmenter fr = new Fragmenter();

        
        StructureExplorer se = new StructureExplorer();
        IAtomContainer molecule = fr.sParser.parseSmiles(molSmiles);
        
       MolecularFormulaManipulator mfm = new MolecularFormulaManipulator();
       IMolecularFormula m = mfm.getMolecularFormula(molecule);
       System.out.println("Natural excat mass: " + AtomContainerManipulator.getNaturalExactMass(molecule));
       System.out.println( "Major isotopse mass: " + MolecularFormulaManipulator.getMajorIsotopeMass(m));
//       System.out.println( "Total exact mass: " +  AtomContainerManipulator.getTotalExactMass(molecule));
        
        
//        LinkedHashMap<String, IAtomContainer> fragments = fr.fragmentMolecule(molecule, ClassName._1_MONOACYLGLYCEROLS, adductType);
//
        
//        for(Map.Entry<String, IAtomContainer> fragment : fragments.entrySet()){
//        	System.out.println(fragment.getKey() + " : " + fr.sGen.create(fragment.getValue()));
//        }        
//        LinkedHashMap<String, Double> masses = fr.computeFragmentMasses(fragments);
//        System.out.println(masses);

        fr.saveSingleCfmidLikeMSPeakList(molecule, fr.bldr, adductType, outputName);
        System.out.println(AtomContainerManipulator.getTotalFormalCharge(molecule));
	}
	
	
	public LinkedHashMap<String, IAtomContainer>  fragmentMolecule(IAtomContainer molecule, StructuralClass.ClassName type) throws Exception{
		return fragmentMolecule(molecule, type, "[M+H]+");
	}
	
	public LinkedHashMap<String, IAtomContainer> fragmentMolecule(IAtomContainer molecule, StructuralClass.ClassName type,
			String adductType) throws Exception{	
		LinkedHashMap<String, IAtomContainer> fragments = new LinkedHashMap<String, IAtomContainer>();

		System.out.println(FPLists.classSpecificFragmentationPatterns);
		System.out.println(FPLists.classSpecificFragmentationPatterns.get(type).get(adductType));
		System.out.println(FPLists.classSpecificFragmentationPatterns.get(type).get(adductType)
				.patterns.keySet());
				
		if(type != StructuralClass.ClassName.NIL && FPLists.classSpecificFragmentationPatterns.get(type).containsKey(adductType)){
//			IAtomContainer stMol = this.sExplorer.standardizeMolecule(molecule);
			SmilesGenerator sg = new SmilesGenerator().isomeric();
			
			for(Map.Entry<String, String[]> pattern : FPLists.classSpecificFragmentationPatterns.get(type).get(adductType)
					.patterns.entrySet()){
				SMIRKSReaction sReaction = this.smrkMan.parse(pattern.getValue()[1]);
				System.out.println("pattern.getKey(): " + pattern.getKey());
				System.out.println(sReaction == null);
				IAtomContainerSet current_set = this.smrkMan.applyTransformationWithSingleCopyForEachPos(
						molecule, null, sReaction);
				
				if(current_set != null){
					
					 System.out.println("Nr of metabolites: " + String.valueOf(current_set.getAtomContainerCount()));	
					 System.out.println("Nr of unique metabolites: " + String.valueOf(this.sExplorer.uniquefy(current_set).getAtomContainerCount()));					 
					 
					 IAtomContainerSet current_set_unique = this.sExplorer.uniquefy(current_set);
					 
					 for(IAtomContainer atc : current_set_unique.atomContainers()){
						 IAtomContainerSet partitions = this.sExplorer.partition(atc);

						 	for(IAtomContainer partition : partitions.atomContainers()){
						 		if(this.sExplorer.containsSmartsPattern(partition, pattern.getValue()[0])){
//								if(this.sExplorer.containsSmartsPattern(this.sExplorer.standardizeMolecule(partition), pattern.getValue()[0])){
									if(pattern.getValue().length==2){
										fragments.put(pattern.getKey(), AtomContainerManipulator.removeHydrogens(partition));
										System.out.println(sg.create(partition) + " was added.");
										
									}else if (pattern.getValue().length==3){
										// take the SMILES of the fragment an add the adduct, whiich would be the 
										// third element in the array.
										IAtomContainer  adjustedPartition = this.sParser.parseSmiles(sg.create(AtomContainerManipulator.removeHydrogens(partition))
												+ "." + pattern.getValue()[2]);	
										fragments.put(pattern.getKey(), adjustedPartition);
										System.out.println(sg.create(adjustedPartition) + " was added.");		
									}

									break;	
								} else{
									System.out.println(sg.create(partition) + " was NOT added.");
								}			
						 	}
						}
					} else{
						
					}
				
				}

				return fragments;
			} else
				System.err.println("WHAT IS GOING ON?");
				return null;
		
		}
	
	public LinkedHashMap<String, Double> computeFragmentMasses(LinkedHashMap<String, IAtomContainer> fragments){
		LinkedHashMap<String, Double> fmasses = new LinkedHashMap<String, Double>();
		
		if(fragments != null){
			for(Map.Entry<String, IAtomContainer> fragment :  fragments.entrySet()){
				System.err.println(fragment.getKey());
				fmasses.put(fragment.getKey(), 
						Math.floor(StructureExplorer.getMajorIsotopeMass(fragment.getValue()) * 100000)/100000);		
			}
		}else{
			fmasses = null;
		}
		
		return fmasses;		
	}
	public LinkedHashMap<String, Double> computeFragmentMassChargeRatios(LinkedHashMap<String, IAtomContainer> fragments) throws CDKException{
		LinkedHashMap<String, Double> fmassToChargeRatios = new LinkedHashMap<String, Double>();
		
		if(fragments != null){
			for(Map.Entry<String, IAtomContainer> fragment :  fragments.entrySet()){
				System.err.println(fragment.getKey());
				System.err.println(this.sGen.create(fragment.getValue()));
				int z = (int) Math.abs(AtomContainerManipulator.getTotalFormalCharge(fragment.getValue()));
				System.err.println("Z = " + z);
				if(z==0){
					fmassToChargeRatios.put(fragment.getKey(), 
							Math.floor(StructureExplorer.getMajorIsotopeMass(fragment.getValue()) * 100000)/100000);						
				} else{
					fmassToChargeRatios.put(fragment.getKey(), 
							Math.floor(StructureExplorer.getMajorIsotopeMass(fragment.getValue()) * 100000)/ (z * 100000));	
				}
	
			}
		}else{
			fmassToChargeRatios = null;
		}
		
		return fmassToChargeRatios;		
	}	
	
	public LinkedHashMap<String, ArrayList<String>>  annotatePeakList(LinkedHashMap<String, IAtomContainer> fragments ,ClassName type, FragmentationCondition fragCondition) throws CDKException{
		
		if(fragments == null || type == null || fragCondition == null){
			return null;
		} else{
			LinkedHashMap<String, ArrayList<String>> results = new LinkedHashMap<String, ArrayList<String>>();
			ArrayList<String> frag_smiles_mass = new ArrayList<String>();
//			LinkedHashMap<Integer,ArrayList<String>> peaks = new LinkedHashMap<Integer,ArrayList<String>>();
			ArrayList<String> peaks = new ArrayList<String>();
			
			String s_peaks = "";
			
//			ArrayList<String> metadata = new ArrayList<String>();
//			ArrayList<String> molName = new ArrayList<String>();
//			ArrayList<String> sAdductType = new ArrayList<String>();
//			sAdductType.add(fragCondition.adductName);
//			String name = molecule.getProperty(CDKConstants.TITLE);
//			
//			if(name == null){
//				molName.add("");
//			}else {
//				molName.add(name);
//			}
			
			
			
			LinkedHashMap<String, Double> fragmentMassChargeRatios = computeFragmentMassChargeRatios(fragments);
//			LinkedHashMap<String, Double> fragmentMasses = computeFragmentMasses(fragments);
			LinkedHashMap<Double, ArrayList<String>> massesToLabel= new LinkedHashMap<Double, ArrayList<String>>();
			fragCondition.print();
			ArrayList<String> labels = new ArrayList<String>(fragmentMassChargeRatios.keySet());
			MSPeakRelativeAbundance mra = MSPRelativeAbundanceList.
					classSpecificRelativeAbundances.get(type).get(fragCondition.adductName + "_" + fragCondition.collisionEnergy);
			
			for(Map.Entry<String, Double> fm : fragmentMassChargeRatios.entrySet()){
				if(mra.getRelativeAbundances().containsKey(fm.getKey())){
					if(massesToLabel.containsKey(fm.getValue())){
						massesToLabel.get(fm.getValue()).add(fm.getKey());
					}else{
						massesToLabel.put(fm.getValue(), new ArrayList<String>());
						massesToLabel.get(fm.getValue()).add(fm.getKey());
					}
				}
			}


			int index = 0;
			mra.print();
			System.out.println("fragmentMasses: " + fragmentMassChargeRatios);
			System.out.println("mra.getRelativeAbundances(): " + mra.getRelativeAbundances());
			
//			for(String label : mra.getRelativeAbundances().keySet()){
//				Double mass = fragmentMasses.get(label.trim());
//				
//				System.out.println(label);
//				System.out.println(labels.indexOf(label));
//				System.out.println(mass);
//				System.out.println(massesToLabel.get(mass).size());
//				
//				
//				if( massesToLabel.get(mass).size() == 1 ){
//					peaks.add(String.format("%.5f %4f %3d (1)",mass, mra.getRelativeAbundances().get(label), labels.indexOf(label)));
//				} else {
//					ArrayList<String> indexes = new ArrayList<String>();
//					ArrayList<String> scores = new ArrayList<String>();
//					for(String x : massesToLabel.get(mass)){
//						indexes.add(String.valueOf(labels.indexOf(label)));
//						scores.add("1");
//					}
//					peaks.add(String.format("%.5f %4f %8s",mass, mra.getRelativeAbundances().get(label), 
//							String.join(" ", indexes)) + " (" + String.join(" ", scores) + ")");
//					
//				}
//				
//				
//			}
			
			for(Double mass : massesToLabel.keySet()){
//				if(){
//				boolean relevantPeaks = true;

					if( massesToLabel.get(mass).size() == 1 && mra.getRelativeAbundances().containsKey(massesToLabel.get(mass).get(0))){
						peaks.add(String.format("%.5f %4f %3d (1)",mass, mra.getRelativeAbundances().get(massesToLabel.get(mass).get(0)), labels.indexOf(massesToLabel.get(mass).get(0))));
					} else{
						
						
						ArrayList<String> indexes = new ArrayList<String>();
						ArrayList<String> scores = new ArrayList<String>();
						for(String x : massesToLabel.get(mass)){
							if(mra.getRelativeAbundances().containsKey(x)){
								indexes.add(String.valueOf(labels.indexOf(x)));
								scores.add("1");
							}
						}
						peaks.add(String.format("%.5f %4f %8s",mass, mra.getRelativeAbundances().get(massesToLabel.get(mass).get(0)), 
								String.join(" ", indexes)) + " (" + String.join(" ", scores) + ")");
					}
//				}
			}
			
			
////			fragCondition.print();
//			switch (type) {
//				case _1_MONOACYLGLYCEROLS:
//					
//					System.out.println(MSPRelativeAbundanceList.
//							classSpecificRelativeAbundances.get(ClassName._1_MONOACYLGLYCEROLS).get("[M+Li]+_10"));
//					
//					MSPeakRelativeAbundance mra_1_mg_li = MSPRelativeAbundanceList.
//							classSpecificRelativeAbundances.get(ClassName._1_MONOACYLGLYCEROLS).get(
//									fragCondition.getAdductName() + "_"+ fragCondition.getCollisionEnergy());				
//					int index = 0;
////					for(String label : mra.getRelativeAbundances().keySet()){
////						System.out.println(label);
////						peaks.add(String.format("%.5f %4f %2d (1)",(fragmentMasses.get(label)), index, mra.getRelativeAbundances().get(label)));
////						System.out.println(String.format("%.5f %4f %2d (1)",(fragmentMasses.get(label)), index,mra.getRelativeAbundances().get(label)));
////						index++;
////					}
//					
//					for(String label : fragments.keySet()){
//						if(mra_1_mg_li.getRelativeAbundances().containsKey(label)){
//							System.out.println(label);
//							
//							if((label == "[M+H]-sn1-H2O" ||  label == "[M+H]-sn2-H2O") && 
//									fragmentMasses.get("[M+H]-sn1-H2O") == fragmentMasses.get("[M+H]-sn2-H2O")) {
//								
//							}
//							peaks.add(String.format("%.5f %4f %2d (1)",(fragmentMasses.get(label)), mra_1_mg_li.getRelativeAbundances().get(label), index));
//							System.out.println(String.format("%.5f %4f %2d (1)",(fragmentMasses.get(label)), mra_1_mg_li.getRelativeAbundances().get(label), index));
//							index++;
//						}
//					}
//					
//				
//				case _2_MONOACYLGLYCEROLS:
//					
//					
//				case PHOSPHATIDIC_ACIDS:
//					
//					MSPeakRelativeAbundance mra_pa_h = MSPRelativeAbundanceList.
//					classSpecificRelativeAbundances.get(ClassName._1_MONOACYLGLYCEROLS).get(
//							fragCondition.getAdductName() + "_"+ fragCondition.getCollisionEnergy());
//					
//					int indexPA = 0;
//					
//					
//					
//					
//					
//					
//				/*
//				 * Add Fragment masses and structures
//				 */
//					
				for(Map.Entry<String, IAtomContainer> frag : fragments.entrySet()){
					frag_smiles_mass.add(String.format("%.5f", fragmentMassChargeRatios.get(frag.getKey())) + " " + this.sGen.create(fragments.get(frag.getKey())));
					System.out.println(String.format("%.5f", fragmentMassChargeRatios.get(frag.getKey())) + " " + this.sGen.create(fragments.get(frag.getKey())));
				}
//					
//
//////			s_peaks = String.format("MW: %.5f", mass + 7.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 7.0) + "\n" + 
//////					"Comment: "  + String.format("Parent= %.5f",  mass + 7.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
//////			metadata.add(s_peaks);	
//	////
//////			break
//			
//			
//			}
						
			
			results.put("peaks_list", peaks);
//			results.put("metadata", metadata);
			results.put("fragments", frag_smiles_mass);
			return results;
		}
		

	}
	
	public LinkedHashMap<String, ArrayList<String>> generateCfmidLikeMSPeakList(IAtomContainer molecule, FragmentationCondition fragCondition) throws Exception{
		LinkedHashMap<String, ArrayList<String>> annotatedPeaks = new LinkedHashMap<String, ArrayList<String>>();
		
		if(fragCondition.getCollisionEnergy() == 10 || fragCondition.getCollisionEnergy() == 20 || fragCondition.getCollisionEnergy() == 40) {
			IAtomContainer standardized_mol = this.sExplorer.standardizeMolecule(molecule);
			System.out.println("The standardized molecule is " + this.sGen.create(standardized_mol));
			StructuralClass.ClassName type = this.sExplorer.findClassName(standardized_mol);
			System.out.println("The type of this molecule is " + String.valueOf(type));
			LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(standardized_mol, type, fragCondition.getAdductName());
			
			return  annotatePeakList(fragments ,type, new FragmentationCondition(fragCondition.getAdductName(), 10));
		} else {
			
			throw new IllegalArgumentException("The collision energy must be either 10 eV, 20 eV, or 40 eV.\nPlease enter a valid collision energy");
		}
	}
	

	
	public LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> generateCfmidLikeMSPeakList(IAtomContainer molecule, String adductType) throws Exception{
		LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> annotatedPeaks = new LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>>();
		
		IAtomContainer standardized_mol = this.sExplorer.standardizeMolecule(molecule);
		StructuralClass.ClassName type = this.sExplorer.findClassName(standardized_mol);
		System.out.println("The type of this molecule is " + String.valueOf(type));
		LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(standardized_mol, type, adductType);
		FragmentationCondition fragCondition_10 =  new FragmentationCondition(adductType, 10);
		FragmentationCondition fragCondition_20 =  new FragmentationCondition(adductType, 20);
		FragmentationCondition fragCondition_40 =  new FragmentationCondition(adductType, 40);
		
		annotatedPeaks.put(10, annotatePeakList(fragments ,type, fragCondition_10));
		annotatedPeaks.put(20,annotatePeakList(fragments ,type, fragCondition_20));
		annotatedPeaks.put(40,annotatePeakList(fragments ,type, fragCondition_40));
		
		return annotatedPeaks;
	}

//	public LinkedHashMap<String, ArrayList<String>> generateCfmidLikeMSPeakList(IAtomContainer molecule, String adductType) throws Exception{
//
//		// System.out.println("THIS IS THE ADDUCT TYPE: " + adductType);
//		LinkedHashMap<String, ArrayList<String>> results = new LinkedHashMap<String, ArrayList<String>>();
//		ArrayList<String> peaks = new ArrayList<String>();
//		String s_peaks = "";
//		ArrayList<String> metadata = new ArrayList<String>();
//		ArrayList<String> molName = new ArrayList<String>();
////		ArrayList<String> adductType = new ArrayList<String>();
//		ArrayList<String> sAdductType = new ArrayList<String>();
//		sAdductType.add(adductType);
//		String name = molecule.getProperty(CDKConstants.TITLE);
//		
//		if(name == null){
//			molName.add("");
//		}else {
//			molName.add(name);
//		}
//			
//		
//		
//		if(adductType == null){
//			adductType = "[M+H]+";
//		}		
//		System.out.println("ADDUCT TYPE: " + adductType);
//        
////    	System.out.println("isFlagApplyStereoTransformation(): " + this.smrkMan.isFlagApplyStereoTransformation());
////    	System.out.println("isFlagFilterEquivalentMappings() : " + this.smrkMan.isFlagFilterEquivalentMappings());
////    	System.out.println("isFlagCheckResultStereo() : " + this.smrkMan.isFlagCheckResultStereo());
////    	System.out.println("isFlagProcessResultStructures() : " + this.smrkMan.isFlagProcessResultStructures());	
////		smrkMan.setFlagApplyStereoTransformation(true); 
//		this.smrkMan.setFlagFilterEquivalentMappings(true);
////		this.smrkMan.setFlagCheckResultStereo(true);
////		smrkMan.setFlagProcessResultStructures(true);
////    	System.out.println("isFlagApplyStereoTransformation(): " + this.smrkMan.isFlagApplyStereoTransformation());
////    	System.out.println("isFlagFilterEquivalentMappings() : " + this.smrkMan.isFlagFilterEquivalentMappings());
////    	System.out.println("isFlagCheckResultStereo() : " + this.smrkMan.isFlagCheckResultStereo());
////    	System.out.println("isFlagProcessResultStructures() : " + this.smrkMan.isFlagProcessResultStructures());
//		
//		
//
//		IAtomContainer standardized_mol = this.sExplorer.standardizeMolecule(molecule);
//		StructuralClass.ClassName type = this.sExplorer.findClassName(standardized_mol);
//		System.out.println("The type of this molecule is " + String.valueOf(type));
//		
////		IAtomContainer cleaned_mol = this.sExplorer.cleanMolecule(molecule, this.smrkMan);
//
//		IAtomContainerSet fragments = fragmentMolecule(standardized_mol,type, adductType);
//		System.out.println("Number of fragments: " + fragments.getAtomContainerCount());
//		// ArrayList of SMILES and mass for each fragment
//		ArrayList<String> frag_smiles_mass = new ArrayList<String>();
//		SmilesGenerator sg = new SmilesGenerator().unique();
//		System.out.println(sg.create(standardized_mol));
//		
////		for(IAtomContainer ac : fragments.atomContainers()){
////			System.out.println(sg.create(ac));
////		}
//		
//
//		double mass = acm.getNaturalExactMass(standardized_mol);
//		
//		double sn1_mass;
//		double sn2_mass;
//		double sn3_mass;
//		double sn4_mass;
//		double sn5_mass;
//		double sn6_mass;
//		double sn7_mass;
//		double sn8_mass;
//		double sn9_mass;
//		double sn10_mass;
//		double sn11_mass;
//		double sn12_mass;
//		double sn13_mass;
//		double electron_mass = 0.00054858026;
//		
//		String c1_smiles, c2_smiles, c3_smiles, c4_smiles, c5_smiles, c6_smiles, c7_smiles, c8_smiles, c9_smiles, c10_smiles,
//		c11_smiles, c12_smiles;
//		
////		NumberFormat mNumberFormat = NumberFormat.getInstance();
////		mNumberFormat.setMinimumFractionDigits(5);
////		mNumberFormat.setMaximumFractionDigits(5);
//		
////		System.out.println(mass);
//		if(fragments.getAtomContainerCount()>0) {
////			switch (type) {
////				case PHOSPHATIDYLCHOLINES:
////									
//////					System.out.println(fragments.getAtomContainerCount());		
//////					System.err.println("Adduct type is...: "  + adductType);
////					// System.err.println("adductType.contentEquals(\"[M+H]+\") " + adductType.contentEquals("[M+H]+"));
////					// System.err.println("adductType == [M+H]+: " + adductType == "[M+H]+");
////					
////					
////					if(adductType.contentEquals("[M+H]+")){
////						
//////						for(int l = 0; l < fragments.getAtomContainerCount(); l++){
//////							System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
//////						}
////						
////						// [M+H]-H2O (-18)
////						sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////		
////						// [M+H]-C3H9N (-59)
////						sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////						
////						// [M+H]-sn2
////						sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////						
////						// [M+H]-sn2-H2O
////						sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////						
////						// [M+H]-sn1
////						sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////						
////						// [M+H]-C5H14NO4P (-183)
////						sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////						
////						// [M+H]-sn1-H2O
////						sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////		
////						frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////						frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////						frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////						frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////						frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////						frag_smiles_mass.add(String.format("%.5f", sn6_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(5)));
////						frag_smiles_mass.add(String.format("%.5f", sn7_mass) + " " + sg.create(fragments.getAtomContainer(6)));					
////						
//////						System.err.println("frag_smiles_mass");
//////						System.err.println(frag_smiles_mass);
////						
////						if(sn3_mass == sn5_mass){
////						
////							peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 400.0/999.0));
////							peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 1.0), 1.0));
////							peaks.add(String.format("%.5f %5f 2 4 (1 1)",sn3_mass, 600.0/999.0));
////							peaks.add(String.format("%.5f %5f 3 6 (1 1)",sn4_mass, 600.0/999.0));
////							peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass, 400.0/999.0));
////						
////						}else{
////							
////							peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 400.0/999.0));
////							//System.out.println(peaks.size());
////							peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass), 1.0));
////							//System.out.println(peaks.size());
////							peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass, 600.0/999.0));
////							//System.out.println(peaks.size());
////							peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass, 600.0/999.0));
////							//System.out.println(peaks.size());
////							peaks.add(String.format("%.5f %5f 4 (1)",sn5_mass, 600.0/999.0));
////							//System.out.println(peaks.size());
////							peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass, 400.0/999.0));
////							//System.out.println(peaks.size());
////							peaks.add(String.format("%.5f %5f 6 (1)",sn7_mass, 600.0/999.0));
////							//System.out.println(peaks.size());
////						
////						
////						}
////		
////						s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////								"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; [M+H]+" + "\n" + String.format("Num Peaks: %2d",peaks.size());
////						
////					} else
////						if(adductType.contentEquals("[M+Na]+")){
//////							for(int l = 0; l < fragments.getAtomContainerCount(); l++){
//////								System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
//////							}
////
////							// [M+Na]-C3H9N (-59)
////							sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////							
////							// [M+Na]-C5H14NO4P (-183)
////							sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////							
////							// [M+Na]-sn1
////							sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////							
////							// [M+Na]-sn2
////							sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////							
////							// [M+Na]-sn1-H2O
////							sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////							
////							// [M+Na]-sn2-H2O
////							sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////							
////							// [M+Na]-59-sn1
////							sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////							
////							// [M+Na]-59-sn2
////							sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////							
////							
////							frag_smiles_mass.add(String.format("%.5f", sn1_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(0)));
////							frag_smiles_mass.add(String.format("%.5f", sn2_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(1)));
////							frag_smiles_mass.add(String.format("%.5f", sn3_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(2)));
////							frag_smiles_mass.add(String.format("%.5f", sn4_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(3)));
////							frag_smiles_mass.add(String.format("%.5f", sn5_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(4)));
////							frag_smiles_mass.add(String.format("%.5f", sn6_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(5)));
////							frag_smiles_mass.add(String.format("%.5f", sn7_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(6)));
////							frag_smiles_mass.add(String.format("%.5f", sn8_mass + 23.0) + " " + sg.create(fragments.getAtomContainer(7)));
////							
////							
////							if(sn3_mass == sn4_mass){
////								
////								peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 23.0), 1.0));
////								peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 23.0), 600.0/999.0));
////								peaks.add(String.format("%.5f %5f 2 3 (1 1)",sn3_mass + 23.0, 20.0/999.0));
////								peaks.add(String.format("%.5f %5f 4 5 (1 1)",sn5_mass + 23.0, 20.0/999.0));
////								peaks.add(String.format("%.5f %5f 6 (1)",sn7_mass + 23.0, 40.0/999.0));
////							
////							}else{
////							
////								peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 23.0), 999.0/999.0));
////								peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 23.0), 600.0/999.0));
////								peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 23.0, 20.0/999.0));
////								peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + 23.0, 20.0/999.0));
////								peaks.add(String.format("%.5f %5f 4 (1)",sn5_mass + 23.0, 20.0/999.0));
////								peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass + 23.0, 20.0/999.0));
////								peaks.add(String.format("%.5f %5f 6 (1)",sn7_mass + 23.0, 40.0/999.0));
////								peaks.add(String.format("%.5f %5f 6 (1)",sn8_mass + 23.0, 40.0/999.0));
////
////							}
////							
////							s_peaks = String.format("MW: %.5f", mass + 23.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 23.0) + "\n" + 
////									"Comment: " + String.format("Parent= %.5f", mass + 23.0)  + "; [M+Na]+" + "\n" + String.format("Num Peaks: %2d",peaks.size());
////						}
////					
////					//System.err.println("s_peaks: \n"+ s_peaks);
////					metadata.add(s_peaks);
////					
////					break;
////	
////				
////				case _1_LYSOPHOSPHATIDYLCHOLINES:
////					
//////					System.out.println(fragments.getAtomContainerCount());
//////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
//////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
//////					}
////					// adductType.add("[M+H]+");
////					// [M+H]-H2O (-18)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					
////					// [M+H]-sn1
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					
////					// [M+H]-sn1-H2O
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					
////					// C5H15NO4P m/z=184
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;			
////	
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////	
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 1.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass, 10.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass, 10.0/999.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass, 300.0/999.0));
////					
////					s_peaks = String.format("MW: %.5f", mass + 1.0  ) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					
////					break;
////	
////	//			case _2_LYSOPHOSPHATIDYLCHOLINES:
////	//				
////	//				sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////	//
////	//				System.out.println(sn2_mass);
////	//
////	//				frag_smiles_mass.add(sn2_mass + "\t" + sg.create(fragments.getAtomContainer(0)));
////	//
////	//				peaks.add(String.format("%.5f %4d \"[M+H]-H2O (-18)\"",(mass - 18.0), 999));
////	//				peaks.add(String.format("%.5f %4d \"[M+H]-sn2\"",(mass - sn2_mass + 1), 10));
////	//				peaks.add(String.format("%.5f %4d \"[M+H]-sn2-H20\"",(mass - sn2_mass - 18 + 1), 10));
////	//				
////	//				 // 184.15087533864045 was generated by creating a molecule of phosphocholine "C[N+](C)(C)CCOP(O)(O)=O"
////	//				// with CDK and calulating the natural exact mass with AtomContainerManipulator
////	//				peaks.add(String.format("%.5f %4d \"C5H15NO4P m/z=184\"",(184.15087533864045), 300));
////	//				
////	//				s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////	//						"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////	//				metadata.add(s_peaks);					
////	//				break;
////	//				
////				case PLASMENYLPHOSPHATIDYLCHOLINES:
////					
//////					System.out.println(fragments.getAtomContainerCount());
//////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
//////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
//////					}
////					
////					// adductType.add("[M+H]+");
////					// [M+H]-C5H14NO4P (-183)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// [M+H]-sn2
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// [M+H]-sn1 (alkenyl ether loss)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					// [M+H]-sn2-H2O
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					// [M+H]-sn1-C3H9N (-59)
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////	
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////					
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 10.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass), 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass), 999.0/999.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass), 10.0/999.0));	
////					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass), 100.0/999.0));
////	
////	
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////					System.out.println(s_peaks);
////	
////					break;
////					
////				case PHOSPHATIDYLETHANOLAMINES:
////					
////					// e.g.:[M+H]+; GPEtn(10:0/11:0)
////					// [M+H]-C2H8NO4P (-141); [M+H]-sn1; [M+H]-sn2; [M+H]-sn1-H2O; [M+H]-sn2-H2O; sn2-O; sn1-O
////					
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M+H]+");
////					// [M+H]-H2O (-18)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// [M+H]-sn1
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					// [M+H]-sn2
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					// [M+H]-sn1-H2O
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					// [M+H]-sn2-H2O
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////					// [M+H]-C2H8NO4P (-141)
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// sn2-O
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////					// sn1-O
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////	
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(5)));	
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass) + " " + sg.create(fragments.getAtomContainer(6)));
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass) + " " + sg.create(fragments.getAtomContainer(7)));
////					
////					// LipidBlast reports conflicting m/z values for [M+H]-sn1, [M+H]-sn1-H2O, [M+H]-sn2, and [M+H]-sn2-H2O
////					
////					if(sn8_mass == sn7_mass){
////						peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 2 (1 1)",(sn2_mass + 1.0), 1.0));
////						peaks.add(String.format("%.5f %5f 3 4 (1 1)",(sn4_mass), 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)",(sn6_mass), 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 7 (1 1)",(sn7_mass), 1.0/999.0));			
////					
////					}else{	
////						peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 1.0), 1.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 1.0, 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + 1.0, 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)",sn5_mass + 1.0, 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass + 1.0, 1.0));
////						peaks.add(String.format("%.5f %5f 6 (1)",sn7_mass, 1.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)",sn8_mass, 1.0/999.0));					
////					}
////					
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass - 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);			
////					break;
////					
////					
////				// e.g.: LysoPE   PE(16:0/0:0)
////				case _1_LYSOPHOSPHATIDYLETHANOLAMINES:
////				
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M+H]+");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000; // [M+H]-H2O (-18)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000; // [M+H]-(C2NH5+H20) (-61)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000; // [M+H]-C3H8NO6P (-141)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000; // [M+H]-sn1
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000; // [M+H]-sn1-H2O
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 1.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 1.0), 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass + 1.0), 1.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass + 1.0), 10.0/999.0));	
////					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass + 1.0), 10.0/999.0));
////	
////	
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////					System.out.println(s_peaks);
////					break;				
////				
////				// e.g.: LysoPE   PE(P-16:0/0:0) or 1-1Z-alkenylphosphoethanolamine OR _1_LYSOPLASMENYL_ETHANOLAMINES
////				case  _1_1Z_ALKYENYL_GLYCEROPHOSPHOETHANOLAMINES:
////				
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M+H]+");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000; // [M+H]-H2O (-18)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000; // [M+H]-(C2NH5+H20) (-61)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000; // [M+H]-C3H8NO6P (-141)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000; // [M+H]-sn1
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000; // [M+H]-sn1-H2O
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 1.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 1.0), 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass + 1.0), 1.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass + 1.0), 10.0/999.0));	
////					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass + 1.0), 10.0/999.0));
////	
////	
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////					System.out.println(s_peaks);
////					break;
////					
////					
////				//	 e.g.: LysoPE   PE(O-16:0/0:0) or 1-1Z-alkylphosphoethanolamine	OR _1_LYSOPLASMANYL_ETHANOLAMINES			
////				case _1_ALKYL_GLYCEROPHOSPHOETHANOLAMINES:
////					
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M+H]+");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000; // [M+H]-H2O (-18)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000; // [M+H]-(C2NH5+H20) (-61)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000; // [M+H]-C3H8NO6P (-141)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000; // [M+H]-sn1
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000; // [M+H]-sn1-H2O
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 1.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 1.0), 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass + 1.0), 1.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass + 1.0), 10.0/999.0));	
////					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass + 1.0), 10.0/999.0));
////	
////	
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////					System.out.println(s_peaks);				
////				break;				
////					
////	//			case _2_LYSOPHOSPHATIDYLETHANOLAMINES:
////	//				
////	//				System.out.println(fragments.getAtomContainerCount());
////	//				for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////	//					System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////	//				}
////	//				// adductType.add("[M+H]+");
////	//				
////	//				sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////	//				sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////	//				sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////	//				sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////	//				sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////	//				
////	//				frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(5)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn5_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(7)));
////	//				
////	//				peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 1.0));
////	//				peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 1.0), 200.0/999.0));
////	//				peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass + 1.0), 999.0/999.0));
////	//				peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass + 1.0), 10.0/999.0));	
////	//				peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass + 1.0), 10.0));
////	//
////	//
////	//				s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////	//						"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////	//				metadata.add(s_peaks);	
////	//				System.out.println(s_peaks);
////	//				break;			
////					
////				// e.g.: plasmenyl-PE 18:0; [M+H]+; PE(P-16:0/2:0)
////				case PLASMENYLPHOSPHATIDYLETHANOLAMINES:
////					
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M+H]+");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000; // [M+H]-H2O (-18)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000; //sn1 ether + C2H8NO3P (+124)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000; // [M+H]-C2H8NO4P (-141)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000; // sn1 ether +C2H8NO3P-H3PO4
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000; // [M+H]-C2H8NO4P-sn2
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////					
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 100.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass), 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass), 30.0/999.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass), 150.0/999.0));	
////					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass), 200.0/999.0));
////	
////	
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////					System.out.println(s_peaks);
////					
////					break;				
////				
////				
////				
////				case PHOSPHATIDYLSERINES:
////					// GPSer(2:0/2:0)
////					
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M+H]+");
////					
////					//[M+H]-sn1
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					//[M+H]-sn1-H2O
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// [M+H]-C3H5NO2 (-87)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////			
////					//[M+H]-sn2
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					//[M+H]-sn2-H2O
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					//[M+H]-C3H8NO6P (-185)
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////	
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 2.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 2.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass) + " " + sg.create(fragments.getAtomContainer(5)));
////					
////					if(sn1_mass == sn4_mass){
////						peaks.add(String.format("%.5f %5f 0 3 (1)",(sn1_mass + 2.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 4 (1)",sn2_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 2.0, 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass, 999.0/999.0));					
////					}else{
////						
////						peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 2.0, 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + 2.0, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)",sn5_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass, 999.0/999.0));					
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);
////					
////					break;
////					
////					
////				case SPHINGOMYELINS:
////					// E.g.: SM(d14:0/2:0)
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M]+");
////					
////					// [M+]-H2O (-18)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// [M+]-C3H9N (-59)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// [M+] (-18 -59)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;			
////					// [M+]-C5H14NO4P (-183)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////		
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + electron_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + electron_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + electron_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + electron_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////	
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + electron_mass), 999.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass + electron_mass, 10.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + electron_mass, 50.0/999.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + electron_mass, 20.0/999.0));
////	
////					s_peaks = String.format("MW: %.5f", mass + electron_mass) + "\n" + String.format("PRECURSORMZ: %.5f",mass + electron_mass) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f", mass + electron_mass)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);
////	
////					break;
////				
////				case PHOSPHATIDIC_ACIDS:
////	
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M+Na2-H]+");
////					
////					// PROBLEM HERE: MASS IS 2 less than reported by lipidblast. m/z are ca. 3 less.
////					// [M+Na2-H]-sn1
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// [M+Na2-H]-sn1-H
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// [M+Na2-H]-sn2
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////					// [M+Na2-H]-sn2-H
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////				
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 46.0 - 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 46.0 - 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass  + 46.0 - 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 46.0 - 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					
////					if(sn1_mass == sn3_mass){
////						peaks.add(String.format("%.5f %5f 0 2 (1)",(sn1_mass + 46.0 - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 3 (1)",(sn2_mass + 46.0 - 1.0), 999.0/999.0));
////				
////					}else{
////						
////						peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 46.0 - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass + 46.0 - 1.0), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass + 46.0 - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass + 46.0 - 1.0), 999.0/999.0));				
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass + 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1.0) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f", mass + 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);
////				
////					break;
////					
////				case _1_LYSOPHOSPHATIDIC_ACIDS:
////	
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////	
////					// adductType.add("[M+H]+");
////					
////					//[M+H]-H2O (-18)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					//[M+H]-sn1
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// frag. C3H7O5P (155)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////					//[M+H]-H2O-C3H7O5P
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////	
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass + 1.0, 50.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 1.0, 999.0/999.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + 1.0, 500.0/999.0));
////					
////					s_peaks = String.format("MW: %.5f", mass + 1) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f",  mass + 1)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					break;
////					
////				case _1_O_ALKENYL_GLYCEROPHOSPHATES:
////	
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////	
////					// adductType.add("[M+H]+");
////					
////					//[M+H]-H2O (-18)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					//[M+H]-sn1
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// frag. C3H7O5P (155)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////					//[M+H]-H2O-C3H7O5P
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////	
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass + 1.0, 50.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 1.0, 999.0/999.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + 1.0, 500.0/999.0));
////					
////					s_peaks = String.format("MW: %.5f", mass + 1) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f",  mass + 1)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					break;				
////					
////	//			case _2_LYSOPHOSPHATIDIC_ACIDS:
////	//				
////	//				sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////	//				sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////	//				
////	//				System.out.println(mass + "\t" + sn1_mass + "\t" + sn2_mass);
////	//				
////	//				frag_smiles_mass.add(sn1_mass + "\t" + sg.create(fragments.getAtomContainer(0)));
////	//				frag_smiles_mass.add(sn2_mass + "\t" + sg.create(fragments.getAtomContainer(1)));
////	//	
////	//				peaks.add(String.format("%.5f %4d \"[M+H]-H20 (-18)\"",(mass + 1 - 18),200));
////	//				peaks.add(String.format("%.5f %4d \"[M+H]-sn2\"",(mass + 1 - sn2_mass),50));
////	//				peaks.add(String.format("%.5f %4d \"frag. C3H7O5P (155)\"", sn1_mass, 999));
////	//				peaks.add(String.format("%.5f %4d \"[M+H]-H2O-C3H7O5P\"",(mass + 1 - 18 - 155),500));
////	//	
////	//				
////	//				
////	//				s_peaks = String.format("MW: %.5f", mass + 1) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1) + "\n" + 
////	//						"Comment: "  + String.format("Parent= %.5f",  mass + 1)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////	//				metadata.add(s_peaks);				
////	//				break;									
////	//				
////				case PHOSPHATIDYLINOSITOLS:
////					
////					// GPIns(2:0/2:0)
////					// [M-H]-sn1; [M-H]-sn2; [M-H]-sn1-H2O; [M-H]-sn2-H2O; [M-H]-sn1-C6H12O6; [M-H]-sn2-C6H12O6; sn2 FA; sn1 FA
////					
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M-H]-");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000; // [M-H]-sn1
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000; // [M-H]-sn2
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000; // [M-H]-sn1-H2O
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000; // [M-H]-sn2-H2O
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000; // [M-H]-sn1-C6H12O6
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000; // [M-H]-sn2-C6H12O6
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000; // sn1 FA
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000; // sn2 FA
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass) + " " + sg.create(fragments.getAtomContainer(5)));
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass) + " " + sg.create(fragments.getAtomContainer(6)));
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass) + " " + sg.create(fragments.getAtomContainer(7)));
////					
////					if(sn7_mass == sn8_mass){
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)",(sn1_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)",sn5_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 5 (1 1)",sn3_mass, 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 7 (1 1)",sn7_mass, 999.0/999.0));					
////	
////					}else{
////						peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass, 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)",sn5_mass, 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass, 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 (1)",sn7_mass, 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)",sn8_mass, 999.0/999.0));
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass - 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					
////					break;
////					
////					
////				case PHOSPHATIDYLGLYCEROLS:
////					
////					// Name: PG 12:0; [M-H]-; GPGro(2:0/10:0)
////					
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M-H]-");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;  // [M-H]-sn1
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;  // [M-H]-sn1-H2O
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;  // [M-H]-sn1-C3H8O3				
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	// [M-H]-sn2		
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;  // [M-H]-sn2-H2O
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000; // [M-H]-sn2-C3H8O3				
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000; // sn2 FA
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;	// sn1 FA
////	
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass) + " " + sg.create(fragments.getAtomContainer(5)));
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass) + " " + sg.create(fragments.getAtomContainer(7)));
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass) + " " + sg.create(fragments.getAtomContainer(6)));
////					
////					if(sn7_mass == sn8_mass){
////						peaks.add(String.format("%.5f %5f 0 3 (1 1)",(sn1_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 4 (1 1)",(sn2_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 5 (1 1)",(sn3_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 7 (1 1)",(sn7_mass), 999.0/999.0));
////					}else{
////						
////						peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)",sn5_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)",sn6_mass, 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 (1)",sn7_mass, 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)",sn8_mass, 999.0/999.0));
////						
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment: " + String.format("Parent= %.5f", mass - 1.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);
////					
////					break;
////							
////					
////				case CARDIOLIPINS:
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M-2H](2-)");
////					
////					// [M-2H](2-)
////					// In most cases, the [M−2H]2− ion is more abundant than the [M−H]− ion. In contrast to ESI-MS analysis, MALDI-MS analysis
////					// in the negative-ion mode shows predominant singly charged deprotonated CL species
////					// Comprehensive Mass Spectrometry of Lipids. Chap 7.8; p.188.
////					// 0 2 4 7 11 14 17 20 26 25 23 24
////					
////									
////					// [M-2H](2-) -sn1 (monolyso)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// [M-2H](2-) -sn3 (monolyso)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;	
////					// [M-2H](2-) -(sn1+sn2 dilyso)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////					// [M-2H](2-) -(sn3+sn4 dilyso)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					// [sn3+C3H6PO4-H](-) (137.00)
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					// [sn4+C3H6PO4-H](-) (137.00)
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;	
////					// [sn1+C3H6PO4-H](-) (137.00)
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;	
////					// [sn2+C3H6PO4-H](-) (137.00)
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////					// sn3 FA
////					sn9_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(8)) * 100000)/100000;
////					// sn4 FA
////					sn10_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(9)) * 100000)/100000;	
////					// sn1 FA
////					sn11_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(10)) * 100000)/100000;	
////					// sn2 FA
////					sn12_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(11)) * 100000)/100000;
////									
////					c1_smiles = sg.create(fragments.getAtomContainer(0));
////					c2_smiles = sg.create(fragments.getAtomContainer(1));
////					c3_smiles = sg.create(fragments.getAtomContainer(2));
////					c4_smiles = sg.create(fragments.getAtomContainer(3));
////					c5_smiles = sg.create(fragments.getAtomContainer(4));
////					c6_smiles = sg.create(fragments.getAtomContainer(5));
////					c7_smiles = sg.create(fragments.getAtomContainer(6));
////					c8_smiles = sg.create(fragments.getAtomContainer(7));
////					c9_smiles = sg.create(fragments.getAtomContainer(8));
////					c10_smiles = sg.create(fragments.getAtomContainer(9));
////					c11_smiles = sg.create(fragments.getAtomContainer(10));
////					c12_smiles = sg.create(fragments.getAtomContainer(11));
////	
////					frag_smiles_mass.add(String.format("%.5f", (sn1_mass)/2.0) + " " + c1_smiles);
////					frag_smiles_mass.add(String.format("%.5f", (sn2_mass)/2.0) + " " + c2_smiles);
////					frag_smiles_mass.add(String.format("%.5f", (sn3_mass)/2.0) + " " + c3_smiles);
////					frag_smiles_mass.add(String.format("%.5f", (sn4_mass)/2.0) + " " + c4_smiles);
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + c5_smiles);
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass) + " " + c6_smiles);
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass) + " " + c7_smiles);
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass) + " " + c8_smiles);
////					frag_smiles_mass.add(String.format("%.5f", sn9_mass) + " " + sn9_mass);
////					frag_smiles_mass.add(String.format("%.5f", sn10_mass) + " " + c10_smiles);
////					frag_smiles_mass.add(String.format("%.5f", sn11_mass) + " " + c11_smiles);
////					frag_smiles_mass.add(String.format("%.5f", sn12_mass) + " " + c12_smiles);
////	//					
////	//				System.out.println("\n<---- ----->\n");
////	//				System.out.println(sn11_mass + "\t" + c11_smiles);
////	//				System.out.println(sn12_mass + "\t" + c12_smiles);
////	//				System.out.println(sn9_mass + "\t" + sn9_mass);
////	//				System.out.println(sn10_mass + "\t" + c10_smiles);
////	//				
////	 				if((sn9_mass == sn10_mass) && (sn9_mass == sn11_mass) && (sn9_mass ==sn12_mass)){
////	//					// e.g. 14:0/14:0/14:0/14:0
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)",((sn1_mass)/2.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)",(sn3_mass)/2.0, 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 5 6 7 (1 1 1 1)",(sn5_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 9 10 11 (1 1 1 1)",(sn9_mass), 999.0/999.0));
////		
////					}
////	 				else if((sn11_mass == sn12_mass) && (sn11_mass == sn9_mass) && (c11_smiles != c10_smiles)){
////						// e.g. 14:0/14:0/14:0/16:0
////				
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)",((sn1_mass)/2.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0));
////	
////						peaks.add(String.format("%.5f %5f 7 (1)",(sn6_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 5 6 (1 1 1 )",(sn5_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 11 (1)",(sn10_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 9 10 (1 1 1)",(sn9_mass), 999.0/999.0));
////			
////					} 
////	 			//  sn1 FA(sn11_mass); sn2 FA (sn12_mass); sn3 FA (sn9_mass); sn4 FA (sn10_mass)
////					else if((sn11_mass == sn12_mass) && (sn11_mass == sn10_mass) && (sn9_mass != sn10_mass)){
////						// e.g. 14:0/14:0/16:0/14:0
////	
////						peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0));				
////						peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 6 7 (1 1 1)",(sn6_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 (1)",(sn9_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 9 10 11 (1 1 1)",(sn10_mass), 999.0/999.0));
////	
////					} 
////	 				else if((sn11_mass != sn12_mass) && (sn11_mass == sn9_mass) && (sn11_mass == sn10_mass)){
////						// e.g. CL(22:6/22:5/22:6/22:6)
////					
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)",((sn1_mass)/2.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0));						
////						peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 5 7 (1 1 1)",(sn8_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 11 (1)",(sn12_mass), 999.0/999.0));		
////						peaks.add(String.format("%.5f %5f 8 9 10 (1 1 1)",(sn9_mass), 999.0/999.0));
////						
////					} 
////	 				else if((sn11_mass != sn12_mass) && (sn12_mass == sn9_mass) && (sn12_mass == sn10_mass)){
////						// e.g. 16:0/14:0/14:0/14:0
////	
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)	
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 4 5 6 (1 1 1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 10 (1)",(sn11_mass), 999.0/999.0)); // sn1 FA
////	 					peaks.add(String.format("%.5f %5f 8 9 11 (1)",(sn9_mass), 999.0/999.0)); // sn3 FA
////	
////					}				
////	 				else if((sn11_mass != sn12_mass) && (sn10_mass == sn9_mass) && (sn11_mass != sn10_mass) && (sn12_mass != sn10_mass)){
////						// e.g. CL(16:1/14:0/18:3/18:3)
////	 					System.out.println("CASE ABCC");
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 4 5 (1 1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 7 (1)",(sn8_mass), 400.0/999.0)); // [sn2+C3H6PO4-H](-) (137.00)61
////	 					peaks.add(String.format("%.5f %5f 8 9 (1 1)",(sn9_mass), 999.0/999.0)); // sn3 FA || sn4 FA
////	 					peaks.add(String.format("%.5f %5f 10 (1)",(sn11_mass), 999.0/999.0)); // sn1 FA
////	 					peaks.add(String.format("%.5f %5f 11 (1)",(sn12_mass), 999.0/999.0)); // sn2 FA
////		
////					} 
////	 				else if((sn12_mass == sn10_mass) && (sn11_mass != sn12_mass) && (sn11_mass != sn9_mass)
////							 && (sn12_mass != sn9_mass)){
////						// e.g. CL(16:0/14:0/20:2/14:0)
////	
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 5 7 (1 1)",(sn8_mass), 400.0/999.0)); // [sn2+C3H6PO4-H](-) (137.00) || [sn4+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 8 (1)",(sn9_mass), 999.0/999.0)); // sn3 FA
////	 					peaks.add(String.format("%.5f %5f 10 (1)",(sn11_mass), 999.0/999.0)); // sn1 FA
////	 					peaks.add(String.format("%.5f %5f 9 11 (1 1)",(sn12_mass), 999.0/999.0)); // sn2 FA || sn4 FA
////		
////						
////					} 
////	 				else if((sn12_mass == sn9_mass) && (sn11_mass != sn12_mass) && (sn11_mass != sn10_mass)
////					 && (sn12_mass != sn10_mass)){
////							// e.g. CL(18:0/16:1/16:1/14:0)
////	 					
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 4 7 (1)",(sn8_mass), 400.0/999.0)); // [sn2+C3H6PO4-H](-) (137.00) || [sn3+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 5 (1)",(sn6_mass), 400.0/999.0)); // [sn4+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 10 (1)",(sn11_mass), 999.0/999.0)); // sn1 FA
////	 					peaks.add(String.format("%.5f %5f 8 11 (1)",(sn12_mass), 999.0/999.0)); // sn2 FA || sn3 FA
////	 					peaks.add(String.format("%.5f %5f 9 (1)",(sn10_mass), 999.0/999.0)); // sn4 FA
////	
////					} 
////	 				else if((sn11_mass == sn12_mass) && (sn11_mass != sn9_mass) && (sn11_mass != sn10_mass) && (sn9_mass != sn10_mass) ){
////	 					// e.g. CL(14:0/14:0/16:0/18:1)
////	 					
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 5 (1)",(sn6_mass), 400.0/999.0)); // [sn4+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 6 7 (1 1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00) || [sn2+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 8 (1)",(sn9_mass), 999.0/999.0)); // sn3 FA
////	 					peaks.add(String.format("%.5f %5f 9 (1)",(sn10_mass), 999.0/999.0)); // sn4 FA 					
////	 					peaks.add(String.format("%.5f %5f 10 11 (1 1)",(sn11_mass), 999.0/999.0)); // sn1 FA || sn2 FA
////	
////	 				} 
////	 				else if((sn11_mass != sn12_mass) && (sn11_mass != sn9_mass) && (sn12_mass != sn9_mass) && (sn3_mass == sn4_mass) ){
////	 					// e.g. CL(18:0/20:4/18:1/18:1)
////	 					
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 4 5 (1 5)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00) || [sn4+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 7 (1)",(sn8_mass), 400.0/999.0)); // [sn2+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 8 9 (1 1)",(sn9_mass), 999.0/999.0)); // sn3 FA || sn4 FA
////	 					peaks.add(String.format("%.5f %5f 10 (1)",(sn11_mass), 999.0/999.0)); // sn1 FA
////	 					peaks.add(String.format("%.5f %5f 11 (1)",(sn12_mass), 999.0/999.0)); // sn2 FA				
////	
////	 				} 
////	 				else if((sn11_mass != sn12_mass) && (sn11_mass != sn9_mass) && (sn9_mass != sn12_mass) && (sn9_mass == sn10_mass)){
////	 					// e.g. CL(16:1/14:0/18:2/18:2)
////	 					
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)			
////	 					peaks.add(String.format("%.5f %5f 4 5 (1 1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00) || [sn4+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 7 (1)",(sn8_mass), 400.0/999.0)); // [sn2+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 8 9 (1 1)",(sn9_mass), 999.0/999.0)); // sn3 FA || sn4 FA
////	 					peaks.add(String.format("%.5f %5f 10 (1)",(sn11_mass), 999.0/999.0)); // sn1 FA
////	 					peaks.add(String.format("%.5f %5f 11 (1)",(sn12_mass), 999.0/999.0)); // sn2 FA					
////	
////	 				} 
////	 				else if((sn11_mass == sn12_mass) && (sn11_mass != sn9_mass) && (sn11_mass != sn10_mass) && (sn9_mass != sn10_mass)){
////	 					// e.g. CL(18:2/18:2/16:1/14:0)
////	 					
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 5 (1)",(sn6_mass), 400.0/999.0)); // [sn4+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 6 7 (1 1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00) || [sn2+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 8 (1)",(sn9_mass), 999.0/999.0)); // sn3 FA
////	 					peaks.add(String.format("%.5f %5f 9 (1)",(sn10_mass), 999.0/999.0)); // sn4 FA
////	 					peaks.add(String.format("%.5f %5f 10 11 (1 1)",(sn11_mass), 999.0/999.0)); // sn1 FA || sn2 FA
////	 					
////	 				} 
////	 				else if((sn1_mass == sn2_mass) && (sn2_mass != sn3_mass) && (sn3_mass == sn4_mass)){
////						// e.g. 16:0/14:0/14:0/16:0
////	
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 4 7 (1 1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00) || [sn2+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 5 6 (1 1)",(sn6_mass), 400.0/999.0)); // [sn4+C3H6PO4-H](-) (137.00) || [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 8 11 (1 1)",(sn9_mass), 999.0/999.0)); // sn3 FA || sn2 FA
////	 					peaks.add(String.format("%.5f %5f 9 10 (1 1)",(sn10_mass), 999.0/999.0)); // sn4 FA || sn1 FA
////		
////	 				} 
////	 				else if((sn11_mass != sn12_mass) && (sn11_mass != sn9_mass) && (sn11_mass != sn10_mass) && 
////							(sn12_mass != sn9_mass) && (sn12_mass != sn10_mass) && (sn9_mass != sn10_mass)){
////						// e.g. CL(16:1/14:0/18:1/18:3)
////	
////	 					peaks.add(String.format("%.5f %5f 0 (1)",((sn1_mass)/2.0), 500.0/999.0)); // [M-2H](2-) -sn1 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 1 (1)",(sn2_mass)/2.0, 500.0/999.0)); // [M-2H](2-) -sn3 (monolyso)
////	 					peaks.add(String.format("%.5f %5f 2 (1)",(sn3_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn1+sn2 dilyso)
////	 					peaks.add(String.format("%.5f %5f 3 (1)",(sn4_mass)/2.0, 100.0/999.0)); // [M-2H](2-) -(sn3+sn4 dilyso)
////	 					peaks.add(String.format("%.5f %5f 4 (1)",(sn5_mass), 400.0/999.0)); // [sn3+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 5 (1)",(sn6_mass), 400.0/999.0)); // [sn4+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 6 (1)",(sn7_mass), 400.0/999.0)); // [sn1+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 7 (1)",(sn8_mass), 400.0/999.0)); // [sn2+C3H6PO4-H](-) (137.00)
////	 					peaks.add(String.format("%.5f %5f 8 (1)",(sn9_mass), 999.0/999.0)); // sn3 FA
////	 					peaks.add(String.format("%.5f %5f 9 (1)",(sn10_mass), 999.0/999.0)); // sn4 FA
////	 					peaks.add(String.format("%.5f %5f 10 (1)",(sn11_mass), 999.0/999.0)); // sn1 FA
////	 					peaks.add(String.format("%.5f %5f 11 (1)",(sn12_mass), 999.0/999.0)); // sn2 FA
////					}
////	
////					s_peaks = String.format("MW: %.5f", (mass - 2.0)/2.0) + "\n" + String.format("PRECURSORMZ: %.5f",(mass - 2.0)/2.0) + "\n" + 
////							"Comment:" + String.format("Parent= %.5f", (mass - 2.0)/2.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);
////				
////					break;
////					
////				case CERAMIDE_1_PHOSPHATES:
//////					System.out.println(fragments.getAtomContainerCount());
//////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
//////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
//////					}
////	
////					// adductType.add("[M+H]+");
////					
////					//[M+H]+
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					//[M+H]-H2O (-18)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// M+H]+ (-HPO3) (-79.96633)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////					// [M+H]+ (-H3PO4) (-97.9769)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					// [M+H]+ (-H3PO4-H20)
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					// sn2 N
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass + 1.0) + " " + sg.create(fragments.getAtomContainer(5)));
////					
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 500.0/500.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass + 1.0, 250.0/500.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 1.0, 200.0/500.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + 1.0, 200.0/500.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn5_mass + 1.0, 200.0/500.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",sn6_mass + 1.0, 300.0/500.0));
////					
//////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 1.0), 200.0/999.0));
//////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass + 1.0, 50.0/999.0));
//////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 1.0, 999.0/999.0));
//////					peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass + 1.0, 500.0/999.0));
//////					peaks.add(String.format("%.5f %5f 2 (1)",sn5_mass + 1.0, 999.0/999.0));
//////					peaks.add(String.format("%.5f %5f 3 (1)",sn6_mass + 1.0, 500.0/999.0));
////					
////					s_peaks = String.format("MW: %.5f", mass + 1) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 1) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f",  mass + 1)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					break;
////					
////				case SULFATIDES:
////	
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////	
////					// adductType.add("[M-H]-");
////					
////					//[M-H]-H2O (-18)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					//[M-H]-sn2
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// [M-H]-sn2-H20
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////					// ion C6H9O8S- (241.00181)
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					// ion SO4H- (96.95956)
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////	
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass), 100.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass, 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass, 200.0/999.0));
////					peaks.add(String.format("%.5f %5f 3 (1)",sn4_mass, 999.0/999.0));
////					peaks.add(String.format("%.5f %5f 4 (1)",sn5_mass, 999.0/999.0));
////					
////					s_peaks = String.format("MW: %.5f", mass - 1) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f",  mass - 1)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					break;
////					
////					
////	//			case GANGLIOSIDES:
////	//	
////	//				double cer_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////	//				System.out.println("cer_mass: " + cer_mass);
////	//				
////	//				frag_smiles_mass.add(cer_mass + "\t" + sg.create(fragments.getAtomContainer(0)));
////	//				
////	//				peaks.add(String.format("%.5f %4d \"[M-H]-\"",(mass - 1.0), 999));
////	//				peaks.add(String.format("%.5f %4d \"[M-H]-NeuAc-H\"",(mass - 1.0 - 293.27091181494916), 500));	
////	//				peaks.add(String.format("%.5f %4d \"[M-H]-Cer-C6H10O5\"",(cer_mass -1 + 164.15672906448447 -1 ), 200)); // The two (-1) represent H atoms rermoval
////	//				peaks.add(String.format("%.5f %4d \"ion ceramide\"",(cer_mass), 200));
////	//				peaks.add(String.format("%.5f %4d \"ion C11H16NO8- (290.08759)\"",290.08759, 400));
////	//				
////	//				s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////	//						"Comment:"  + String.format("Parent= %.5f", mass - 1.0) + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////	//				metadata.add(s_peaks);					
////	//				
////	//				break;
////	//			
////				
////	//			case CHOLESTERYL_ESTERS:
////	//
////	//				System.out.println(fragments.getAtomContainerCount());
////	//				for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////	//					System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////	//				}
////	//
////	//				// adductType.add("[M-NH4]+");
////	//				peaks.add(String.format("%.5f %4d \"sterol frag. (369.352)\"",369.35213, 999));
////	//				peaks.add(String.format("%.5f %4d \"sterol frag. (175.148)\"",175.14868, 300));
////	//				peaks.add(String.format("%.5f %4d \"sterol frag. (161.133)\"",161.13303, 350));
////	//				peaks.add(String.format("%.5f %4d \"sterol frag. (147.117)\"",147.11738, 400));
////	//				peaks.add(String.format("%.5f %4d \"sterol frag. (135.117)\"",135.11738, 400));
////	//
////	//				s_peaks = String.format("MW: %.5f", mass + 18.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 18.0) + "\n" + 
////	//							"Comment:" + String.format("Parent= %.5f", mass + 18.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////	//					metadata.add(s_peaks);	
////	//				
////	//				
////	//				break;
////	//				
////	//			case CHOLESTERYL_ESTERS_LANOSTERYL_OR_ZYMOSTERYL:
////	//				
////	//				peaks.add(String.format("%.5f %4d \"sterol frag. (369.352)\"",369.35213, 999));
////	//				peaks.add(String.format("%.5f %4d \"sterol frag. C10H15  (135.117) || sterol frag. C11H15 (147.117) || sterol frag. C12H17 (161.133) || sterol frag. C13H19 (175.148) \"",0.00000, 300));
////	//
////	//				s_peaks = String.format("MW: %.5f", mass + 18.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 18.0) + "\n" + 
////	//							"Comment:" + String.format("Parent= %.5f", mass + 18.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////	//					metadata.add(s_peaks);	
////	//
////	//				break;		
////	
////				case  _1_MONOACYLGLYCEROLS:				
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////	
////					// adductType.add("[M+Li]+");
////					
////					// [M+Li]+
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// C3H8O3Li (99.06)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// C3H8O3Li-H20 (81.05)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////	
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 7.0) + " " + (sg.create(fragments.getAtomContainer(0)) + ".[Li]" ));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 7.0) + " " + (sg.create(fragments.getAtomContainer(1)) + ".[Li]" ));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 7.0) + " " + (sg.create(fragments.getAtomContainer(2)) + ".[Li]"));
////	
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 7.0), 600.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass + 7.0, 800.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 7.0, 999.0/999.0));
////					
////					s_peaks = String.format("MW: %.5f", mass + 7.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 7.0) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f",  mass + 7.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////	
////					break;			
////	
////				case  _2_MONOACYLGLYCEROLS:
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////	
////					// adductType.add("[M+Li]+");
////					
////					// [M+Li]+
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// C3H8O3Li (99.06)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// C3H8O3Li-H20 (81.05)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;	
////	
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 7.0) + " " + (sg.create(fragments.getAtomContainer(0)) + ".[Li]" ));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 7.0) + " " + (sg.create(fragments.getAtomContainer(1)) + ".[Li]" ));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 7.0) + " " + (sg.create(fragments.getAtomContainer(2)) + ".[Li]"));
////	
////					peaks.add(String.format("%.5f %5f 0 (1)",(sn1_mass + 7.0), 600.0/999.0));
////					peaks.add(String.format("%.5f %5f 1 (1)",sn2_mass + 7.0, 800.0/999.0));
////					peaks.add(String.format("%.5f %5f 2 (1)",sn3_mass + 7.0, 999.0/999.0));
////	
////					
////					s_peaks = String.format("MW: %.5f", mass + 7.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 7.0) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f",  mass + 7.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////	
////					break;
////					
////				case DIACYLGLYCEROLS:
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println(Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))* 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					// adductType.add("[M+Li]+");
////	
////	//				[M+Li]+ 600
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////	//				[M+Li]-H20 (-18) 200
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////	//				[M+Li]-sn1 400
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////	//				[M+Li]-sn1+H20 400
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////	//				[M+Li]-sn2 400
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////	//				[M+Li]-sn2+H20 400
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////	//				sn1+Li 999
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////	//				sn2+Li 999
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(0)) + ".[Li]");
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(1)) + ".[Li]");
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(2)) + ".[Li]");
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(3)) + ".[Li]");
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(4)) + ".[Li]");
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(5)) + ".[Li]");
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(6)) + ".[Li]");
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(7)) + ".[Li]");
////					
////					if(sn7_mass == sn8_mass){
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass + 7.0), 600.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass + 7.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)", (sn3_mass + 7.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 5 (1 1)", (sn5_mass + 7.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 7 (1 1)", (sn7_mass + 7.0), 999.0/999.0));
////						
////					}else{
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass + 7.0), 600.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass + 7.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass + 7.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)", (sn4_mass + 7.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass + 7.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass + 7.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass + 7.0), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass + 7.0), 999.0/999.0));
////					}
////	
////					
////					s_peaks = String.format("MW: %.5f", mass + 7.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 7.0) + "\n" + 
////							"Comment: "  + String.format("Parent= %.5f",  mass + 7.0)  + "; " + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////				break;
////				
////				case TRIACYLGLYCEROLS:
////					
////					System.out.println(fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println((Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) + 7.0)  * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M+Li]+");
////					
////					// Watch out. Here, using the partition function, sn2 and sn3 seem to be interchanged.
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass + 7.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					
////					if(sn1_mass != sn2_mass && sn1_mass != sn3_mass && sn2_mass != sn3_mass){
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass + 7.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass + 7.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass + 7.0), 500.0/999.0));					
////					}else if(sn1_mass == sn2_mass && sn1_mass != sn3_mass){
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)", (sn1_mass + 7.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass + 7.0), 500.0/999.0));
////					}else if(sn1_mass != sn2_mass && sn2_mass == sn3_mass){
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass + 7.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 2 (1 1)", (sn2_mass + 7.0), 500.0/999.0));
////						
////					} else if(sn1_mass == sn2_mass && sn1_mass == sn3_mass){
////						peaks.add(String.format("%.5f %5f 0 1 2 (1 1 1)", (sn1_mass + 7.0), 500.0/999.0));
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass + 7) + "\n" + String.format("PRECURSORMZ: %.5f",mass + 7) + "\n" + 
////							"Comment:"  + String.format("Parent= %.5f" ,  mass + 7.0) + "; "  + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////					
////					break;
////					
////				case MONOGALACTOSYLDIACYLGLYCEROLS:
////					// MGDG
////					
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println((Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) - 1.0)  * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M-H]-");
////					
////					// Watch out. Here, using the partition function, sn2 and sn3 seem to be interchanged.
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					
////					if(sn1_mass == sn2_mass){
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)", (sn1_mass - 1.0), 1.0));
////	
////					}else{
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass  - 1.0), 1.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 1.0));					
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment:" +  String.format("Parent= %.5f" ,  mass - 1.0) + "; "  + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					
////					break;
////					
////				case DIGALACTOSYLDIACYLGLYCEROLS:
////					// DGDG
////					
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println((Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) - 1.0)  * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M-H]-");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					
////					if(sn1_mass == sn2_mass){
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)", (sn1_mass - 1.0), 1.0));
////	
////					}else{
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass  - 1.0), 1.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 1.0));					
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment:" + String.format("Parent= %.5f" ,  mass - 1.0) + "; "  + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);						
////					break;
////					
////					
////				case SULFOQUINOVOSYLDIACYLGLYCEROLS:
////					// SQDG
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println((Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)) - 1.0)  * 100000)/100000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M-H]-");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					
////					if(sn1_mass == sn2_mass){
////						peaks.add(String.format("%.5f %5f 0 1 (1 1)", (sn1_mass - 1.0), 1.0));
////	
////					}else{
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass  - 1.0), 1.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 1.0));					
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment:" + String.format("Parent= %.5f" ,  mass - 1.0) + "; "  + adductType + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////					break;
////					
////				case DIACYLATED_PHOSPHATIDYLINOSITOL_MONOMANNOSIDES:
////					//  Ac2PIM1
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println( String.format("%.5f", (Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))  * 10000000)/10000000)) + "\t" + sg.create(fragments.getAtomContainer(l)) );
////					}
////					
////					// adductType.add("[M-H]-");
////					
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////					sn9_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(8)) * 100000)/100000;
////	
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass) + " " + sg.create(fragments.getAtomContainer(5)));
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass) + " " + sg.create(fragments.getAtomContainer(6)));
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass) + " " + sg.create(fragments.getAtomContainer(7)));
////					frag_smiles_mass.add(String.format("%.5f", sn9_mass) + " " + sg.create(fragments.getAtomContainer(8)));
////	
////					if(sn8_mass == sn9_mass){
////						
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 (1)", (sn3_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass), 600.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 6 (1)", (sn6_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 8 (1)", (sn8_mass), 400.0/999.0));
////	
////					}else{
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)", (sn4_mass), 300.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass), 600.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass), 400.0/999.0));	
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment:" +  String.format("Parent= %.5f" ,  mass - 1.0) + "; "  + adductType +"\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////	
////					break;
////	//				
////				case DIACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES:
////					// Ac2PIM2
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println((Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l)))  * 10000000)/10000000 + "\t" + sg.create(fragments.getAtomContainer(l)));
////					}
////					
////					// adductType.add("[M-H]-");
////		
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////	
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(5)));
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(6)));
////	
////					if(sn6_mass == sn7_mass){				
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)", (sn3_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass), 600.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 6 (1 1)", (sn6_mass), 400.0/999.0));
////					}else{
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)", (sn4_mass), 300.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass), 600.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass), 400.0/999.0));	
////					}
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment:" +  String.format("Parent= %.5f" ,  mass - 1.0) + "; "  + adductType +"\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);	
////			
////					break;
////					
////					
////				case TRIACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES:
////					// Ac3PIM2
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println( String.format("%.5f", (Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))  * 10000000)/10000000)) + "\t" + sg.create(fragments.getAtomContainer(l)) );
////					}
////					
////	
////					// adductType.add("[M-H]-");
////		
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////					sn9_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(8)) * 100000)/100000;
////					sn10_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(9)) * 100000)/100000;
////					sn11_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(10)) * 100000)/100000;
////					sn12_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(11)) * 100000)/100000;
////					sn13_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(12)) * 100000)/100000;
////					
////	//				System.out.println("sn1: " + sn11_mass);
////	//				System.out.println("sn2: " + sn12_mass);
////	//				System.out.println("sn3: " + sn13_mass);
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(5)));
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(6)));
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(7)));
////					frag_smiles_mass.add(String.format("%.5f", sn9_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(8)));
////					frag_smiles_mass.add(String.format("%.5f", sn10_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(9)));
////					frag_smiles_mass.add(String.format("%.5f", sn11_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(10)));
////					frag_smiles_mass.add(String.format("%.5f", sn12_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(11)));
////					frag_smiles_mass.add(String.format("%.5f", sn13_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(12)));
////	
////					
////					
////					if((sn11_mass == sn12_mass) && (sn11_mass == sn13_mass)){	
////						// Ac3PIM2(15:0/15:0/15:0)
////						System.out.println("CASE 1");
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 50.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 5 (1 1 1)", (sn3_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 6 7 (1 1 1)", (sn6_mass - 1.0), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 10 11 12 (1 1 1)", (sn11_mass - 1.0), 400.0/999.0));
////					}
////					
////					else if ((sn11_mass != sn12_mass) && (sn11_mass == sn13_mass)){
////						// Ac3PIM2(15:0/14:0/15:0)
////						System.out.println("CASE 2");
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 50.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)", (sn3_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 6 (1 1)", (sn5_mass - 1.0), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 11 12 (1 1)", (sn12_mass - 1.0), 400.0/999.0));
////						
////					}
////					else if((sn11_mass == sn12_mass) && (sn13_mass != sn12_mass)){
////						// Ac3PIM2(15:0/14:0/14:0)
////						System.out.println("CASE 3");
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 50.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 5 (1 1)", (sn4_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 7 (1 1)", (sn7_mass - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass - 1.0), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 12 (1)", (sn13_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 10 11 (1 1)", (sn11_mass - 1.0), 400.0/999.0));
////		
////					}
////					else if((sn11_mass != sn12_mass) && (sn13_mass == sn12_mass)){
////						// Ac3PIM2(14:0/14:0/16:0)
////						System.out.println("CASE 4");
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 50.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)", (sn3_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 6 (1 1)", (sn5_mass - 1.0), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 11 12 (1 1)", (sn12_mass - 1.0), 400.0/999.0));
////	
////					}
////					else if((sn11_mass != sn12_mass) && (sn11_mass != sn13_mass) && (sn12_mass != sn13_mass)){
////						// Ac3PIM2(14:0/15:0/16:0)
////						System.out.println("CASE 5");
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 50.0/999.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 3 (1)", (sn4_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass - 1.0), 999.0/999.0));
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 500.0/999.0));
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 200.0/999.0));
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 11 (1)", (sn12_mass - 1.0), 400.0/999.0));
////						peaks.add(String.format("%.5f %5f 12 (1)", (sn13_mass - 1.0), 400.0/999.0));
////					
////					}
////	
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment:" +  String.format("Parent= %.5f" ,  mass - 1.0) + "; "  + adductType +"\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////	//				
////					break;				
////	//				
////	//				
////				case TETRACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES:
////					// Ac4PIM2
////					
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println( String.format("%.5f", (Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))  * 10000000)/10000000)) + "\t" + sg.create(fragments.getAtomContainer(l)) );
////					}
////	
////					// adductType.add("[M-H]-");
////		
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////					sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////					sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////					sn9_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(8)) * 100000)/100000;
////					sn10_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(9)) * 100000)/100000;
////					sn11_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(10)) * 100000)/100000;
////					
////	//				System.out.println("sn1: " + sn11_mass);
////	//				System.out.println("sn2: " + sn12_mass);
////	//				System.out.println("sn3: " + sn13_mass);
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(5)));
////					frag_smiles_mass.add(String.format("%.5f", sn7_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(6)));
////					frag_smiles_mass.add(String.format("%.5f", sn8_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(7)));
////					frag_smiles_mass.add(String.format("%.5f", sn9_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(8)));
////					frag_smiles_mass.add(String.format("%.5f", sn10_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(9)));
////					frag_smiles_mass.add(String.format("%.5f", sn11_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(10)));
////	
////					if ((sn2_mass == sn3_mass) && (sn2_mass == sn4_mass) && (sn2_mass == sn5_mass)){
////						// Ac4PIM2(16:0/16:0/16:0/16:0)
////						System.out.println("CASE 1");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1 || [M-H]-sn2 || [M-H]-sn3 || [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 1 2 3 4 (1 1 1 1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn1-sn2 || [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 5 7 (1 1)", (sn6_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn1-sn2-C3H5O || [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 6 8 (1 1)", (sn7_mass - 1.0), 300.0/999.0));
////						// sn1 FA || sn2 FA
////						peaks.add(String.format("%.5f %5f 9 10 (1 1)", (sn11_mass - 1.0), 400.0/999.0));
////						
////					}
////					else if ((sn3_mass == sn4_mass) && (sn4_mass == sn5_mass) && (sn2_mass != sn3_mass)){
////						// Ac4PIM2(14:0/15:0/15:0/15:0)
////						System.out.println("CASE 2");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn2 || [M-H]-sn3 || [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 2 3 4 (1 1 1)", (sn3_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 100.0/999.0));					
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 300.0/999.0));
////						// sn2 FA
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						// sn1 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////						
////					}
////					else if ((sn2_mass == sn3_mass) && (sn2_mass == sn4_mass) && (sn2_mass != sn5_mass)){
////						// Ac4PIM2(16:0/16:0/16:0/14:0)
////						System.out.println("CASE 3");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn1 || [M-H]-sn2 || [M-H]-sn3
////						peaks.add(String.format("%.5f %5f 1 2 3 (1 1 1)", (sn2_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 100.0/999.0));	
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 300.0/999.0));
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// sn1 FA || sn2 FA
////						peaks.add(String.format("%.5f %5f 9 10 (1 1)", (sn11_mass - 1.0), 400.0/999.0));
////						
////					}
////					else if((sn2_mass == sn4_mass) && (sn3_mass == sn5_mass) && (sn2_mass != sn3_mass)){	
////						// Ac4PIM2(14:0/15:0/14:0/15:0)
////						System.out.println("CASE 4");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1
////						peaks.add(String.format("%.5f %5f 1 3 (1 1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn2
////						peaks.add(String.format("%.5f %5f 2 4 (1 1)", (sn3_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 7 (1 1)", (sn6_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 8 (1 1)", (sn7_mass - 1.0), 300.0/999.0));
////						// sn2 FA
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						// sn1 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////					}
////					
////					else if ((sn3_mass == sn4_mass) && (sn4_mass == sn5_mass) && (sn2_mass != sn3_mass)){
////						// Ac4PIM2(14:0/16:0/16:0/14:0)
////						System.out.println("CASE 5");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1 || [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 1 4 (1 1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn2 || [M-H]-sn3
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)", (sn3_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn1-sn2 || [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 5 7 (1 1)", (sn6_mass - 1.0), 100.0/999.0));					
////						// [M-H]-sn1-sn2-C3H5O || [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 6 8 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// sn2 FA
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						// sn1 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////						
////					}				
////	
////					else if((sn3_mass == sn4_mass) && (sn2_mass != sn3_mass) && (sn4_mass != sn5_mass) && (sn2_mass != sn5_mass)){
////						// Ac4PIM2(14:0/15:0/15:0/16:0)
////						System.out.println("CASE 6");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn2 || [M-H]-sn3
////						peaks.add(String.format("%.5f %5f 2 3 (1 1)", (sn3_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 100.0/999.0));					
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 300.0/999.0));
////						// sn2 FA
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						// sn1 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////		
////					}
////					else if((sn2_mass != sn3_mass) && (sn2_mass != sn4_mass) && (sn3_mass != sn4_mass) && (sn4_mass == sn5_mass)){
////						// Ac4PIM2(14:0/15:0/16:0/16:0)
////						System.out.println("CASE 7");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn2
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn3 || [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 3 4 (1 1)", (sn4_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 100.0/999.0));					
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 300.0/999.0));
////						// sn2 FA
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						// sn1 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////		
////					}
////					else if((sn2_mass == sn3_mass) && (sn2_mass != sn4_mass) && (sn2_mass != sn5_mass) && (sn4_mass != sn5_mass)){
////						// Ac4PIM2(16:0/16:0/17:0/14:0)
////						System.out.println("CASE 8");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn1 || [M-H]-sn2
////						peaks.add(String.format("%.5f %5f 1 2 (1 1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn3
////						peaks.add(String.format("%.5f %5f 3 (1)", (sn4_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 100.0/999.0));	
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 300.0/999.0));
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// sn1 FA || sn2 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////		
////					}
////					else if((sn2_mass != sn3_mass) && (sn2_mass == sn4_mass) && (sn2_mass == sn5_mass)){
////						// Ac4PIM2(14:0/16:0/14:0/14:0)
////						System.out.println("CASE 9");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1 || [M-H]-sn3 || [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 1 3 4 (1 1 1)", (sn2_mass - 1.0), 999.0/999.0));	
////						// [M-H]-sn2
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 100.0/999.0));
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));	
////						// [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 300.0/999.0));
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// sn2 FA
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						// sn1 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////						
////					}
////	
////					else if((sn2_mass != sn3_mass) && (sn2_mass != sn4_mass) && (sn2_mass != sn5_mass) && (sn3_mass != sn4_mass) && (sn3_mass != sn5_mass) && (sn4_mass != sn5_mass) ){
////						// Ac4PIM2(14:0/15:0/17:0/18:2)
////						System.out.println("CASE 10");
////						// [M-H]-
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass - 1.0), 200.0/999.0));
////						// [M-H]-sn1
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass - 1.0), 999.0/999.0));
////						// [M-H]-sn2
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn3
////						peaks.add(String.format("%.5f %5f 3 (1)", (sn4_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn4
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass - 1.0), 500.0/999.0));
////						// [M-H]-sn1-sn2
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass - 1.0), 100.0/999.0));				
////						// [M-H]-sn1-sn2-C3H5O
////						peaks.add(String.format("%.5f %5f 6 (1)", (sn7_mass - 1.0), 300.0/999.0));
////						// [M-H]-sn3-sn4
////						peaks.add(String.format("%.5f %5f 7 (1)", (sn8_mass - 1.0), 100.0/999.0));	
////						// [M-H]-sn3-sn4-C3H5O
////						peaks.add(String.format("%.5f %5f 8 (1)", (sn9_mass - 1.0), 300.0/999.0));
////						// sn2 FA
////						peaks.add(String.format("%.5f %5f 9 (1)", (sn10_mass - 1.0), 400.0/999.0));
////						// sn1 FA
////						peaks.add(String.format("%.5f %5f 10 (1)", (sn11_mass - 1.0), 400.0/999.0));
////					
////					}
////	
////	
////					s_peaks = String.format("MW: %.5f", mass - 1.0) + "\n" + String.format("PRECURSORMZ: %.5f",mass - 1.0) + "\n" + 
////							"Comment:" + "\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);				
////					break;
////	//				
////	////				
////	////			case DIPHOSPHORYLATED_HEXAACYL_LIPID_A:
////	////				
////	////				break;
////	//		
////				case ACYL_CARNITINES:
////					
////					System.out.println("Nr. of fragments: " + fragments.getAtomContainerCount());
////					for(int l = 0; l < fragments.getAtomContainerCount(); l++){
////						System.out.println( String.format("%.5f", (Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(l))  * 10000000)/10000000)) + "\t" + sg.create(fragments.getAtomContainer(l)) );
////					}
////					
////					// adductType.add("[M]+");
////					
////					// C3H10N fragment	(m/z 60)
////					sn1_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(0)) * 100000)/100000;
////					// [M+H] - FA - N(CH3)3 - H2O  (m/z 85)
////					sn2_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(1)) * 100000)/100000;
////					// [M+H] - FA (McLafferty rearrangement: m/z 144)
////					sn3_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(2)) * 100000)/100000;
////					// [M+H] - carnitine - H20
////					sn4_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(3)) * 100000)/100000;
////					// [M+H] - carnitine  (-161)
////					sn5_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(4)) * 100000)/100000;
////					// [M+H] - N(CH3)3 -H20
////					sn6_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(5)) * 100000)/100000;
////	//				sn7_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(6)) * 100000)/100000;
////	//				sn8_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(7)) * 100000)/100000;
////	//				sn9_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(8)) * 100000)/100000;
////	//				sn10_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(9)) * 100000)/100000;
////	//				sn11_mass = Math.floor(acm.getNaturalExactMass(fragments.getAtomContainer(10)) * 100000)/100000;
////	
////					
////					frag_smiles_mass.add(String.format("%.5f", sn1_mass) + " " + sg.create(fragments.getAtomContainer(0)));
////					frag_smiles_mass.add(String.format("%.5f", sn2_mass) + " " + sg.create(fragments.getAtomContainer(1)));
////					frag_smiles_mass.add(String.format("%.5f", sn3_mass) + " " + sg.create(fragments.getAtomContainer(2)));
////					frag_smiles_mass.add(String.format("%.5f", sn4_mass) + " " + sg.create(fragments.getAtomContainer(3)));
////					frag_smiles_mass.add(String.format("%.5f", sn5_mass) + " " + sg.create(fragments.getAtomContainer(4)));
////					frag_smiles_mass.add(String.format("%.5f", sn6_mass) + " " + sg.create(fragments.getAtomContainer(5)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn7_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(6)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn8_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(7)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn9_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(8)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn10_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(9)));
////	//				frag_smiles_mass.add(String.format("%.5f", sn11_mass - 1.0) + " " + sg.create(fragments.getAtomContainer(10)));
////					
////					SmartsPattern smp = SmartsPattern.create("[#8;A;X2,X1-][#6;X3](-[#6,#1])=[O;X1]", bldr);
////					int nrOfCarboxylGroups = StructureExplorer.fragmentCount(smp, molecule);
////					System.out.println("nrOfCarboxylGroups: " + nrOfCarboxylGroups);
////					
////	//				if(nrOfCarboxylGroups == 2){
////						
////						peaks.add(String.format("%.5f %5f 0 (1)", (sn1_mass), 7.0/100.0));
////						peaks.add(String.format("%.5f %5f 1 (1)", (sn2_mass), 100.0/100.0));
////						peaks.add(String.format("%.5f %5f 2 (1)", (sn3_mass), 15.0/100.0));
////						peaks.add(String.format("%.5f %5f 3 (1)", (sn4_mass), 10.0/100.0));
////						peaks.add(String.format("%.5f %5f 4 (1)", (sn5_mass), 10.0/100.0));
////						peaks.add(String.format("%.5f %5f 5 (1)", (sn6_mass), 10.0/100.0));
////	//				}
////	
////					s_peaks = String.format("MW: %.5f", mass) + "\n" + String.format("PRECURSORMZ: %.5f",mass) + "\n" + 
////							"Comment:" +  String.format("Parent= %.5f" ,  mass) + "; "  + adductType +"\n" + String.format("Num Peaks: %2d",peaks.size());
////					metadata.add(s_peaks);
////					break;
////					
////	//			case NIL:
////	//				
////	//				break;
////					
////			}
////		}
//			
//		results.put("molName", molName);
//		results.put("peaks_list", peaks);
//		results.put("metadata", metadata);
//		results.put("fragments", frag_smiles_mass);
//		results.put("adductType", sAdductType);
//		
//		return results;
//	}
	
//	public void saveSingleCfmidLikeMSPeakListFromSDF(String sourceFilename, IChemObjectBuilder bldr) throws Exception{
//		saveSingleCfmidLikeMSPeakListFromSDF(sourceFilename, bldr, "[M+H]+");
//	}
//	
//	
//	public void saveSingleCfmidLikeMSPeakListFromSDF(String sourceFilename, IChemObjectBuilder bldr, String adductType) throws Exception{
//		
//		Path path = Paths.get(sourceFilename);
//		
//		String directory = path.getParent().toAbsolutePath().toString();
//		
//		FileWriter fwerr = new FileWriter(directory + "/missing.log");
//		BufferedWriter bwerr = new BufferedWriter(fwerr);		
//		
//		System.out.println("PATH: " + directory);
//		IteratingSDFReader sdfReader = new IteratingSDFReader(new FileReader(sourceFilename), bldr);
//		int mol_nr = 0;
//		int fmol = 0;
//		while (sdfReader.hasNext()) {
//			IAtomContainer molecule = sdfReader.next();
//			LinkedHashMap<String,ArrayList<String>> peaksResults = generateCfmidLikeMSPeakList(molecule, adductType);
//			ArrayList<String> peaks = peaksResults.get("peaks_list");
//			ArrayList<String> frags = peaksResults.get("fragments");
//			
//	
//			mol_nr++;
//			String name = peaksResults.get("molName").get(0);
//			if(name == null || name.trim().length()==0){
//				name = "molecule_" + String.valueOf(mol_nr);
//				
//			}
//
//			if(peaks.size() != 0){
//				fmol ++;
//				FileWriter fw = new FileWriter(directory + "/" + name.replace("/", "_") + ".log");
//				BufferedWriter bw = new BufferedWriter(fw);
//				
//				bw.write("energy2");
//				bw.newLine();
//				for(int p = 0; p < peaks.size();p++){
//					bw.write(peaks.get(p));
//					bw.newLine();
//				}
//			
//				bw.newLine();
//				
//				// saving the fragment list
//				for(int k = 0; k < frags.size(); k++){
//					bw.write(k +" " + frags.get(k));
//					bw.newLine();
//				}
//				bw.newLine();
//				bw.write(peaksResults.get("adductType").get(0));
//				bw.newLine();
//				
//				bw.close();	
//
//				
//			} else{
//				bwerr.write(name);
//				bwerr.newLine();
//
//			}
//			
//		}
//		sdfReader.close();
//		fwerr.close();
//		
//		System.out.println("\nTotal nr. of molecules: " + String.valueOf(mol_nr));
//		System.out.println("Total nr. of fragmented molecules: " + String.valueOf(fmol));
//	}
//	
//	
//	
	public void saveSingleCfmidLikeMSPeakList(IAtomContainer molecule, IChemObjectBuilder bldr, String adductType, String outputname) throws Exception{

		FileWriter fwerr = new FileWriter("data/missing.log");
		BufferedWriter bwerr = new BufferedWriter(fwerr);		

		int fmol = 0;

//			LinkedHashMap<String,ArrayList<String>> peaksResults = generateCfmidLikeMSPeakList(molecule, adductType);
			LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> peaksResults = generateCfmidLikeMSPeakList(molecule, adductType);
//			ArrayList<String> peaks = peaksResults.get("peaks_list");
//			ArrayList<String> frags = peaksResults.get("fragments");

//			String name = peaksResults.get("molName").get(0);
//			if(name == null || name.trim().length()==0){
//				name = "molecule_" + String.valueOf(mol_nr);
//				
//			}
			
			
			if(peaksResults != null) {
				
				FileWriter fw = new FileWriter(outputname);
//				FileWriter fw = new FileWriter("data/" + outputname.replace("/", "_") + ".log");
				BufferedWriter bw = new BufferedWriter(fw);
				int level = 0;
				int fragIndex = 0;
				
				for(Map.Entry<Integer, LinkedHashMap<String, ArrayList<String>>> peaks : peaksResults.entrySet()){
					System.out.println("Energy: " + peaks.getKey());
					System.out.println("peaks: " + peaks.getValue());
					System.out.println("Nr. of peaks: " + peaks.getValue().get("peaks_list").size());
					bw.write("energy" + level);
					bw.newLine();
					level++;
					for(String s : peaks.getValue().get("peaks_list")){
						System.out.println();
						bw.write(s);
						bw.newLine();
					}
				}
				
				bw.newLine();
				for( String fg : peaksResults.get(10).get("fragments")){
					bw.write( fragIndex + " " + fg);
					bw.newLine();
					fragIndex++;
				}

				bw.newLine();
				bw.write(adductType);
				bw.newLine();
				
				bw.close();	
				fw.close();
			}
			

		
			
			
//			if(peaks.size() != 0){
//
//				
////				FileWriter fw = new FileWriter("data/" + outputname.replace("/", "_") + ".log");
//				
////				System.out.println("Path : " + file.getAbsolutePath());
//				
//				
//
//			
//				bw.newLine();
//				
//				// saving the fragment list
//				for(int k = 0; k < frags.size(); k++){
//					bw.write(k +" " + frags.get(k));
//					bw.newLine();
//				}
//				bw.newLine();
//				bw.write(peaksResults.get("adductType").get(0));
//				bw.newLine();
//				
//				bw.close();	
//				fw.close();
//				fmol ++;
//			} else{
//				bwerr.write(name);
//				bwerr.newLine();
//
//			}

		fwerr.close();

//		System.out.println("Total nr. of fragmented molecules: " + String.valueOf(fmol));
	}
//
//	
//    public static void main( String[] args ) throws Exception
//    {
////    	AtomContainerManipulator acm = new AtomContainerManipulator();
////    	SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
////    	IAtomContainer li = sp.parseSmiles("[Li]");
////    	System.out.println(Math.floor(acm.getNaturalExactMass(li) * 100000)/100000);
////    	System.out.println(acm.getNaturalExactMass(li));
//    	
//    	
//    	
////    	String molSmiles 	= args[0];
////    	String adductType 	= args[1];
////    	String outputName	= args[2];
//    	
////    	String molSmilesPC_16_0-18_1 	= "C[N+](C)(C)CCOP(O)(=O)OCC(CO(C(=O)CCCCCCCCCCCCCCC))O(C(=O)CCCCCCCC=CCCCCCCCC)";
////    	String molSmilesPI_18_1_9Z = "O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(=O)OC[C@@H](CO(C(=O)CCCCCCCC=CCCCCCCCC))O(C(=O)CCCCCCCC=CCCCCCCCC))[C@H](O)[C@@H]1O";    	
////    	String molSmilesLPC_14_0 = "[H][C@@](O)(COC(=O)CCCCCCCCCCCCC)COP([O-])(=O)OCC[N+](C)(C)C";
//    	
////    	String molSmilesPlamalogen_18_0_18_1_9Z.log = "CCCCCCCCCCCCCCCCC=COC[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCCCCCCCC";
//    	
//    	String molSmiles = "CCCCCCCCCCCCCC=C[C@@H](O)[C@H](COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCC=CCCCCCCCC";
//    	String adductType 	= "M+";
//    	String outputName	= "data/SM18_1_9Z.log";
//    	
//    	Fragmenter fr = new Fragmenter();
//
//        SmilesParser sp =  new SmilesParser(SilentChemObjectBuilder.getInstance());
//        SmilesGenerator sg = new SmilesGenerator().unique();
//        StructureExplorer se = new StructureExplorer();
//        IAtomContainer molecule = sp.parseSmiles(molSmiles);
//        
//        System.out.println("Molecule standardized: " + sg.create(se.standardizeMolecule(molecule)) + "\n\n");
//        
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);        
//        IChemObjectBuilder bldr = DefaultChemObjectBuilder.getInstance(); 
//        fr.saveSingleCfmidLikeMSPeakList(molecule, bldr, adductType, outputName);        
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/test3.sdf",bldr,"[M+H]+");
//
//    }
	
	
}

