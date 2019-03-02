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
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

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
		// DG(16:0/18:1) - "[H][C@](CO)(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCC=CCCCCCCCC"
		// MG(18:2/0:0/0:0) - "CCCCCC=CCC=CCCCCCCCC(=O)OCC(O)CO"
		// MG(16:0/0:0/0:0) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO"
		// TG(16:0/16:0/16:0) - "CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
		// TG(20:0/20:1/20:0) - "[H]C([H])(OC(=O)CCCCCCCCCCCCCCCCCCC)C([H])(OC(=O)CCCCCCCCC\C=C/CCCCCCCC)C([H])([H])OC(=O)CCCCCCCCCCCCCCCCCCC"
		// TG(22:6/22:6/22:6) - "CCC=CCC=CCC=CCC=CCC=CCC=CCCC(=O)OCC(COC(=O)CCC=CCC=CCC=CCC=CCC=CCC=CCC)OC(=O)CCC=CCC=CCC=CCC=CCC=CCC=CCC"		
		// PA(16:0/18:1(9Z)) - "OP(O)(=O)OCC(CO(C(=O)CCCCCCCCCCCCCCC))O(C(=O)CCCCCCCC=CCCCCCCCC)"
		// PE(16:0/16:0) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OCCN)OC(=O)CCCCCCCCCCCCCCC"
		// PE(16:0/18:1(9Z)) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OCCN)OC(=O)CCCCCCCC=CCCCCCCCC"
		// PE(17:0/14:1(9Z)) - "NCCOP(O)(=O)OC[C@@H](CO(C(=O)CCCCCCCCCCCCCCCC))O(C(=O)CCCCCCCC=CCCCC)"
		// PE(16:0/20:4(5Z,8Z,11Z,14Z)) - "CCCCCC=CCC=CCC=CCC=CCCCC(=O)O[C@@H](COP(O)(=O)OCCN)COC(=O)CCCCCCCCCCCCCCC"
		// PE(17:0/20:4(5Z,8Z,11Z,14Z)) - "CCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCCN)OC(=O)CCCC=CCC=CCC=CCC=CCCCCC"
		// PE(18:0/20:4(5Z,8Z,11Z,14Z)) - "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCCN)OC(=O)CCCC=CCC=CCC=CCC=CCCCCC"		
		// PC(16:0/18:1(9Z)) - "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCCCCCCCC"
		// PC(18:0/20:4) - "C[N+](C)(C)CCOP([O-])(=O)OC[C@@H](COC(=O)CCCCCCCCCCCCCCCCC)OC(=O)CCCC=CCC=CCC=CCC=CCCCCC"
		// PC(16:0/0:0) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCC[N+](C)(C)C"
		// PC(18:1(9Z):0) - "CCCCCCCCC=CCCCCCCCC(=O)OCC(O)COP(O)(=O)OCC[N+](C)(C)C"
		// PA(18:1(9Z):0) - "CCCCCCCC/C=C\CCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O"
		 // PA(12:0/12:0)- "CCCCCCCCCCCC(=O)OCC(COP(=O)(O)O)OC(=O)CCCCCCCCCCC"
		// PS(16:0/18:1(9Z)) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OC[C@@H](N)C(O)=O)OC(=O)CCCCCCCC=CCCCCCCCC"
		// PS(16:0/20:4(5Z,8Z,11Z,14Z)) - "CCCCCC=CCC=CCC=CCC=CCCCC(=O)O[C@@H](COC(=O)CCCCCCCCCCCCCCC)COP(O)(=O)OC[C@@H](N)C(O)=O"
		// Cer(18:1/12:0) - "CCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)/C=C/CCCCCCCCCCCCC"
		// Cer(18:1/24:1) - "CCCCCCCCCCCCCC=C[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCC=CCCCCCCCCCCC"
		// Cer(18:1(9Z)/16:0) - "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCC)=O)CO"
		// Cer(18:1(9Z)/OH-18:0) - "CCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\C=C\CCCCCCCCCCCCC"
		// SM(18:1(9Z):16:0) - "CCCCCCCCCCCCCCCC(=O)N[C@H](COP(O)(=O)OCC[N+](C)(C)C)[C@@H](O)C=CCCCCCCCCCCCCC"
		// CL(1'-[18:1(9Z)/18:1(9Z)],3'-[18:1(9Z)/18:1(9Z)]) - "[H]C(O)(COP(O)(=O)OC[C@@]([H])(COC(=O)CCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCC=CCCCCCCCC)COP(O)(=O)OC[C@@]([H])(COC(=O)CCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCC=CCCCCCCCC"
		// Carnitine 16 - "CCCCCCCCCCCCCCCC(=O)OC(CC(O)=O)C[N+](C)(C)C"
		// PG(16:0/18:1(9Z)) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OCC(O)CO)OC(=O)CCCCCCCC=CCCCCCCCC"
		// PG(16:0) - "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@@H](O)CO"
		// PC(O-16:0/18:1(9Z)) - "C[N+](C)(C)CCOP(O)(=O)OC[C@@H](COCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCC=CCCCCCCCC"
		// PC(O-16:0/20:3(8Z,11Z,14Z)) - "C[N+](C)(C)CCOP([O-])(=O)OC[C@@H](COCCCCCCCCCCCCCCCC)OC(=O)CCCCCCC=CCC=CCC=CCCCCC"
		// PC(O-18:0/2:0) - "CCCCCCCCCCCCCCCCCCOC[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(C)=O"
		// PC(P-18:0/20:4(5Z,8Z,11Z,14Z)) - "C[N+](C)(C)CCOP(O)(=O)OC[C@@H](COC=CCCCCCCCCCCCCCCCC)OC(=O)CCCC=CCC=CCC=CCC=CCCCCC"
		// PC(P-18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)) - "CCCCCCCCCCCCCCCCC=COC[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCC=CCC=CCC=CCC=CCC=CCC=CCC"
		// PC(P-18:0/18:1(9Z)) - "C[N+](C)(C)CCOP([O-])(=O)OC[C@@H](COC=CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCC=CCCCCCCCC"
		// PC(O-18:0/0:0) - "CCCCCCCCCCCCCCCCCCOC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C"
		// PC(P-16:0/0:0) - "CCCCCCCCCCCCCCC=COC[C@@H](O)COP(O)(=O)OCC[N+](C)(C)C"
//		 PI(17:0/14:1(9Z)) - "CCCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)O[C@H]1C(O)[C@@H](O)C(O)C(O)C1O)OC(=O)CCCCCCCC=CCCCC"
		// Oleoyl_sn_LPC - "CCCCCCCCC=CCCCCCCCC(=O)OCC(O)COP(O)(=O)OCC[N+](C)(C)C"
		// PC(P-16:0/20:4) - "C(CCCCCC)=CCC=CCC=CCC=CCCC(OC(COP(O)(=O)OCCN)COC=CCCCCCCCCCCCCCC)=O"
		// PC(16:0/22:6) - "OP(=O)(OCCN)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCC=CCC=CCC=CCC=CCC=CCC=CCC"
		// PA(14:0/0:0)- "CCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)O)O"
		// N-Palmitoyl-D-erythro-Sphingosine - "CCCCCCCCCCCCCCCC(=O)NC(CO)C(C=CCCCCCCCCCCCCC)O"
		// Palmitoylcarnitine - ""
		// PC(18:0/18:2 - "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC=CCC=CCCCCCC"
		// PC(6:0/16:0) - "C[N+](C)(C)CCOP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
		// PC(6:0/16:0) - "C[N+](C)(C)CCOP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCC)OC(=O)CCCCCCCCCCC"
		// PC(16:0/16:0) - "C[N+](C)(C)CCOP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
		// PS(18:0/20:4) - "CCCCCC=CCC=CCC=CCC=CCCCC(=O)O[C@@H](COC(=O)CCCCCCCCCCCCCCCCC)COP(O)(=O)OC[C@@H](N)C(O)=O"
		// DG(15:0/0:0/15:0) (d5) - "[2H]C(O)(C([2H])([2H])OC(=O)CCCCCCCCCCCCCC)C([2H])([2H])OC(=O)CCCCCCCCCCCCCC"
		// DG(16:0/16:0/0:0)- "CCCCCCCCCCCCCCCC(=O)O[C@@H](CO)COC(=O)CCCCCCCCCCCCCCC"
		// PG(14:0/14:0) - "CCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OC[C@H](O)CO)OC(=O)CCCCCCCCCCCCC"
		// PG(18:0/18:0) - "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OC[C@H](O)CO)OC(=O)CCCCCCCCCCCCCCCCC"
		// PG(18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)) - "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OC[C@H](O)CO)OC(=O)CCC=CCC=CCC=CCC=CCC=CCC=CCC"
		// SM(18:1/16:0) - "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCCCC)=O)COP(OCC[N+](C)(C)C)([O-])=O"
		
		
    	String molSmiles = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCCCCCCCC";
    	String adductType 	= "[M+H]+";
    	String outputName	= "data/PC(18-0_20-4)";
    	
    	Fragmenter fr = new Fragmenter();

        
        StructureExplorer se = new StructureExplorer();
        IAtomContainer molecule = fr.sParser.parseSmiles(molSmiles.replace("[O-]", "O"));
        
//       MolecularFormulaManipulator mfm = new MolecularFormulaManipulator();
//       IMolecularFormula m = mfm.getMolecularFormula(molecule);
//       System.out.println("Natural excat mass: " + AtomContainerManipulator.getNaturalExactMass(molecule));
//       System.out.println( "Major isotopse mass: " + MolecularFormulaManipulator.getMajorIsotopeMass(m));
//       System.out.println( "Total exact mass: " +  AtomContainerManipulator.getTotalExactMass(molecule));
        
        
//        LinkedHashMap<String, IAtomContainer> fragments = fr.fragmentMolecule(molecule, ClassName._1_MONOACYLGLYCEROLS, adductType);
//
        
//        for(Map.Entry<String, IAtomContainer> fragment : fragments.entrySet()){
//        	System.out.println(fragment.getKey() + " : " + fr.sGen.create(fragment.getValue()));
//        }        
//        LinkedHashMap<String, Double> masses = fr.computeFragmentMasses(fragments);
//        System.out.println(masses);

        long t1 = System.nanoTime();        
        fr.saveSingleCfmidLikeMSPeakList(molecule, fr.bldr, adductType, outputName);
        long t2 = System.nanoTime();
        System.out.println(AtomContainerManipulator.getTotalFormalCharge(molecule));
        System.out.println("Execution time: " + ((t2 - t1) * 1e-6) + " milliseconds");
	}
	
	
	public LinkedHashMap<String, IAtomContainer>  fragmentMolecule(IAtomContainer molecule, StructuralClass.ClassName type) throws Exception{
		return fragmentMolecule(molecule, type, "[M+H]+");
	}
	
	public boolean isConditionValidForChemClass(StructuralClass.ClassName type,
			String adductType){
		boolean valid = false;
		
		if(FPLists.classSpecificFragmentationPatterns.get(type).containsKey(adductType)){
			valid = true;
		}
		
		return valid;
	}
	
	public LinkedHashMap<String, IAtomContainer> fragmentMolecule(IAtomContainer molecule, StructuralClass.ClassName type,
			String adductType) throws Exception{	
		LinkedHashMap<String, IAtomContainer> fragments = new LinkedHashMap<String, IAtomContainer>();

//		System.out.println(FPLists.classSpecificFragmentationPatterns);
//		System.out.println(FPLists.classSpecificFragmentationPatterns.get(type).get(adductType));
		
		
		
//		System.out.println("KEY SET");
//		System.out.println(FPLists.classSpecificFragmentationPatterns.get(type).get(adductType));
//		System.out.println(FPLists.classSpecificFragmentationPatterns.get(type).get(adductType)
//				.patterns.keySet());
				
		if(type != StructuralClass.ClassName.NIL && FPLists.classSpecificFragmentationPatterns.get(type).containsKey(adductType)){
			new SmilesGenerator();
			//			IAtomContainer stMol = this.sExplorer.standardizeMolecule(molecule);
			SmilesGenerator sg = SmilesGenerator.isomeric();
			
			for(Map.Entry<String, String[]> pattern : FPLists.classSpecificFragmentationPatterns.get(type).get(adductType)
					.patterns.entrySet()){
				SMIRKSReaction sReaction = this.smrkMan.parse(pattern.getValue()[1]);
//				System.out.println("pattern.getKey(): " + pattern.getKey());
//				System.out.println(sReaction == null);
				IAtomContainerSet current_set = this.smrkMan.applyTransformationWithSingleCopyForEachPos(
						molecule, null, sReaction);
				
				if(current_set != null){
					
//					 System.out.println("Nr of metabolites: " + String.valueOf(current_set.getAtomContainerCount()));	
//					 System.out.println("Nr of unique metabolites: " + String.valueOf(this.sExplorer.uniquefy(current_set).getAtomContainerCount()));					 
					 
					 IAtomContainerSet current_set_unique = StructureExplorer.uniquefy(current_set);
					 
					 for(IAtomContainer atc : current_set_unique.atomContainers()){
						 IAtomContainerSet partitions = StructureExplorer.partition(atc);

						 	for(IAtomContainer partition : partitions.atomContainers()){
						 		if(StructureExplorer.containsSmartsPattern(partition, pattern.getValue()[0])){
//								if(this.sExplorer.containsSmartsPattern(this.sExplorer.standardizeMolecule(partition), pattern.getValue()[0])){
									if(pattern.getValue().length==2){
										fragments.put(pattern.getKey(), AtomContainerManipulator.removeHydrogens(partition));
//										System.out.println(sg.create(partition) + " was added.");
										
									}else if (pattern.getValue().length==3){
										// take the SMILES of the fragment an add the adduct, whiich would be the 
										// third element in the array.
										IAtomContainer  adjustedPartition = this.sParser.parseSmiles(sg.create(AtomContainerManipulator.removeHydrogens(partition))
												+ "." + pattern.getValue()[2]);	
										fragments.put(pattern.getKey(), adjustedPartition);
//										System.out.println(sg.create(adjustedPartition) + " was added.");		
									}

									break;	
								} else{
//									System.out.println(sg.create(partition) + " was NOT added.");
								}			
						 	}
						}
					} else{
						// TO DO
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
				//System.err.println(fragment.getKey());
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
//				System.err.println(fragment.getKey());
//				System.err.println(this.sGen.create(fragment.getValue()));
				int z = (int) Math.abs(AtomContainerManipulator.getTotalFormalCharge(fragment.getValue()));
//				System.err.println("Z = " + z);
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
//			fragCondition.print();
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
//			mra.print();
//			System.out.println("fragmentMasses: " + fragmentMassChargeRatios);
//			System.out.println("mra.getRelativeAbundances(): " + mra.getRelativeAbundances());
			
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
						peaks.add(String.format("%.5f %.1f %3d (1)",mass, mra.getRelativeAbundances().get(massesToLabel.get(mass).get(0)), labels.indexOf(massesToLabel.get(mass).get(0))));
					} else{
						
						
						ArrayList<String> indexes = new ArrayList<String>();
						ArrayList<String> scores = new ArrayList<String>();
						for(String x : massesToLabel.get(mass)){
							if(mra.getRelativeAbundances().containsKey(x)){
								indexes.add(String.valueOf(labels.indexOf(x)));
								scores.add("1");
							}
						}
						peaks.add(String.format("%.5f %.1f %8s",mass, mra.getRelativeAbundances().get(massesToLabel.get(mass).get(0)), 
								String.join(" ", indexes)) + " (" + String.join(" ", scores) + ")");
					}
//				}
			}

				/*
				 * Add Fragment masses and structures
				 */
					
				for(Map.Entry<String, IAtomContainer> frag : fragments.entrySet()){
					frag_smiles_mass.add(String.format("%.5f", fragmentMassChargeRatios.get(frag.getKey())) + " " + this.sGen.create(fragments.get(frag.getKey())));
//					System.out.println(String.format("%.5f", fragmentMassChargeRatios.get(frag.getKey())) + " " + this.sGen.create(fragments.get(frag.getKey())));
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
//			System.out.println("The standardized molecule is " + this.sGen.create(standardized_mol));
			StructuralClass.ClassName type = StructureExplorer.findClassName(standardized_mol);
//			System.out.println("The type of this molecule is " + String.valueOf(type));
			LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(standardized_mol, type, fragCondition.getAdductName());
			
			return  annotatePeakList(fragments ,type, new FragmentationCondition(fragCondition.getAdductName(), 10));
		} else {
			
			throw new IllegalArgumentException("The collision energy must be either 10 eV, 20 eV, or 40 eV.\nPlease enter a valid collision energy");
		}
	}
	

	
	public LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> generateCfmidLikeMSPeakList(IAtomContainer molecule, String adductType) throws Exception{
		LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> annotatedPeaks = new LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>>();
		
		IAtomContainer standardized_mol = this.sExplorer.standardizeMolecule(molecule);
		StructuralClass.ClassName type = StructureExplorer.findClassName(standardized_mol);
//		System.out.println("The type of this molecule is " + String.valueOf(type));
		LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(standardized_mol, type, adductType);
		FragmentationCondition fragCondition_10 =  new FragmentationCondition(adductType, 10);
		FragmentationCondition fragCondition_20 =  new FragmentationCondition(adductType, 20);
		FragmentationCondition fragCondition_40 =  new FragmentationCondition(adductType, 40);
		
		annotatedPeaks.put(10, annotatePeakList(fragments ,type, fragCondition_10));
		annotatedPeaks.put(20,	annotatePeakList(fragments ,type, fragCondition_20));
		annotatedPeaks.put(40, annotatePeakList(fragments ,type, fragCondition_40));
		
		return annotatedPeaks;
	}

	public void saveSingleCfmidLikeMSPeakList(IAtomContainer molecule, IChemObjectBuilder bldr, String outputname, ArrayList<String> adduct_types) throws Exception{	
		IAtomContainer standardized_mol = this.sExplorer.standardizeMolecule(molecule);
		StructuralClass.ClassName type = StructureExplorer.findClassName(standardized_mol);
//		System.out.println("The type of this molecule is " + String.valueOf(type));

        if(FPLists.classSpecificFragmentationPatterns.containsKey(type)){
            saveSingleCfmidLikeMSPeakList(standardized_mol, bldr, outputname, type, adduct_types);  	
        } 
        else{
        	
        	if(type == ClassName.GLYCEROLIPIDS || type == ClassName.GLYCEROPHOSPHOLIPIDS || type == ClassName.SPHINGOLIPIDS 
        			|| type == ClassName.CERAMIDE_1_PHOSPHATES || type == ClassName.DIPHOSPHORYLATED_HEXAACYL_LIPID_A ||
        			type == ClassName.SULFATIDES || type == ClassName.FATTY_ACID_ESTERS_OF_HYDROXYL_FATTY_ACIDS ||
        			type == ClassName.ETHER_LIPIDS
        			){            	
    			System.out.println("\nSTATUS REPORT = 4\nThe compound belongs to the lipid class of " + type + ", which is not covered "
    					+ "in the current version of the fragmenter.");       		
        	}
        	else{
    			System.out.println("\nSTATUS REPORT = 5\nInvalid chemical class. The query compound does not belong to any of the classes covered "
    					+ "in the current version of the fragmenter.");         		
        	}
        		


        }
	}
	
	public void saveSingleCfmidLikeMSPeakList(IAtomContainer molecule, IChemObjectBuilder bldr, String outputname, StructuralClass.ClassName type, ArrayList<String> adduct_types) throws Exception{
		ArrayList<String> adduct_types_valid = new ArrayList<String>();
		ArrayList<String> adduct_types_invalid = new ArrayList<String>();
		for(String adduct : adduct_types){
			if(FPLists.classSpecificFragmentationPatterns.get(type).keySet().contains(adduct)){
				adduct_types_valid.add(adduct);
			}else{
				adduct_types_invalid.add(adduct);
			}
		}
		System.out.println(adduct_types_valid.size() < adduct_types.size());
		if(adduct_types_valid.size() < adduct_types.size()){
			if(adduct_types.size() == 1){			
				System.out.println(
					"\nSTATUS REPORT = 3\nThe following adducts are not covered for the class of " + type + ": " +
							Arrays.toString(adduct_types_invalid.toArray()) + "\n"
									+ "Next Step: Predicting MS spectra for covered adduct types: "
									+ Arrays.toString(FPLists.classSpecificFragmentationPatterns.get(type).keySet().toArray())
						);
				saveSingleCfmidLikeMSPeakList(molecule, bldr, type, outputname, false);
			}else{
				System.out.println(
						"\nSTATUS REPORT = 2\nThe following adducts are not covered for the class of " + type + ": " +
								Arrays.toString(adduct_types_invalid.toArray()) + "\n"
										+ "Next Step: Predicting MS spectra for remaining adduct types: "
										+ Arrays.toString(adduct_types_valid.toArray())
				);
				for(String adt : adduct_types_valid){
					System.out.println("Generating peak list for the adduct type: " + String.valueOf(adt));
					LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> annotatedPeaks = new LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>>();
					LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(molecule, type, adt);
					FragmentationCondition fragCondition_10 =  new FragmentationCondition(adt, 10);
					FragmentationCondition fragCondition_20 =  new FragmentationCondition(adt, 20);
					FragmentationCondition fragCondition_40 =  new FragmentationCondition(adt, 40);
					
					annotatedPeaks.put(10, annotatePeakList(fragments ,type, fragCondition_10));
					annotatedPeaks.put(20,	annotatePeakList(fragments ,type, fragCondition_20));
					annotatedPeaks.put(40, annotatePeakList(fragments ,type, fragCondition_40));
					
					String[] outputnameSplit = outputname.split(Pattern.quote("."));
//					System.err.println(outputnameSplit.length);
					String destination = outputnameSplit[0];
					if(outputnameSplit.length>0){
						destination = outputnameSplit[0] + "_" + adt + "." + outputnameSplit[1];
					}
					else{
						destination = outputnameSplit[0] + "_" + adt;
					}
					
					saveSingleCfmidLikeMSAnnotatedPeakList(annotatedPeaks, adt, destination);
					System.out.println("Saving spectrum to " + destination);
				}				
			}
			
		} else{
				System.out.println("\nSTATUS REPORT = 1\nEach specified adduct is covered for the class of " + type + ". "
						+ "Next Step: Predicting MS spectra for the following adduct types: "
						+ Arrays.toString(adduct_types_valid.toArray())
				);
				for(String adduct : adduct_types){
					System.out.println("Generating peak list for the adduct type: " + String.valueOf(adduct));
					LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> annotatedPeaks = new LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>>();
					LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(molecule, type, adduct);
//					System.out.println(fragments.keySet().size());
//					for(String x : fragments.keySet()){
//						System.out.println(x);
//						System.out.println(fragments.get(x)==null);
//					}
					FragmentationCondition fragCondition_10 =  new FragmentationCondition(adduct, 10);
					FragmentationCondition fragCondition_20 =  new FragmentationCondition(adduct, 20);
					FragmentationCondition fragCondition_40 =  new FragmentationCondition(adduct, 40);
					
					annotatedPeaks.put(10, annotatePeakList(fragments ,type, fragCondition_10));
					annotatedPeaks.put(20,	annotatePeakList(fragments ,type, fragCondition_20));
					annotatedPeaks.put(40, annotatePeakList(fragments ,type, fragCondition_40));
					
					String[] outputnameSplit = outputname.split(Pattern.quote("."));
//					System.err.println(outputnameSplit.length);
					String destination = outputnameSplit[0];
					if(outputnameSplit.length>0){
						destination = outputnameSplit[0] + "_" + adduct + "." + outputnameSplit[1];
					}
					else{
						destination = outputnameSplit[0] + "_" + adduct;
					}					
					saveSingleCfmidLikeMSAnnotatedPeakList(annotatedPeaks, adduct, destination);
					System.out.println("Saving spectrum to " + destination );
				}				
		}
	}
	
	
	

	public void saveSingleCfmidLikeMSPeakList(IAtomContainer molecule, IChemObjectBuilder bldr, String outputname) throws Exception{
		IAtomContainer standardized_mol = this.sExplorer.standardizeMolecule(molecule);
		StructuralClass.ClassName type = StructureExplorer.findClassName(standardized_mol);
		System.out.println("The type of this molecule is " + String.valueOf(type));

		saveSingleCfmidLikeMSPeakList(standardized_mol, bldr, type, outputname, false);
		
//		if(FPLists.classSpecificFragmentationPatterns.containsKey(type)){			
//			for(String adduct : FPLists.classSpecificFragmentationPatterns.get(type).keySet()){
//				System.out.println("Generating peak list for the adduct type: " + String.valueOf(adduct));
//				LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> annotatedPeaks = new LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>>();
//				LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(standardized_mol, type, adduct);
//				FragmentationCondition fragCondition_10 =  new FragmentationCondition(adduct, 10);
//				FragmentationCondition fragCondition_20 =  new FragmentationCondition(adduct, 20);
//				FragmentationCondition fragCondition_40 =  new FragmentationCondition(adduct, 40);
//				
//				annotatedPeaks.put(10, annotatePeakList(fragments ,type, fragCondition_10));
//				annotatedPeaks.put(20,	annotatePeakList(fragments ,type, fragCondition_20));
//				annotatedPeaks.put(40, annotatePeakList(fragments ,type, fragCondition_40));
//				
//				saveSingleCfmidLikeMSAnnotatedPeakList(annotatedPeaks, adduct, (outputname + "_" + adduct));
//				System.out.println("Saving spectrum to " + (outputname + "_" + adduct) );
//			}			
//		}
	}

	
	public void saveSingleCfmidLikeMSPeakList(IAtomContainer molecule, IChemObjectBuilder bldr, StructuralClass.ClassName type, String outputname, boolean standardize) throws Exception{
		
		IAtomContainer standardized_mol = molecule;
		if(standardize){
			standardized_mol = this.sExplorer.standardizeMolecule(molecule);
		}
//		System.out.println("HEY");
		if(FPLists.classSpecificFragmentationPatterns.containsKey(type)){
			for(String adduct : FPLists.classSpecificFragmentationPatterns.get(type).keySet()){
				System.out.println("Generating peak list for the adduct type: " + String.valueOf(adduct));
				LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>> annotatedPeaks = new LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>>();
				LinkedHashMap<String, IAtomContainer> fragments = fragmentMolecule(standardized_mol, type, adduct);
				FragmentationCondition fragCondition_10 =  new FragmentationCondition(adduct, 10);
				FragmentationCondition fragCondition_20 =  new FragmentationCondition(adduct, 20);
				FragmentationCondition fragCondition_40 =  new FragmentationCondition(adduct, 40);
				
				annotatedPeaks.put(10, annotatePeakList(fragments ,type, fragCondition_10));
				annotatedPeaks.put(20,	annotatePeakList(fragments ,type, fragCondition_20));
				annotatedPeaks.put(40, annotatePeakList(fragments ,type, fragCondition_40));
				
				String[] outputnameSplit = outputname.split(Pattern.quote("."));
//				System.err.println(outputnameSplit.length);
				String destination = outputnameSplit[0];
				if(outputnameSplit.length>0){
					destination = outputnameSplit[0] + "_" + adduct + "." + outputnameSplit[1];
				}
				else{
					destination = outputnameSplit[0] + "_" + adduct;
				}

				saveSingleCfmidLikeMSAnnotatedPeakList(annotatedPeaks, adduct, destination);

				System.out.println("Saving spectrum to " + destination);
			}			
		}
		else{
        	if(type == ClassName.GLYCEROLIPIDS || type == ClassName.GLYCEROPHOSPHOLIPIDS || type == ClassName.SPHINGOLIPIDS 
       			|| type == ClassName.CERAMIDE_1_PHOSPHATES || type == ClassName.DIPHOSPHORYLATED_HEXAACYL_LIPID_A ||
        			type == ClassName.SULFATIDES || type == ClassName.FATTY_ACID_ESTERS_OF_HYDROXYL_FATTY_ACIDS ||
        			type == ClassName.ETHER_LIPIDS
        			){            	
    			System.out.println("\nSTATUS REPORT = 4\nThe compound belongs to the lipid class of " + type + ", which is not covered "
    					+ "in the current version of the fragmenter.");       		
        	}
        	else{
    			System.out.println("\nSTATUS REPORT = 5\nInvalid chemical class. The query compound does not belong to any of the classes covered "
    					+ "in the current version of the fragmenter.");         		
       	}
        }
		
		

	}
	
	
	
	public void saveSingleCfmidLikeMSAnnotatedPeakList(LinkedHashMap<Integer, LinkedHashMap<String, ArrayList<String>>>annotatedPeaks, String adductType, String outputname) throws IOException{
		if(annotatedPeaks != null) {
			
			FileWriter fw = new FileWriter(outputname);
//			FileWriter fw = new FileWriter("data/" + outputname.replace("/", "_") + ".log");
			BufferedWriter bw = new BufferedWriter(fw);
			int level = 0;
			int fragIndex = 0;
			
			for(Map.Entry<Integer, LinkedHashMap<String, ArrayList<String>>> peaks : annotatedPeaks.entrySet()){
//				System.out.println("Energy: " + peaks.getKey());
//				System.out.println("peaks: " + peaks.getValue());
//				System.out.println("Nr. of peaks: " + peaks.getValue().get("peaks_list").size());
				bw.write("energy" + level);
				bw.newLine();
				level++;
				for(String s : peaks.getValue().get("peaks_list")){
//					System.out.println();
					bw.write(s);
					bw.newLine();
				}
			}
			
			bw.newLine();
			for( String fg : annotatedPeaks.get(10).get("fragments")){
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
	}
	
	
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
//					System.out.println("Energy: " + peaks.getKey());
//					System.out.println("peaks: " + peaks.getValue());
//					System.out.println("Nr. of peaks: " + peaks.getValue().get("peaks_list").size());
					bw.write("energy" + level);
					bw.newLine();
					level++;
					for(String s : peaks.getValue().get("peaks_list")){
//						System.out.println();
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

