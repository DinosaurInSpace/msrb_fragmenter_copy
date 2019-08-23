package wishartlab.cfmid_plus.molecules;

/**
 * @author Yannick Djoumbou Feunang
 *
 */

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
//import org.openscience.cdk.smiles.SmilesGenerator;
//import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
//import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;


public class StructureExplorer extends AtomContainer  {
	ArrayList<SMIRKSReaction> standardizationReactions = new ArrayList<SMIRKSReaction>();
	SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
//	SmilesGenerator smiGen = new SmilesGenerator().isomeric();
	
	public StructureExplorer(){
		setUpStandardizer();
		this.smrkMan.setFlagFilterEquivalentMappings(true);		
	}
	
	
	
	public static boolean containsSmartsPattern(IAtomContainer molecule, String pattern) throws SMARTSException, CDKException{
		SmartsPatternCDK smartsPattern = new SmartsPatternCDK(pattern);
		return smartsPattern.hasSMARTSPattern(molecule) > 0;
		
//		SMARTSQueryTool sqt = new SMARTSQueryTool(pattern, SilentChemObjectBuilder.getInstance());
//		return sqt.matches(molecule);

	}
	
	// partition the return products, as the products of each reaction are returned as a 
	// single fragment
	

	public static IAtomContainerSet partition(IAtomContainer molecule)
			throws CDKException {
		IAtomContainerSet partitions = new AtomContainerSet();
		IAtomContainerSet ms = ConnectivityChecker.partitionIntoMolecules(molecule);
		// System.out.println("MS: " + ms.getAtomContainerCount() + " molecule(s).");

		for (int k = 0; k < ms.getAtomContainerCount(); k++) {
			IAtomContainer current_partition = ms.getAtomContainer(k);
			AtomContainerManipulator
					.percieveAtomTypesAndConfigureAtoms(current_partition);
			for (IAtom atom : current_partition.atoms())
				if (atom.getImplicitHydrogenCount() == null) {
					atom.setImplicitHydrogenCount(0);
				}
//			AtomContainerManipulator.suppressHydrogens(current_partition);
			partitions.addAtomContainer(current_partition);
		}
		//System.out.println("MS PARTITIONS: " + partitions.getAtomContainerCount()
		//		+ " partition(s).");
		return partitions;
	}
	
	public static StructuralClass.ClassName findClassName(IAtomContainer molecule) throws SMARTSException, CDKException{
		StructuralClass.ClassName type = StructuralClass.ClassName.NIL;
		
		if(isFAHFA(molecule)){
			type = StructuralClass.ClassName.FATTY_ACID_ESTERS_OF_HYDROXYL_FATTY_ACIDS;			
		}

		else{
			for(Map.Entry<StructuralClass.ClassName, String > s : StructuralClass.backbones.entrySet()){
//				System.out.println(s);
				if(containsSmartsPattern(molecule, s.getValue())){
//					System.err.println("THIS MOLECULE BELONGS TO THE " + s.getKey() +"\n");
					type = s.getKey();
					break;
				}
			}			
		}
		return type;
	}
	
	public static boolean isFAHFA(IAtomContainer molecule) throws SMARTSException, CDKException{
		boolean isFahfa = false;
//		String pattern1 = "";
//		String pattern2 = "";
//		String carboxyl_head = "[H][#6;R0]-[#6;R0]([H])-[#6;R0]([H])-[#6;X3]([#8;A;X1-,X2H1])=[O;X1]";
//		String ester_head    = "[H][#6;R0]-[#6;R0]([H])-[#6;R0]([H])-[#6;X3](=[O;X1])[#8;A;X2R0][C;R0]([H])([#6;R0])[#6;R0][H]";

		// Minimum C-10 acyl chains
		String fatty_acid_head  = "[H][#6;R0]-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-[#6;X3]([#8;A;X1-,X2H1])=[O;X1]";
		String fatty_ester_head = "[H][#6;R0]-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-,=[#6;R0]([H])-[#6;X3](=[O;X1])[#8;A;X2R0][C;R0]([H])([#6;R0])[#6;A;X4R0;H2,H3]";

		// nr fragments = 1
		// nr of oxygen atoms = 4
		
//		System.out.println(!isMixture(molecule));
//		System.out.println(containsSmartsPattern(molecule, fatty_acid_head));
//		System.out.println(containsSmartsPattern(molecule, fatty_ester_head));
//		System.out.println(numberOfAtomWithAtomicNumber(molecule, 16));
//		for(int l = 0; l <molecule.getAtomCount(); l++){
//			System.out.println(molecule.getAtom(l).getAtomTypeName() + "\t" + String.valueOf(molecule.getAtom(l).getAtomicNumber()));
//		}
		
		if( (!isMixture(molecule) ) && containsSmartsPattern(molecule, fatty_acid_head) && containsSmartsPattern(molecule, fatty_ester_head) &&
				numberOfAtomWithAtomicNumber(molecule, 8) == 4 ){
			isFahfa = true;
		}
		
		return isFahfa;
	}

	
	public static boolean isMixture(IAtomContainer molecule) throws CDKException{
		// compound is not a mixture (checkConnectivity returns 2 or more atomContainers)
		boolean mixture = ConnectivityChecker.partitionIntoMolecules(molecule).getAtomContainerCount()>1;
		return mixture;	
	}
	
	/**
	 * 
	 * @param molecule
	 * @param atomicNumber
	 * @return the number of atoms with the given atomic number 
	 */
	public static int numberOfAtomWithAtomicNumber(IAtomContainer molecule,int atomicNumber){
		int atCount = 0;
		
		for(int l = 0; l <molecule.getAtomCount(); l++){
			if(molecule.getAtom(l).getAtomicNumber() == atomicNumber)
				atCount++;
		}
		return atCount;
	}
	
	/**
	 * Given two molecules, find out whether they are equal. This function uses
	 * InChI strings generated by CDK's InChIGeneratorFactory class as a
	 * comparator.
	 * 
	 * @param mol1
	 * @param mol2
	 * @return a boolean value that specifies whether the two molecules are
	 *         equal.
	 * @throws CDKException
	 */
	public static boolean inchiEqualityHolds(IAtomContainer mol1, IAtomContainer mol2)
			throws CDKException {
//		boolean equal = false;

		// Generate factory - throws CDKException if native code does not load
		InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
		// Get InChIGenerator
		InChIGenerator gen1 = factory.getInChIGenerator(mol1);
		InChIGenerator gen2 = factory.getInChIGenerator(mol2);

		String inchikey1 = gen1.getInchiKey();
		String inchikey2 = gen2.getInchiKey();
		
		return inchikey1.equals(inchikey2);
	}

	
	/**
	 * Given a molecule (or AtomContainer) M and an AtomContainerSet S, check
	 * whether M is included in the set S.
	 * 
	 * @param molecules
	 * @param mol
	 * @return
	 * @throws CDKException
	 */
	public static boolean atomContainerInclusionHolds(IAtomContainerSet molecules,
			IAtomContainer mol) throws CDKException {
		boolean inclusion = false;
		for (int i = 0; i < molecules.getAtomContainerCount(); i++) {

			if (inchiEqualityHolds(mol, molecules.getAtomContainer(i))) {
				inclusion = true;
				break;
			}
		}
		return inclusion;
	}
	
	/**
	 * Given an AtomContainerSet S, return an AtomContainerSet Su with only
	 * unique compounds of S
	 * 
	 * @param molecules
	 * @return A subset of the input AtomContainerSet S with only unique
	 *         compounds
	 * @throws Exception
	 */
	
	public static IAtomContainerSet uniquefy(IAtomContainerSet molecules)
			throws Exception {
		if (molecules != null && (!molecules.isEmpty()) && molecules.getAtomContainerCount() > 1) {
			
			IAtomContainerSet uniqueContainer = DefaultChemObjectBuilder.getInstance()
					.newInstance(IAtomContainerSet.class);
			
			uniqueContainer.addAtomContainer(molecules.getAtomContainer(0));
			
			for (int i = 1; i < molecules.getAtomContainerCount(); i++) {
				if (! ( (molecules.getAtomContainer(i) == null) || atomContainerInclusionHolds(uniqueContainer,
						molecules.getAtomContainer(i) ))) {
					uniqueContainer.addAtomContainer(molecules.getAtomContainer(i));
				}
			}

			return uniqueContainer;
		}

		else
			return molecules;

	}

	public static IAtomContainer cleanMolecule(IAtomContainer molecule, SMIRKSManager smrkMan) throws Exception{
		IAtomContainer cleaned_mol =  molecule;
		StructuralClass.ClassName type = findClassName(molecule);
		String cleaningSmirks;
		
		IAtomContainer mol_no_hydro = AtomContainerManipulator.removeHydrogens(molecule);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol_no_hydro);
//		SmilesGenerator sg = new SmilesGenerator().unique();
//		String molecule_c = sg.create(mol_no_hydro);
//		System.out.println(molecule_c);
		
		switch (type) {
		case SPHINGOMYELINS:
			cleaningSmirks = "[H][#8:1]P(=O)([#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#8;R0][#6;A;H2X4][#6;A;H1X4]([#7;A;H1X3][#6;R0](=O)[#6;A;H2X4][#6;A;H2X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H3X4]>>[#6;A;H2X4][#6;A;H2X4][#6;R0](=O)[#7;A;H1X3][#6;A;H1X4]([#6;A;H2X4][#8;R0]P([#8;X1-:1])(=O)[#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H3X4]";
			SMIRKSReaction sReaction = smrkMan.parse(cleaningSmirks);
			IAtomContainerSet cleaned = smrkMan.applyTransformationWithSingleCopyForEachPos(mol_no_hydro, null, sReaction);
			if(cleaned.getAtomContainerCount()>1){
				cleaned_mol = cleaned.getAtomContainer(0);
			}
		}

		return cleaned_mol;

	}

	public IAtomContainer standardizeMolecule(IAtomContainer molecule) throws Exception{
//		System.out.println("Before standardization: " + this.smiGen.create(molecule));
		IAtomContainer stMol =  molecule.clone();
		

		
//		System.out.println("Before standardization: " + this.smiGen.create(stMol));

		for(int i = 0; i < this.standardizationReactions.size(); i++){
//			System.out.println(this.standardizationReactions.get(i).reactantsSmarts);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(stMol);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(stMol);
			while(StructureExplorer.compoundMatchesReactionConstraints(this.standardizationReactions.get(i), stMol)){
				this.smrkMan.applyTransformation(stMol, this.standardizationReactions.get(i));
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(stMol);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(stMol);
			}
		}		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(stMol);
		CDKHydrogenAdder.getInstance(stMol.getBuilder()).addImplicitHydrogens(stMol);
		
		return stMol;
	}
	
	
//	public IAtomContainer standardizeMolecule(IAtomContainer molecule) throws Exception{
//		IAtomContainer stMol =  molecule.clone();
//		
//
//		
//		System.out.println("Before standardization: " + this.smiGen.create(stMol));
//
//		for(int i = 0; i < this.standardizationReactions.size(); i++){
//			System.out.println(this.standardizationReactions.get(i).reactantsSmarts);
//			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(stMol);
//			AtomContainerManipulator.convertImplicitToExplicitHydrogens(stMol);
//			
////			if(StructureExplorer.compoundMatchesReactionConstraints(this.standardizationReactions.get(i), stMol)){
////				this.smrkMan.apply
////				.applyTransformationWithCombinedOverlappedPos(stMol, null,this.standardizationReactions.get(i));
////				
////			}
//			
//			IAtomContainerSet s = this.smrkMan.applyTransformationWithCombinedOverlappedPos(stMol, null,this.standardizationReactions.get(i));
//			
//			for(IAtomContainer a : s.atomContainers()){
//				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(a);
//				System.err.println(this.smiGen.create(a));
//			}
//		}		
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(stMol);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(stMol);
//		
//		System.out.println("After standardization: " + this.smiGen.create(stMol));
//		
//		return stMol;
//	}
	
	
	
	
	public static boolean compoundMatchesReactionConstraints(SMIRKSReaction reaction,
			IAtomContainer molecule) throws SMARTSException, CDKException, IOException {
		boolean match = true;

		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//		CDKHydrogenAdder.getInstance(molecule.getBuilder()).addImplicitHydrogens(molecule);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
//		System.out.println(reaction.reactantsSmarts.replaceAll("\\:[0-9]+\\]", "\\]"));
		
		if (match) {				
//				System.out.println(reaction.reactantsSmarts.replaceAll("\\:[0-9]+\\]", "\\]"));
				Pattern smp = SmartsPattern.create(reaction.reactantsSmarts.replaceAll("\\:[0-9]+\\]", "\\]"), bldr);
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
				CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
				adder.addImplicitHydrogens(molecule);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
//				SmilesGenerator smiGen = new SmilesGenerator().isomeric();
//				System.err.println(smiGen.create(molecule));
				boolean status = smp.matches(molecule);
//				System.out.println("status: " + status);				
//				SmartsPatternCDK smartsPattern = new SmartsPatternCDK(reaction.reactantsSmarts.replaceAll("\\:[0-9]+\\]", "\\]"));
//				System.out.println("status via SmartsPatternCDK: " + smartsPattern.match(molecule));
//				boolean status = (smartsPattern.match(molecule) > 0);
				if (!status) {
					match = false;
				}
		}
//		System.out.println("PARTIAL MATCH: " + match);

		return match;

	}
	
	public static int fragmentCount(Pattern smp, IAtomContainer ac){
		int fc = 0;
		
		if (smp.matches(ac)){
			fc = smp.matchAll(ac).count();
		}
		
		return fc;
	}
	
	
	private void setUpStandardizer(){		
		SMIRKSReaction oxygenTransform = smrkMan.parse("[#8;X1-:2]!@-[*:1]>>[H][#8;X2:2]!@-[*:1]");
//		SMIRKSReaction phosphateTransform = smrkMan.parse("[#8;X1-:1][P;X4:2](=[O;X1:3])([#8;X2:4]-[*,#1:5])[#8;X2:6]-[*,#1:7]>>[H][#8;X2:1][P;X4:2](=[O;X1:3])([#8;X2:4]-[*,#1:5])[#8;X2:6]-[*,#1:7]");
//		SMIRKSReaction phosphateTransform = smrkMan.parse("[#6:1]-[#8;X2:2][P;X4:3]([#8;X1-:4])(=[O;X1:5])[#8;X2:6]-[#6:7]>>[H][#8;X2:4][P;X4:3](=[O;X1:5])([#8;X2:2]-[#6:1])[#8;X2:6]-[#6:7]");
//		SMIRKSReaction phosphateTransform2 = smrkMan.parse("[#1:7]-[#8;X2:1][P;X4:2]([#8;X1-:4])(=[O;X1:5])[#8;X2:3]-[#6:6]>>[H][#8;X2:4][P;X4:2](=[O;X1:5])([#8;X2:1]-[#1:7])[#8;X2:3]-[#6:6]");				
//		SMIRKSReaction sphingomyelinTransform = smrkMan.parse("[H][#8:1]P(=O)([#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#8;R0][#6;A;H2X4][#6;A;H1X4]([#7;A;H1X3][#6;R0]"
//				+ "(=O)[#6;A;H2X4][#6;A;H2X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]"
//				+ "[#6;A;H2X4][#6;A;H2X4][#6;A;H3X4]>>[#6;A;H2X4][#6;A;H2X4][#6;R0](=O)[#7;A;H1X3][#6;A;H1X4]([#6;A;H2X4][#8;R0]P([#8;X1-:1])(=O)[#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])"
//				+ "([#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4]"
//				+ "[#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H3X4]");
		
		SMIRKSReaction sphingomyelinTransform = smrkMan.parse("[H][#6;A;X4]([H])[#6;A;X4]([H])([H])[#6;R0](=O)[#7;A;X3]([H])[#6;A;X4]([H])([#6;A;X4]([H])([H])[#8;R0]P([#8;X1-:1])(=O)[#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#6;A;X4]([H])([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[C;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[H]>>[H][#8:1]P(=O)([#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#8;R0][#6;A;X4]([H])([H])[#6;A;X4]([H])([#7;A;X3]([H])[#6;R0](=O)[#6;A;X4]([H])([H])[#6;A;X4]([H])[H])[#6;A;X4]([H])([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[C;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[#6;A;X4]([H])([H])[H]");
		SMIRKSReaction carboxylTransform = smrkMan.parse("[#8;X1-:1]-[#6:2]([#6,#1;A:3])=[O;X1:4]>>[H][#8;X2:1]-[#6:2]([#6,#1;A:3])=[O;X1:4]");
		SMIRKSReaction sulfonylTransform = smrkMan.parse("[#6:1][S;X4:2]([#8;X1-:5])(=[O;X1:3])=[O;X1:4]>>[H][#8;X2:5][S;X4:2]([#6:1])(=[O;X1:3])=[O;X1:4]");
		SMIRKSReaction carnitineTransform = smrkMan.parse("[#6;A;H3X4:1][N;X4+:2]([#6;A;H3X4:3])([#6;A;H3X4:4])[#6;A;H2X4:5][#6;A;H1X4:6]([#6;A;H2X4:7][#6:8](-[#8;X1-:9])=[O;X1:10])[#8;X2:11]-[#6:12]([#6,#1;A:13])=[O;X1:14]>>[H][#8;X2:9]-[#6:8](=[O;X1:10])[#6;A;H2X4:7][#6;A;H1X4:6]([#6;A;H2X4:5][N;X4+:2]([#6;A;H3X4:1])([#6;A;H3X4:3])[#6;A;H3X4:4])[#8;X2:11]-[#6:12]([#6,#1;A:13])=[O;X1:14]");
		
		
//		System.out.println(phosphateTransform.reactantsSmarts);
//		System.out.println(phosphateTransform2.reactantsSmarts);
		this.standardizationReactions.add(oxygenTransform);
//		this.standardizationReactions.add(phosphateTransform);
//		this.standardizationReactions.add(phosphateTransform2);
//		this.standardizationReactions.add(carboxylTransform);
		this.standardizationReactions.add(sulfonylTransform);
		this.standardizationReactions.add(sphingomyelinTransform);
//		this.standardizationReactions.add(carnitineTransform);
		
	}

	public static double getMajorIsotopeMass(IAtomContainer molecule){
//		Double mass;		
//        MolecularFormulaManipulator mfm = new MolecularFormulaManipulator();
        IMolecularFormula m = MolecularFormulaManipulator.getMolecularFormula(molecule);		
		return MolecularFormulaManipulator.getMajorIsotopeMass(m);
	}
	
	
}

