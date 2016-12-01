package wishartlab.cfmid_plus.molecules;

/**
 * @author Yannick Djoumbou Feunang
 *
 */

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

import org.openscience.cdk.AtomContainer;
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
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.smiles.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;


public class StructureExplorer extends AtomContainer  {
	ArrayList<SMIRKSReaction> standardizationReactions = new ArrayList<SMIRKSReaction>();
	SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
	
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
		System.out.println("MS: " + ms.getAtomContainerCount() + " molecule(s).");

		for (int k = 0; k < ms.getAtomContainerCount(); k++) {
			IAtomContainer current_partition = ms.getAtomContainer(k);
			AtomContainerManipulator
					.percieveAtomTypesAndConfigureAtoms(current_partition);
			for (IAtom atom : current_partition.atoms())
				if (atom.getImplicitHydrogenCount() == null) {
					atom.setImplicitHydrogenCount(0);
				}
			AtomContainerManipulator.suppressHydrogens(current_partition);
			partitions.addAtomContainer(current_partition);
		}
		System.out.println("MS PARTITIONS: " + partitions.getAtomContainerCount()
				+ " partition(s).");
		return partitions;
	}
	
	public static StructuralClass.ClassName findClassName(IAtomContainer molecule) throws SMARTSException, CDKException{
		StructuralClass.ClassName type = StructuralClass.ClassName.NIL;
		for(Map.Entry<StructuralClass.ClassName, String > s : StructuralClass.backbones.entrySet()){
			System.out.println(s);
			if(containsSmartsPattern(molecule, s.getValue())){
				System.err.println("THIS MOLECULE BELONGS TO THE " + s.getKey());
				type = s.getKey();
				break;
			}
		}
		return type;
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
		boolean equal = false;

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
		SmilesGenerator sg = new SmilesGenerator().unique();
		String molecule_c = sg.create(mol_no_hydro);
		System.out.println(molecule_c);
		
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
		IAtomContainer stMol =  molecule.clone();
		
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(stMol);

		for(int i = 0; i < this.standardizationReactions.size(); i++){
			while(StructureExplorer.compoundMatchesReactionConstraints(this.standardizationReactions.get(i), stMol)){
				this.smrkMan.applyTransformation(stMol, this.standardizationReactions.get(i));
			}
		}
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(stMol);
		CDKHydrogenAdder.getInstance(stMol.getBuilder()).addImplicitHydrogens(stMol);
		
		return stMol;
	}
	
	public static boolean compoundMatchesReactionConstraints(SMIRKSReaction reaction,
			IAtomContainer molecule) throws SMARTSException, CDKException, IOException {
		boolean match = true;

		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder.getInstance(molecule.getBuilder()).addImplicitHydrogens(molecule);
		
		if (match) {				
				Pattern smp = SmartsPattern.create(reaction.reactantsSmarts.replaceAll("\\:[0-9]+\\]", "\\]"), bldr);				
				boolean status = smp.matches(molecule);
//				System.out.println("status: " + status);				
//				SmartsPatternCDK smartsPattern = new SmartsPatternCDK(reaction.getReactantSMARTS().get(j));
//				System.out.println("status via SmartsPatternCDK: " + smartsPattern.match(molecule));
				
				if (!status) {
					match = false;
				}
		}
//		System.out.println("PARTIAL MATCH: " + match);

		return match;

	}
	
	
	private void setUpStandardizer(){
		
		SMIRKSReaction phosphateTransform = smrkMan.parse("[#6:7]-[#8;X2:1][P;X4:2]([#8;X1-:4])(=[O;X1:5])[#8;X2:3]-[#6:6]>>[H][#8;X2:4][P;X4:2](=[O;X1:5])([#8;X2:1]-[#6:7])[#8;X2:3]-[#6:6]");
		SMIRKSReaction phosphateTransform2 = smrkMan.parse("[#1:7]-[#8;X2:1][P;X4:2]([#8;X1-:4])(=[O;X1:5])[#8;X2:3]-[#6:6]>>[H][#8;X2:4][P;X4:2](=[O;X1:5])([#8;X2:1]-[#1:7])[#8;X2:3]-[#6:6]");				
		SMIRKSReaction sphingomyelinTransform = smrkMan.parse("[H][#8:1]P(=O)([#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#8;R0][#6;A;H2X4][#6;A;H1X4]([#7;A;H1X3][#6;R0]"
				+ "(=O)[#6;A;H2X4][#6;A;H2X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]"
				+ "[#6;A;H2X4][#6;A;H2X4][#6;A;H3X4]>>[#6;A;H2X4][#6;A;H2X4][#6;R0](=O)[#7;A;H1X3][#6;A;H1X4]([#6;A;H2X4][#8;R0]P([#8;X1-:1])(=O)[#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])"
				+ "([#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6]-,=[#6][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4]"
				+ "[#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H3X4]");
		//		System.out.println(phosphateTransform.reactantsSmarts);
//		System.out.println(phosphateTransform2.reactantsSmarts);
		this.standardizationReactions.add(phosphateTransform);
		this.standardizationReactions.add(phosphateTransform2);
		this.standardizationReactions.add(sphingomyelinTransform);
	}

}

