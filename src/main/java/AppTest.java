
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
//import ambit2.

import org.openscience.cdk.smiles.SmilesGenerator;

import wishartlab.cfmid_plus.fragmentation.Fragmenter;
import wishartlab.cfmid_plus.molecules.StructureExplorer;

/**
 * Hello world!
 *
 */
public class AppTest 
{
    public static void main( String[] args ) throws Exception
    {
        SmilesParser sp =  new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator sg = new SmilesGenerator().unique();
        StructureExplorer se = new StructureExplorer();
        
////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O[C@H]3[C@@H]([C@H]([C@H]([C@H](O3)CO)O)O)NC(=O)C)O[C@@]4(C[C@@H]([C@H](C(O4)[C@@H]([C@@H](CO)O)O)NC(=O)C)O)C(=O)O)O)O)O)[C@@H](/C=C/CCCCCCCCCCCCC)O");
////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@H]1C([C@H]([C@@H](C(O1)CO)O[C@H]2C([C@H]([C@H](C(O2)CO)O)O[C@@]3(CC([C@H](C(O3)[C@@H]([C@@H](CO)O)O)NC(=O)C)O)C(=O)O)O)O)O)[C@@H](/C=C/CCCCCCCCCCCCC)O");
////        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1OC(CO)[C@H](O)[C@H](O[C@@]2(CC(O)[C@@H](NC(C)=O)C(O2)[C@H](O)[C@H](O)CO)C(O)=O)C1O)[C@H](O)C=CCCCCCCCCCCCCC");
//
////        IAtomContainer molecule = sp.parseSmiles("OCC(O)COP(O)(=O)OCC(CO(C(=O)C))O(C(=O)CC)");
//        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@H](CO)O)OC(=O)CCCCCCCCC/C=C\\CCCCCC");
//        IAtomContainer molecule = sp.parseSmiles("CCCCC(=O)OC[C@H](COP(O)(=O)OC[C@@H](O)CO)OC(=O)CCCC");
//        IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCC(O)COP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC");
        	IAtomContainer molecule = sp.parseSmiles("CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC(O)COP([O-])(=O)OC[C@@H](COC(=O)CCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC");
        	IAtomContainer molecule2 = sp.parseSmiles("CCCCCCCC\\C=C/CCCCCCCC(=O)O[C@H]1CC[C@]2(C)C3CC[C@]4(C)[C@H](CC[C@H]4C3CCC2=C1)[C@H](C)CCCC(C)C");
        	
        	IAtomContainer molecule3 = sp.parseSmiles("[O-]C(CCCCCCCCCCCCCCC)=O");
        	
        	System.out.println("Molecule3 standardized: " + sg.create(se.standardizeMolecule(molecule3)) + "\n\n");
        	
        	
//         System.out.println("[#6:7]-[#8;X2:1][P;X4:2]([#8;X1-:4])(=[O;X1:5])[#8;X2:3]-[#6:6]".replaceAll("\\:[0-9]+\\]", "\\]"));
//         IAtomContainer molecule_c = se.standardizeMolecule(molecule);
//         AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule_c);
//         CDKHydrogenAdder.getInstance(molecule_c.getBuilder()).addImplicitHydrogens(molecule_c);
//         String molecule_cs = sg.create(molecule_c);
//         System.out.println("Standardized molecule: "+ molecule_cs);
        	
        // CCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCC(O)COP(O)(=O)OC[C@@H](COC(=O)CCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCC
        //        
//        String molecule_c = sg.create(StructureExplorer.cleanMolecule(molecule,new SMIRKSManager(SilentChemObjectBuilder.getInstance())));
//        System.out.println("Cleaned molecule: "+ molecule_c);
//        System.out.println(AtomContainerManipulator.getNaturalExactMass(molecule));
//        System.out.println("Type of this molecule: " + StructureExplorer.fin dClassName(molecule));
        Fragmenter fr = new Fragmenter();

    	

        
////        IAtomContainerSet fragments = fr.fragmentMolecule(molecule);
////        
////        for (int i = 0; i < fragments.getAtomContainerCount(); i++){
////        	System.out.println(sg.create(fragments.getAtomContainer(i)));
////        }
//        LinkedHashMap<String,ArrayList<String>> peaksResults = fr.generateMSPeakList(molecule2);
//        ArrayList<String> peaks  =  peaksResults;
//        fr.generateGraphFromPeakListAndSaveToDir(peaksResults);   
//        System.out.println("WRITING  PEAKS");
//        for(Map.Entry<String, ArrayList<String>> p : peaksResults.entrySet() ){
//        	for (int k =0; k < p.getValue().size(); k++){
//        		System.out.println(p.getValue().get(k));
//        	}
//        }
        
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		
		System.out.println("MOLECULE WITH EXPLICIT HYDROGENS: " + sg.create(molecule));      
        
        IChemObjectBuilder bldr = DefaultChemObjectBuilder.getInstance();
//        ArrayList<LinkedHashMap<String,ArrayList<String>>> allPeaks = fr.generatePeakListFromSDF("data/lipids.sdf", bldr);
//        fr.generateGraphFromPeakListAndSaveToDir(allPeaks.get(0));
              
//         fr.generateSpectraFromPeakListAndSaveToDir(allPeaks);

        
//      IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
//      fr.savePeakListFromSDF("data/ymdb_test.sdf",bldr,"data/ymdb_spectra.msp"); 
		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/test3.sdf",bldr);
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/DG_IDS/DG_IDS.sdf",bldr);
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/TG_IDS/TG_IDS.sdf",bldr);
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PA_IDS/PA_IDS.sdf",bldr);
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/Lyso-PC_IDS/Lyso-PC_IDS.sdf",bldr);
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/Lyso-PE_IDS/Lyso-PE_IDS.sdf",bldr);
//		fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/MG_IDS/MG_IDS.sdf",bldr); 
//	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PI_IDS/PI_IDS.sdf",bldr);
//	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PG_IDS/PG_IDS.sdf",bldr);
//	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PE_IDS/PE_IDS.sdf",bldr);
//	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PS_IDS/PS_IDS.sdf",bldr);   
//	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/PC_IDS/PC_IDS.sdf",bldr);
//	    fr.saveSingleCfmidLikeMSPeakListFromSDF("data/LipidIds/CL_IDS/CL_IDS.sdf",bldr);
      
      
      
       


      
        
//        System.out.println(f.size());
//        System.out.println(f.get(f.size()-1));
//        
//        
//        System.out.println("\nH: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("H")));
//        System.out.println("Li: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("[Li]")));
//        System.out.println("H2O: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("[H]O[H]")));
//        System.out.println("O: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("O")));
//        System.out.println("C: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("C")));
//        System.out.println(sp.parseSmiles("C"));
//        System.out.println(AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("CC(=O)N[C@@H]1[C@@H](O)C[C@@H](O[C@H]1[C@H](O)[C@H](O)CO)C(O)=O")));
//        System.out.println(AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("OCC(O)CO")));
//        
//        //196.41083318690977
////        10.0794075382579
//        
//        
//        System.out.println("C6H10O5: " +  AtomContainerManipulator.getNaturalExactMass(sp.parseSmiles("OC1CC(O)C(O)C(O)C1O")));
//        

        
//        fr.generateMSSpecImage();
    }
}

