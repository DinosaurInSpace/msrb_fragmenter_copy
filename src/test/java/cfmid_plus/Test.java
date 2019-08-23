package cfmid_plus;

import org.openscience.cdk.silent.SilentChemObjectBuilder;
//import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import wishartlab.cfmid_plus.molecules.StructureExplorer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
//import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;

//import org.apache.commons.io.FilenameUtils;
//import org.openscience.cdk.AtomContainer;
//import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.interfaces.IAtomContainer;
//import org.openscience.cdk.io.SDFWriter;


public class Test {

	public static void add_mass(String tsvFileName, String outputFileName) throws Exception{
        SmilesParser sp =  new SmilesParser(SilentChemObjectBuilder.getInstance());
//        SmilesGenerator sg = new SmilesGenerator().unique();
        
		BufferedReader bRead = new BufferedReader(new FileReader(tsvFileName));
		int counter = 0;
		String line = null;

		@SuppressWarnings("resource")
		BufferedWriter bw0 = new BufferedWriter(new FileWriter(outputFileName));
		bw0.write("	challengename	PRECURSOR_MZ	ION_MODE	RT	nPeaks	NAME	SMILES	IUPAC	INCHIKEY	CSID	PC_CID	MajorIsotopicMass");
		
		
		while((line = bRead.readLine()) !=null){
			counter++;
			System.out.println(counter);
			String[] sline = line.split("\t");
			
			if(! (line.contains("challengename") || sline[7].contentEquals("NULL") )){
				IAtomContainer atc = sp.parseSmiles(sline[7].trim());
				double mass = StructureExplorer.getMajorIsotopeMass(atc);
				bw0.newLine();
				bw0.write(line + "\t" + String.valueOf(mass));				
			}
		}
		
		bw0.close();
		bRead.close();
		
	}
	public static void main(String[] args) throws Exception{
//		add_mass("/Users/yandj/Programming/Projects/SpectraPrediction/data/CASMI_Compounds/CASMI2016_Cat2and3_Challenge.txt", 
//				"/Users/yandj/Programming/Projects/SpectraPrediction/data/CASMI_Compounds/CASMI2016_Cat2and3_Challenge_with_mass.txt");
        
        
        
        
	}
}
