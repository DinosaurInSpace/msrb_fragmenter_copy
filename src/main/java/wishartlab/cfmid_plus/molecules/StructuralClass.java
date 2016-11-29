/**
 * 
 */
package wishartlab.cfmid_plus.molecules;

/**
 * @author yandj
 *
 */

import java.util.LinkedHashMap;

public class StructuralClass {
	public enum ClassName{
		PHOSPHATIDYLCHOLINES, _1_LYSOPHOSPHATIDYLCHOLINES, _2_LYSOPHOSPHATIDYLCHOLINES, PLASMENYLPHOSPHATIDYLCHOLINES,
		PHOSPHATIDYLETHANOLAMINES, _1_LYSOPHOSPHATIDYLETHANOLAMINES, _1_LYSOPLASMENYL_ETHANOLAMINES, _1_LYSOPLASMANYL_ETHANOLAMINES,
		_2_LYSOPHOSPHATIDYLETHANOLAMINES, PLASMENYLPHOSPHATIDYLETHANOLAMINES,
		PHOSPHATIDYLSERINES, SPHINGOMYELINS, PHOSPHATIDIC_ACIDS, _1_LYSOPHOSPHATIDIC_ACIDS, _2_LYSOPHOSPHATIDIC_ACIDS,
		PHOSPHATIDYLINOSITOLS, PHOSPHATIDYLGLYCEROLS, CARDIOLIPINS, CERAMIDE_1_PHOSPHATES,
		N_ACYL_SPHINGOSINES, SULFATIDES, GANGLIOSIDES, CHOLESTERYL_ESTERS,
		_1_MONOACYLGLYCEROLS, _2_MONOACYLGLYCEROLS, DIACYLGLYCEROLS, TRIACYLGLYCEROLS, MONOGALACTOSYLDIACYLGLYCEROLS,
		DIGALACTOSYLDIACYLGLYCEROLS, SULFOQUINOVOSYLDIACYLGLYCEROLS,
		DIACYLATED_PHOSPHATIDYLINOSITOL_MONOMANNOSIDES, DIACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES,
		TRIACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES, TETRACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES, 
		DIPHOSPHORYLATED_HEXAACYL_LIPID_A, ACYL_CARNITINES, ACYL_CHOLINES, NIL, CHOLESTERYL_ESTERS_LANOSTERYL_OR_ZYMOSTERYL
	}
	
	public static final LinkedHashMap<ClassName, String>backbones;
	static{
		backbones = new LinkedHashMap<ClassName, String>();
		backbones.put(ClassName.PHOSPHATIDYLCHOLINES,"[#6;A;X4][#6;X3](=[O;X1])-[#8][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H2X4][N;X4+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#8;X2]-[#6;X3]([#6;A;X4])=[O;X1]");
		backbones.put(ClassName._1_LYSOPHOSPHATIDYLCHOLINES,"[#6;A;X4;H1,H2,H3][#6](!@=[O;X1])!@-[#8]!@-[#6;A;H2X4]!@-[#6;A;H1X4]([#8;A;H1X2])!@-[#6;A;H2X4]!@-[#8;X2][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[N;X4+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4]");
		backbones.put(ClassName._2_LYSOPHOSPHATIDYLCHOLINES,"[#6;A;X4;H1,H2,H3][#6](=[O;X1])-[#8;X2][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H2X4][#8;X2]P([#8;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H2X4][N;X4+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4]");
		backbones.put(ClassName.PLASMENYLPHOSPHATIDYLCHOLINES,"[#6;H3X4][N+]([#6;H3X4])([#6;H3X4])[#6;H2X4]-[#6;H2X4]-[#8;X2]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2]-[#6;H2X4]-[#6;H1X4](-[#6;H2X4]-[#8]-[#6;H1X3]=[#6;X3])-[#8;X2]-[#6]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.PHOSPHATIDYLETHANOLAMINES,"[#6;A;X4][#6;R0](=[O;X1R0])-[#8;X2R0][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2R0]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2R0][#6;A;H2X4][#6;A;H2X4][#7;A;H2X3])[#8;X2R0]-[#6;R0]([#6;A;X4])=[O;R0]");
		backbones.put(ClassName._1_LYSOPHOSPHATIDYLETHANOLAMINES,"[#6;X4]-[#6;R0](=[O;X1R0])-[#8;X2R0][#6;A;H2X4][#6;A;H1X4]([#8;H1X2R0])[#6;A;H2X4][#8;X2R0]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2R0][#6;A;H2X4][#6;A;H2X4][#7;A;X3H2,X4H3+]");
		
		backbones.put(ClassName._1_LYSOPLASMENYL_ETHANOLAMINES,"[H]-[#6;X3]=[#6;R0](-[H])-[#8;X2R0][#6;A;H2X4][#6;A;H1X4]([#8;H1X2R0])[#6;A;H2X4][#8;X2R0]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2R0][#6;A;H2X4][#6;A;H2X4][#7;A;X3H2,X4H3+]");

		backbones.put(ClassName._1_LYSOPLASMANYL_ETHANOLAMINES,"[H][#6;X4]([H])[C;R0]([H])([H])[#8;X2R0][#6;A;H2X4][#6;A;H1X4]([#8;H1X2R0])[#6;A;H2X4][#8;X2R0]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2R0][#6;A;H2X4][#6;A;H2X4][#7;A;X3H2,X4H3+]");
		
		backbones.put(ClassName._2_LYSOPHOSPHATIDYLETHANOLAMINES,"[#6;X4]-[#6;R0](=[O;X1R0])-[#8;X2R0][#6;A;H1X4]([#6;A;H2X4][#8;H1X2R0])[#6;A;H2X4][#8;X2R0]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2R0][#6;A;H2X4][#6;A;H2X4][#7;A;X3H2,X4H3+]");
		backbones.put(ClassName.PLASMENYLPHOSPHATIDYLETHANOLAMINES,"[#7;A;H2X3][#6;H2X4]-[#6;H2X4]-[#8;X2]P([#8;A;X2H1,X1-])(=[O;X1])[#8;X2]-[#6;H2X4]-[#6;H1X4](-[#6;H2X4]-[#8]-[#6;H1X3]=[#6;X3]-[#6,#1;A])-[#8;X2]-[#6]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.PHOSPHATIDYLSERINES,"[#7;A;H2X3][#6;H1X4]([#6;A;H2X4][#8]P([#8;A;X2H1,X1-])(=O)[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2R0]-[#6;R0]([#6,#1;A])=[O;R0])[#8;R0]-[#6;R0]([#6,#1;A])=[O;R0])-[#6]([#8;A;X2H1,X1-])=O");

		// orginial from ClassyFire
		backbones.put(ClassName.SPHINGOMYELINS,"[$([#6;A;H2X4][#6;A;H2X4][#6;R0](=O)[#7;A;H1X3][#6;A;H1X4]([#6;A;H2X4][#8;R0]P([#8;X2H1,X1-])(=O)[#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6;H1X3R0]=[#6;H1X3R0][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H3X4]),$([#6;A;H2X4][#6;A;H2X4][#6;R0](=O)[#7;A;H1X3][#6;A;H1X4]([#6;A;H2X4][#8;R0]P([#8;X2H1,X1-])(=O)[#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H3X4])]");
		// Updated SMARTS
		// backbones.put(ClassName.SPHINGOMYELINS,"[H][#6;R0]([#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H3X4])=[#6;R0]([H])-[#6;R0](-[#8])[C;R0]([H])([#6;A;H2X4][#8;R0]P([!#1!#6;$(OH),$([O-])])(=O)[#8][#6;A;H2X4][#6;A;H2X4][N+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4])[#7]-[#6;R0](-[#6;R0])=O");
		
		backbones.put(ClassName.PHOSPHATIDIC_ACIDS,"[#8;A;X2H1,X1-][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;H1X4R0]([#6;A;H2X4][#8;A;X2][#6;R0]([#6,#1;A])=[O;X1])[#8;A;X2][#6;R0]([#6,#1;A])=[O;X1]");
 		backbones.put(ClassName._1_LYSOPHOSPHATIDIC_ACIDS,"[#6;X4]-[#6;R0](=[O;X1])[#8;A;X2][#6;A;H2X4][#6;H1X4R0]([#8;A;H1X2])[#6;A;H2X4][#8;X2][P;X4]([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1]");
		backbones.put(ClassName._2_LYSOPHOSPHATIDIC_ACIDS,"[#8;A;H1X2][#6;A;H2X4][#6;H1X4R0]([#6;A;H2X4][#8;X2][P;X4]([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1])[#8;A;X2][#6;R0]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.PHOSPHATIDYLINOSITOLS,"[#8;A;H1X2][#6;A;H1X4]1[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;X2][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]1[#8;A;H1X2]");
		backbones.put(ClassName.PHOSPHATIDYLGLYCEROLS,"[#8;A;H1X2][#6;A;H2X4][#6;A;H1X4]([#8;A;H1X2])[#6;A;H2X4][#8;X2][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.CARDIOLIPINS,"[#6;A;X4;H1,H2,H3][#6;X3](=O)-[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]P([#8;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#8;A;H1X2])[#6;A;H2X4][#8;X2]P([#8;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]-[#6;X3]([#6;A;X4;H1,H2,H3])=[O;X1])[#8;X2]-[#6;X3]([#6;A;X4;H1,H2,H3])=O)[#8;X2]-[#6;X3]([#6;A;X4;H1,H2,H3])=[O;X1]");
		backbones.put(ClassName.CERAMIDE_1_PHOSPHATES,"[#6;A;X4;H2,H3][#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;A;H1X3]=[#6;A;H1X3]-[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#6;H2X4]-[#8;X2][P;X4]([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1])[#7;A;H1X3][#6;X3]([#6,#1;A])=[O;X1]");		
//		backbones.put(ClassName.N_ACYL_SPHINGOSINES,"");
//		backbones.put(ClassName.SULFATIDES,"[#6;A;X4;H2,H3][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X3]=[#6;A;H1X3]-[#6;A@@H;H1X4]([#8;A;H1X2])[#6;A@H;H1X4]([#6;A;H2X4][#8][#6;A@@H;H1X4]1[#8][#6;A@H;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A@H;H1X4]([#8;A;H1X2])[#6;A@H;H1X4]([#8]S([#8;A;X2H1,X1-])(=O)=O)[#6;A@H;H1X4]1[#8;A;H1X2])[#7;A;H1X3][#6]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.SULFATIDES,"[#6;A;X4;H2,H3][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X3]=[#6;A;H1X3][#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#6;A;H2X4][#8][#6;A;H1X4]1[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;X2][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])[#6;A;H1X4]1[#8;A;H1X2])[#7;A;H1X3][#6]([#6,#1;A])=[O;X1]");
		
//		backbones.put(ClassName.GANGLIOSIDES,"[H][C@;X4]1([#8][C@@;X4]([#6;H2X4]-[#6@H;H1X4]([#8;A;H1X2])-[#6@H;H1X4]1-[#7;H1X3]-[#6](-[#6;H3X4])=O)([#8]-[#6@H;H1X4]-1-[#6@@H;H1X4]([#8;A;H1X2])-[#6@@H;H1X4](-[#6;H2X4]-[#8;H1X2])-[#8]-[#6@@H;H1X4](-[#8]-[#6;H2X4]-[#6@H;H1X4](-[#7;H1X3]-[#6]([#6,#1;A])=[O;X1])-[#6@H;H1X4]([#8;A;H1X2])\\[#6;H1X3]=[#6;H1X3]\\[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H1X4][#6;A;X4;H2,H3])-[#6@@H;H1X4]-1-[#8;H1X2])[#6](-[#8;H1X2])=[O;X1])[#6@H](-[#8;H1X2])-[#6@H](-[#8;H1X2])-[#6;H2X4]-[#8;H1X2]");
		
		// GM2, GM3, GM4
		backbones.put(ClassName.GANGLIOSIDES,"[$([#6;A;X4][#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;A;H1X3]=[#6;A;H1X3]-[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#6]-[#8][#6;A;H1X4]1[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8][#6;A;H1X4]2[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8][#6;A;H1X4]3[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]3[#7;A;H1X3][#6]([#6;A;H3X4])=O)[#6;A;H1X4]([#8]C3([#6;A;H2X4][#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#7;A;H1X3][#6](-[#6])=O)[#6;A;H1X4]([#8]3)[#6]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H2X4][#8;A;H1X2])[#6]([#8;A;X2H1,X1-])=O)[#6;A;H1X4]2[#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]1[#8;A;H1X2])[#7;A;H1X3][#6]([#6,#1;A])=O),$([#6;A;X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X3]=[#6;A;H1X3][#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#6;A;H2X4][#8][#6;A;H1X4]1[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8][#6;A;H1X4]2[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8]C3([#6;A;H2X4][#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#7;A;H1X3][#6](-[#6])=[O;X1])[#6;A;H1X4]([#8]3)[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H2X4][#8;A;H1X2])[#6]([#8;A;X2H1,X1-])=[O;X1])[#6;A;H1X4]2[#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]1[#8;A;H1X2])[#7;A;H1X3][#6]([#6,#1;A])=[O;X1]),$([#6;A;X4][#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;H2X4]-[#6;A;H1X3]=[#6;A;H1X3]-[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#6;A;H2X4][#8][#6;A;H1X4]1[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8]C2([#6;A;H2X4][#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#7;A;H1X3][#6]([#6;A;H3X4])=O)[#6;A;H1X4]([#8]2)[#6]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H2X4][#8;A;H1X2])[#6]([#8;A;X2H1,X1-])=O)[#6;A;H1X4]1[#8;A;H1X2])[#7;A;H1X3][#6]([#6,#1;A])=O)]");
		
		backbones.put(ClassName.CHOLESTERYL_ESTERS,"[#6;A;H3X4][#6;A;H1X4]([#6;A;H3X4])[#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#6;A;H3X4])[#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]2[#6;A;H1X4]3[#6;A;H2X4][#6;A;H1X3]=[#6]4[#6;A;H2X4][#6]([#6;A;H2X4][#6;A;H2X4]C4([#6;A;H3X4])[#6;A;H1X4]3[#6;A;H2X4][#6;A;H2X4]C12[#6;A;H3X4])-[#8;X2]-[#6;X3](-[#6])=[O;X1]");
//		backbones.put(ClassName.CHOLESTERYL_ESTERS_LANOSTERYL_OR_ZYMOSTERYL,"[#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X3]=[#6;A;X3](-[#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4][#6;A;X4]2([*,#1;#1,CX4H3])[#6;A;X3]3=[#6;A;X3]([#6;A;H2X4][#6;A;H2X4]C12[#6;A;H3X4])C1([#6;A;H3X4])[#6;A;H2X4][#6;A;H2X4][#6;H1X4](-[#8;X2]-[#6;X3](-[#6])=[O;X1])[#6;A;X4]([*,#1;#1,CX4H3])([*,#1;#1,CX4H3])[#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4]3");
//		backbones.put(ClassName.CHOLESTERYL_ESTERS_LANOSTERYL_OR_ZYMOSTERYL,"[#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X3]=[#6;A;X3](-[#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4][#6;A;X4;$([CX4](-[CX4H3,H])(C)(C)C)]2[#6;A;X3]3=[#6;A;X3]([#6;A;H2X4][#6;A;H2X4]C12[#6;A;H3X4])C1([#6;A;H3X4])[#6;A;H2X4][#6;A;H2X4][#6;H1X4]([#6;A;X4;$([CX4](-[H])(-[H])(C)C),$([CX4](-[CX4H3])(-[CX4H3])(C)C)][#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4]3)-[#8;X2]-[#6;X3](-[#6])=[O;X1]");
//		backbones.put(ClassName.CHOLESTERYL_ESTERS_LANOSTERYL_OR_ZYMOSTERYL,"[#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X3]=[#6;A;X3](-[#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4][#6;CX4H3,$([CX4]-[CX4H3])]2[#6;A;X3]3=[#6;A;X3]([#6;A;H2X4][#6;A;H2X4]C12[#6;A;H3X4])C1([#6;A;H3X4])[#6;A;H2X4][#6;A;H2X4][#6;H1X4](-[#6][#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4]3)-[#8;X2]-[#6;X3](-[#6])=[O;X1]");
		backbones.put(ClassName.CHOLESTERYL_ESTERS_LANOSTERYL_OR_ZYMOSTERYL,"[#6;A;H3X4][#6;A;H1X4]([#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X3]=[#6;A;X3](-[#6;A;H3X4])[#6;A;H3X4])[#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4][#6]2[#6;A;X3]3=[#6;A;X3]([#6;A;H2X4][#6;A;H2X4]C12[#6;A;H3X4])C1([#6;A;H3X4])[#6;A;H2X4][#6;A;H2X4][#6;H1X4](-[#6][#6;A;H1X4]1[#6;A;H2X4][#6;A;H2X4]3)-[#8;X2]-[#6;X3](-[#6])=[O;X1]");
		
		
		backbones.put(ClassName._1_MONOACYLGLYCEROLS,"[#8;A;H1X2][#6;A;H2X4][#6;H1X4R0]([#8;A;H1X2])[#6;A;H2X4][#8;X2]-[#6;R0]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName._2_MONOACYLGLYCEROLS,"[#8;A;H1X2][#6;A;H2X4][#6;H1X4R0]([#6;A;H2X4][#8;A;H1X2])-[#8;X2]-[#6;R0]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.DIACYLGLYCEROLS,"[#8;A;H1X2][#6;A;H2X4][#6;H1X4R0]([#6;A;H2X4][#8;X2]-[#6;R0]([#6,#1;A])=[O;X1])[#8;A;X2][#6;R0]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.TRIACYLGLYCEROLS,"[#6,#1;A][#6;R0](=[O;X1])-[#8;X2][#6;A;H2X4][#6;H1X4R0]([#6;A;H2X4][#8;A;X2][#6;R0]([#6,#1;A])=[O;X1])[#8;A;X2][#6;R0]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.MONOGALACTOSYLDIACYLGLYCEROLS,"[#8;A;H1X2][#6;A;H2X4][#6;A@H;H1X4]1[#8;X2][#6;A@@H;H1X4]([#8;X2][#6;A;H2X4][#6;A@@H;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A@H;H1X4]1[#8;A;H1X2]");
		backbones.put(ClassName.DIGALACTOSYLDIACYLGLYCEROLS,"[#8;A;H1X2][#6;A;H2X4][#6;A@H;H1X4]1[#8][#6;A@H;H1X4]([#8;X2][#6;A;H2X4][#6;A@H;H1X4]2[#8][#6;A@@H;H1X4]([#8;X2][#6;A;H2X4][#6;A@@H;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A@H;H1X4]([#8;A;H1X2])[#6;A@@H;H1X4]([#8;A;H1X2])[#6;A@H;H1X4]2[#8;A;H1X2])[#6;A@H;H1X4]([#8;A;H1X2])[#6;A@@H;H1X4]([#8;A;H1X2])[#6;A@H;H1X4]1[#8;A;H1X2]");
		backbones.put(ClassName.SULFOQUINOVOSYLDIACYLGLYCEROLS,"[#8;A;H1X2][#6;A;H1X4]1[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#6;A;H2X4][S;X4]([#8;A;X2H1,X1-])(=[O;X1])=[O;X1])[#8][#6;A;H1X4]([#8;X2][#6;A;H2X4][#6]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])-[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]1[#8;A;H1X2]");
		backbones.put(ClassName.DIACYLATED_PHOSPHATIDYLINOSITOL_MONOMANNOSIDES,"[#8;A;H1X2][#6;A;H2X4][#6;A;H1X4]1[#8][#6;A;H1X4]([#8;X2][#6;A;H1X4]2[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;X4]2[#8;X2][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]1[#8;A;H1X2]");
		backbones.put(ClassName.DIACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES,"[#8;A;H1X2][#6;A;H2X4][#6;A;H1X4]1[#8][#6;A;H1X4]([#8;X2][#6;A;H1X4]2[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;X2][#6;A;H1X4]3[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]3[#8;A;H1X2])[#6;A;X4]2[#8;X2][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]1[#8;A;H1X2]");
		backbones.put(ClassName.TRIACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES,"[#8;A;H1X2][#6;A;H2X4][#6;A;H1X4]1[#8][#6;A;H1X4]([#8;X2][#6;A;H1X4]2[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;X2][#6;A;H1X4]3[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;X2][#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]3[#8;A;H1X2])[#6;A;X4]2[#8;X2][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]1[#8;A;H1X2]");
		backbones.put(ClassName.TETRACYLATED_PHOSPHATIDYLINOSITOL_DIMANNOSIDES,"[#8;A;H1X2][#6;A;H2X4][#6;A;H1X4]1[#8][#6;A;H1X4]([#8;X2][#6;A;H1X4]2[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;X2][#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;X2][#6;A;H1X4]3[#8][#6;A;H1X4]([#6;A;H2X4][#8;A;X2][#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]3[#8;A;H1X2])[#6;A;X4]2[#8;X2][P;X4]([#8;A;X2H1,X1-])(=[O;X1])[#8;X2][#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#8;X2]-[#6]([#6,#1;A])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]([#8;A;H1X2])[#6;A;H1X4]1[#8;A;H1X2]");
//		backbones.put(ClassName.DIPHOSPHORYLATED_HEXAACYL_LIPID_A,"");
		backbones.put(ClassName.ACYL_CARNITINES,"[#6;A;H3X4][N;X4+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H2X4][#6;A;H1X4]([#6;A;H2X4][#6]([#8;A;X1-,X2H1])=[O;X1])[#8;X2]-[#6]([#6,#1;A])=[O;X1]");
		backbones.put(ClassName.ACYL_CHOLINES,"[#6;A;H3X4]!@-[N+](!@-[#6;A;H3X4])(!@-[#6;A;H3X4])!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#8]!@-[#6]([#6,#1;A])!@=O");
//		backbones.put(ClassName.LYSOPHOSPHATIDYLETHANOLAMINES,"");
//		backbones.put(ClassName.LYSOPHOSPHATIDYLETHANOLAMINES,"");
//		backbones.put(ClassName.LYSOPHOSPHATIDYLETHANOLAMINES,"");
//		backbones.put(ClassName.LYSOPHOSPHATIDYLETHANOLAMINES,"");
//		backbones.put(ClassName.LYSOPHOSPHATIDYLETHANOLAMINES,"");
//		backbones.put(ClassName.LYSOPHOSPHATIDYLETHANOLAMINES,"");
//		backbones.put(ClassName.NIL,"");
		
	}
	
	
}
