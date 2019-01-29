package wishartlab.cfmid_plus.fragmentation;

import java.util.LinkedHashMap;

import wishartlab.cfmid_plus.molecules.StructuralClass.ClassName;

public class MSPRelativeAbundanceList {

//	public static LinkedHashMap<ClassName, LinkedHashMap<FragmentationCondition,MSPeakRelativeAbundance>> classSpecificRelativeAbundances;
	public static LinkedHashMap<ClassName, LinkedHashMap<String,MSPeakRelativeAbundance>> classSpecificRelativeAbundances;

		static{
//			classSpecificRelativeAbundances = new LinkedHashMap<ClassName, LinkedHashMap<FragmentationCondition,MSPeakRelativeAbundance>>();
			classSpecificRelativeAbundances = new LinkedHashMap<ClassName, LinkedHashMap<String,MSPeakRelativeAbundance>>();

			setUpRealtiveAbundanceLibrary();
		}
		
		private static void setUpRealtiveAbundanceLibrary(){
			/*
			 * 1-Monoacylglycerols
			 */
	
			LinkedHashMap<String,MSPeakRelativeAbundance> _1_mg = new LinkedHashMap<String,MSPeakRelativeAbundance>();
			LinkedHashMap<String, Double> _1_MG_RA_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _1_MG_RA_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _1_MG_RA_40 = new LinkedHashMap<String, Double>();			
			_1_MG_RA_10.put("[M+Li]+", 1.0);		
			_1_MG_RA_20.put("[M+Li]+", 1.0);
			_1_MG_RA_40.put("[M+Li]+", 1.0);
			_1_MG_RA_40.put("[M+Li]-sn1", 1.0);
			_1_MG_RA_40.put("[M+Li]-sn1-H2O", 1.0);
			
			
			FragmentationCondition _1_MG_10 = new FragmentationCondition("[M+Li]+", 10);
			FragmentationCondition _1_MG_20 = new FragmentationCondition("[M+Li]+", 20);
			FragmentationCondition _1_MG_40 = new FragmentationCondition("[M+Li]+", 40);
					
			_1_mg.put(
					"[M+Li]+_10",				
					new MSPeakRelativeAbundance(
							ClassName._1_MONOACYLGLYCEROLS, _1_MG_10, _1_MG_RA_10)
					);
			_1_mg.put(
					"[M+Li]+_20",				
					new MSPeakRelativeAbundance(
							ClassName._1_MONOACYLGLYCEROLS, _1_MG_20, _1_MG_RA_20)
					);			
			_1_mg.put(
					"[M+Li]+_40",				
					new MSPeakRelativeAbundance(
							ClassName._1_MONOACYLGLYCEROLS, _1_MG_40, _1_MG_RA_40)
					);			
			
			
			LinkedHashMap<String, Double> _1_MG_NH4_RA_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _1_MG_NH4_RA_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _1_MG_NH4_RA_40 = new LinkedHashMap<String, Double>();
			_1_MG_NH4_RA_10.put("[M+NH4]+", 1.0);
			_1_MG_NH4_RA_20.put("[M+NH4]+", 1.0);
			_1_MG_NH4_RA_40.put("[M+NH4]+", 1.0);
			_1_MG_NH4_RA_40.put("[M-NH4+H]+ (-18)", 1.0);
			_1_MG_NH4_RA_40.put("[M+NH4]-NH4-H2O", 1.0);
			_1_MG_NH4_RA_40.put("Sn1-FA-H", 1.0);
			_1_MG_NH4_RA_40.put("Sn1-H", 1.0);
			
			FragmentationCondition _1_MG_NH4_10 = new FragmentationCondition("[M+NH4]+", 10);
			FragmentationCondition _1_MG_NH4_20 = new FragmentationCondition("[M+NH4]+", 20);
			FragmentationCondition _1_MG_NH4_40 = new FragmentationCondition("[M+NH4]+", 40);
			
			_1_mg.put(
					"[M+NH4]+_10",				
					new MSPeakRelativeAbundance(
							ClassName._1_MONOACYLGLYCEROLS, _1_MG_NH4_10, _1_MG_NH4_RA_10)
					);
			_1_mg.put(
					"[M+NH4]+_20",				
					new MSPeakRelativeAbundance(
							ClassName._1_MONOACYLGLYCEROLS, _1_MG_NH4_20, _1_MG_NH4_RA_20)
					);
			_1_mg.put(
					"[M+NH4]+_40",				
					new MSPeakRelativeAbundance(
							ClassName._1_MONOACYLGLYCEROLS, _1_MG_NH4_40, _1_MG_NH4_RA_40)
					);

			classSpecificRelativeAbundances.put(ClassName._1_MONOACYLGLYCEROLS,
					_1_mg);
		
			/*
			 * 2-Monoacylglycerols
			 */
			LinkedHashMap<String,MSPeakRelativeAbundance> _2_mg = new LinkedHashMap<String,MSPeakRelativeAbundance>();
			
			LinkedHashMap<String, Double> _2_MG_H_RA_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _2_MG_H_RA_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _2_MG_H_RA_40 = new LinkedHashMap<String, Double>();			
			_2_MG_H_RA_10.put("[M+H]+", 1.0);
			_2_MG_H_RA_20.put("[M+H]+", 1.0);
			_2_MG_H_RA_20.put("Sn1 FA-H", 1.0);	
			_2_MG_H_RA_20.put("[M+H]-(H2O)2", 1.0);	
			_2_MG_H_RA_40.put("[M+H]+", 1.0);
			_2_MG_H_RA_40.put("[M-H2O+H]+", 1.0);
			_2_MG_H_RA_40.put("[M+H]-(H2O)2", 1.0);
			_2_MG_H_RA_40.put("Sn1 FA-H", 1.0);
			_2_MG_H_RA_40.put("Sn1-H", 1.0);
			FragmentationCondition _2_MG_H_10 = new FragmentationCondition("[M+H]+", 10);
			FragmentationCondition _2_MG_H_20 = new FragmentationCondition("[M+H]+", 20);
			FragmentationCondition _2_MG_H_40 = new FragmentationCondition("[M+H]+", 40);
			
			_2_mg.put(
					"[M+H]+_10",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_H_10, _2_MG_H_RA_10));
			_2_mg.put(
					"[M+H]+_20",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_H_20, _2_MG_H_RA_20));			
			_2_mg.put(
					"[M+H]+_40",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_H_40, _2_MG_H_RA_40));			
			
			
			LinkedHashMap<String, Double> _2_MG_NH4_RA_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _2_MG_NH4_RA_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _2_MG_NH4_RA_40 = new LinkedHashMap<String, Double>();
			_2_MG_NH4_RA_10.put("[M+NH4]+", 1.0);
			_2_MG_NH4_RA_20.put("[M+NH4]+", 1.0);
			_2_MG_NH4_RA_20.put("[M+NH4]-H2O", 1.0);
			_2_MG_NH4_RA_40.put("[M+NH4]+", 1.0);
			_2_MG_NH4_RA_40.put("[M+NH4]-H2O", 1.0);
			_2_MG_NH4_RA_40.put("Sn1-FA + NH4+", 1.0);
//			_2_MG_NH4_RA_40.put("Sn1-FA", 1.0);			
			
			FragmentationCondition _2_MG_NH4_10 = new FragmentationCondition("[M+NH4]+", 10);
			FragmentationCondition _2_MG_NH4_20 = new FragmentationCondition("[M+NH4]+", 20);
			FragmentationCondition _2_MG_NH4_40 = new FragmentationCondition("[M+NH4]+", 40);
			
			_2_mg.put(
					"[M+NH4]+_10",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_NH4_10, _2_MG_NH4_RA_10));
			_2_mg.put(
					"[M+NH4]+_20",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_NH4_20, _2_MG_NH4_RA_20));
			_2_mg.put(
					"[M+NH4]+_40",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_NH4_40, _2_MG_NH4_RA_40));	
			
			
			LinkedHashMap<String, Double> _2_MG_NA_RA_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _2_MG_NA_RA_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> _2_MG_NA_RA_40 = new LinkedHashMap<String, Double>();	
			
			_2_MG_NA_RA_10.put("[M+Na]+", 1.0);
			_2_MG_NA_RA_20.put("[M+Na]+", 1.0);
			_2_MG_NA_RA_40.put("[M+Na]+", 1.0);
			_2_MG_NA_RA_40.put("Sn1 FA-H + [Na+]", 1.0);
			_2_MG_NA_RA_40.put("[M+Na]-(H2O)2-Na", 1.0);
			
			FragmentationCondition _2_MG_NA_10 = new FragmentationCondition("[M+Na]+", 10);
			FragmentationCondition _2_MG_NA_20 = new FragmentationCondition("[M+Na]+", 20);
			FragmentationCondition _2_MG_NA_40 = new FragmentationCondition("[M+Na]+", 40);
			_2_mg.put(
					"[M+Na]+_10",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_NA_10, _2_MG_NA_RA_10));
			_2_mg.put(
					"[M+Na]+_20",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_NA_20, _2_MG_NA_RA_20));
			_2_mg.put(
					"[M+Na]+_40",				
					new MSPeakRelativeAbundance(
							ClassName._2_MONOACYLGLYCEROLS, _2_MG_NA_40, _2_MG_NA_RA_40));
			
			classSpecificRelativeAbundances.put(ClassName._2_MONOACYLGLYCEROLS,
					_2_mg);
			
			/*
			 * Diacylglycerols
			 */
			
//			dg_frag_m_nh4
			LinkedHashMap<String,MSPeakRelativeAbundance> dg = new LinkedHashMap<String,MSPeakRelativeAbundance>();
			LinkedHashMap<String, Double> DG_NH4_RA_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> DG_NH4_RA_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> DG_NH4_RA_40 = new LinkedHashMap<String, Double>();
			
			
			DG_NH4_RA_10.put("[M+NH4]+", 100.0);
			DG_NH4_RA_20.put("[M+NH4]+", 15.0);
			DG_NH4_RA_20.put("[M+NH4]-NH4", 15.0);
			DG_NH4_RA_20.put("[M+NH4]-NH4-H2O (-17-18)", 60.0);
			DG_NH4_RA_20.put("[M+NH4]-NH4-sn1-H2O", 90.0);
			DG_NH4_RA_20.put("[M+NH4]-NH4-sn2-H2O", 90.0);
			DG_NH4_RA_40.put("[M+NH4]+", 60.0);
			DG_NH4_RA_40.put("[M+NH4]-NH4", 15.0);
			DG_NH4_RA_40.put("[M+NH4]-NH4-H2O (-17-18)", 15.0);
			DG_NH4_RA_40.put("[M+NH4]-NH4-sn1-H2O", 90.0);
			DG_NH4_RA_40.put("[M+NH4]-NH4-sn2-H2O", 90.0);
			
			dg.put("[M+NH4]+_10",				
					new MSPeakRelativeAbundance(
							ClassName._12_DIACYLGLYCEROLS, new FragmentationCondition("[M+NH4]+", 10), DG_NH4_RA_10)
					);	
			
			dg.put("[M+NH4]+_20",				
					new MSPeakRelativeAbundance(
							ClassName._12_DIACYLGLYCEROLS, new FragmentationCondition("[M+NH4]+", 20), DG_NH4_RA_20)
					);	
			
			dg.put("[M+NH4]+_40",				
					new MSPeakRelativeAbundance(
							ClassName._12_DIACYLGLYCEROLS, new FragmentationCondition("[M+NH4]+", 40), DG_NH4_RA_40)
					);
			
			classSpecificRelativeAbundances.put(ClassName._12_DIACYLGLYCEROLS, dg);

			
			
			
//			dg_frag_m_na
//			LinkedHashMap<String, Double> DG_NA_RA_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> DG_NA_RA_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> DG_NA_RA_40 = new LinkedHashMap<String, Double>();
//			
//			
//			DG_NA_RA_10.put("[M+Na]+", 100.0);
//			DG_NA_RA_20.put("[M+Na]+", 100.0);
//			DG_NA_RA_40.put("[M+Na]+", 100.0);
//			DG_NA_RA_40.put("[M+Na]-sn1-H2O", 100.0);
//			DG_NA_RA_40.put("[M+Na]-sn2-H2O", 100.0);
//			
//			dg.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._12_DIACYLGLYCEROLS, new FragmentationCondition("[M+Na]+", 10), DG_NA_RA_10)
//					);	
//			
//			dg.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._12_DIACYLGLYCEROLS, new FragmentationCondition("[M+Na]+", 20), DG_NA_RA_20)
//					);	
//			
//			dg.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._12_DIACYLGLYCEROLS, new FragmentationCondition("[M+Na]+", 40), DG_NA_RA_40)
//					);
			
			classSpecificRelativeAbundances.put(ClassName._12_DIACYLGLYCEROLS, dg);
			
			
			
			
			
			LinkedHashMap<String,MSPeakRelativeAbundance> tg = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> TG_NA_RA_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> TG_NA_RA_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> TG_NA_RA_40 = new LinkedHashMap<String, Double>();
//			
//			TG_NA_RA_10.put("[M+Na]+", 1.0);
//			TG_NA_RA_20.put("[M+Na]+", 1.0);
//			TG_NA_RA_40.put("[M+Na]+", 1.0);
//			TG_NA_RA_40.put("[M+Na]-sn1-H2O", 1.0);
//			TG_NA_RA_40.put("[M+Na]-sn2-H2O", 1.0);
//			TG_NA_RA_40.put("[M+Na]-sn3-H2O", 1.0);
//			
//			tg.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+Na]+", 10), TG_NA_RA_10));				
//			tg.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+Na]+", 20), TG_NA_RA_20));				
//			tg.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+Na]+", 40), TG_NA_RA_40));		
			
			LinkedHashMap<String, Double> TG_NH4_RA_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> TG_NH4_RA_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> TG_NH4_RA_40 = new LinkedHashMap<String, Double>();
			
			TG_NH4_RA_10.put("[M+NH4]+", 100.0);
			TG_NH4_RA_20.put("[M+NH4]+", 100.0);
			TG_NH4_RA_40.put("[M+NH4]+", 15.0);
			TG_NH4_RA_40.put("[M+NH4]-17", 60.0);
			TG_NH4_RA_40.put("[M+NH4]-NH4-sn1-H2O", 90.0);
			TG_NH4_RA_40.put("[M+NH4]-NH4-sn2-H2O", 90.0);
			TG_NH4_RA_40.put("[M+NH4]-NH4-sn3-H2O", 90.0);
			
			tg.put("[M+NH4]+_10",				
					new MSPeakRelativeAbundance(
							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+NH4]+", 10), TG_NH4_RA_10));				
			tg.put("[M+NH4]+_20",				
					new MSPeakRelativeAbundance(
							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+NH4]+", 20), TG_NH4_RA_20));				
			tg.put("[M+NH4]+_40",				
					new MSPeakRelativeAbundance(
							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+NH4]+", 40), TG_NH4_RA_40));				

			
			
//			LinkedHashMap<String, Double> TG_LI_RA_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> TG_LI_RA_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> TG_LI_RA_40 = new LinkedHashMap<String, Double>();
//			
//			TG_LI_RA_10.put("[M+Li]+", 1.0);
//			TG_LI_RA_20.put("[M+Li]+", 1.0);
//			TG_LI_RA_40.put("[M+Li]+", 1.0);
//			TG_LI_RA_40.put("[M+Li]-sn1-H2O", 1.0);
//			TG_LI_RA_40.put("[M+Li]-sn2-H2O", 1.0);
//			TG_LI_RA_40.put("[M+Li]-sn3-H2O", 1.0);
//			
//			tg.put("[M+Li]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+Li]+", 10), TG_LI_RA_10));				
//			tg.put("[M+Li]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+Li]+", 20), TG_LI_RA_20));				
//			tg.put("[M+Li]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.TRIACYLGLYCEROLS, new FragmentationCondition("[M+Li]+", 40), TG_LI_RA_40));				

			
			
			classSpecificRelativeAbundances.put(ClassName.TRIACYLGLYCEROLS,
					tg);
			
			
			
			
			
			/*
			 * Phosphatidic acids ([M+H]+, [M+Na]+, [M-H]-)
			 * 
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=671&LM_ID=LMGP10010023&TRACK_ID=254
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=699&LM_ID=LMGP10010036&TRACK_ID=256
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=723&LM_ID=LMGP10010024&TRACK_ID=260
			 */
			
//			LinkedHashMap<String,MSPeakRelativeAbundance> pa= new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			
//			LinkedHashMap<String, Double> PA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();			
//			PA_RelativeAbundances_10.put("[M+H]+", 1.0);
//			PA_RelativeAbundances_10.put("[M-H2O+H]+", 1.0);
//			PA_RelativeAbundances_10.put("[M+H]-H3PO4 (-98)", 1.0);
//			PA_RelativeAbundances_20.put("[M+H]+", 1.0);
//			PA_RelativeAbundances_20.put("[M-H2O+H]+", 1.0);
//			PA_RelativeAbundances_20.put("[M+H]-H3PO4 (-98)", 1.0);
//			PA_RelativeAbundances_40.put("[M+H]+", 1.0);
//			PA_RelativeAbundances_40.put("[M-H2O+H]+", 1.0);
//			PA_RelativeAbundances_40.put("[M+H]-H3PO4 (-98)", 1.0);
//			PA_RelativeAbundances_40.put("[M+H]-sn1-H2O", 1.0);
//			PA_RelativeAbundances_40.put("[M+H]-sn2-H2O", 1.0);			
//			
//			
//			pa.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+H]+", 10), PA_RelativeAbundances_10)
//					);			
//			pa.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+H]+", 20), PA_RelativeAbundances_20)
//					);
//			pa.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+H]+", 40), PA_RelativeAbundances_40)
//					);
//					
//		
//			LinkedHashMap<String, Double> PA_Na_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PA_Na_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PA_Na_RelativeAbundances_40 = new LinkedHashMap<String, Double>();		
//			PA_Na_RelativeAbundances_10.put("[M+Na]+", 1.0);
//			PA_Na_RelativeAbundances_20.put("[M+Na]+", 1.0);
//			PA_Na_RelativeAbundances_20.put("[M+Na]-H3PO4", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]+", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]-H2O", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]-H3PO4", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]-H2PO4Na", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]-sn1", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]-sn2", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]-sn1-H2O", 1.0);
//			PA_Na_RelativeAbundances_40.put("[M+Na]-sn2-H2O", 1.0);
//		
//			pa.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+Na]+", 10), PA_Na_RelativeAbundances_10)
//					);
//			pa.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+Na]+", 20), PA_Na_RelativeAbundances_20)
//					);		
//			pa.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+Na]+", 40), PA_Na_RelativeAbundances_40)
//					);		
//		
//		
//			LinkedHashMap<String, Double> PA_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PA_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PA_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PA_H_RelativeAbundances_10.put("[M-H]-", 1.0);
//			PA_H_RelativeAbundances_20.put("[M-H]-", 1.0);
//			PA_H_RelativeAbundances_20.put("[M-H]-sn1", 1.0);
//			PA_H_RelativeAbundances_20.put("[M-H]-sn2", 1.0);
//			PA_H_RelativeAbundances_20.put("[M-H]-sn1-H2O", 1.0);
//			PA_H_RelativeAbundances_20.put("[M-H]-sn2-H2O", 1.0);
//			PA_H_RelativeAbundances_20.put("Sn1 FA-H", 1.0);
//			PA_H_RelativeAbundances_20.put("Sn2 FA-H", 1.0);
//			PA_H_RelativeAbundances_20.put("Ion fragment C3H6O4[O‐]P", 1.0);
//			PA_H_RelativeAbundances_20.put("Oxophosphinate", 1.0);	
//			PA_H_RelativeAbundances_40.put("[M-H]-", 1.0);
//			PA_H_RelativeAbundances_40.put("[M-H]-sn1", 1.0);
//			PA_H_RelativeAbundances_40.put("[M-H]-sn2", 1.0);
//			PA_H_RelativeAbundances_40.put("[M-H]-sn1-H2O", 1.0);
//			PA_H_RelativeAbundances_40.put("[M-H]-sn2-H2O", 1.0);
//			PA_H_RelativeAbundances_40.put("Sn1 FA-H", 1.0);
//			PA_H_RelativeAbundances_40.put("Sn2 FA-H", 1.0);
//			PA_H_RelativeAbundances_40.put("Ion fragment C3H6O4[O‐]P", 1.0);
//			PA_H_RelativeAbundances_40.put("Oxophosphinate", 1.0);
//						
//			pa.put("[M-H]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M-H]-", 10), PA_H_RelativeAbundances_10));
//			pa.put("[M-H]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M-H]-", 20), PA_H_RelativeAbundances_20));		
//			pa.put("[M-H]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDIC_ACIDS, new FragmentationCondition("[M-H]-", 40), PA_H_RelativeAbundances_40));
//
//			classSpecificRelativeAbundances.put(ClassName.PHOSPHATIDIC_ACIDS, pa);
//			
			
			
			
			/*
			 * Lysophosphatidic acids ([M+H]+, [M+Na]+, [M-H]-)
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=435&LM_ID=LMGP10050005&TRACK_ID=263
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=437&LM_ID=LMGP10050008&TRACK_ID=264
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=457&LM_ID=LMGP10050013&TRACK_ID=515
			 * 
			 */
			LinkedHashMap<String,MSPeakRelativeAbundance> _1_lpa= new LinkedHashMap<String,MSPeakRelativeAbundance>();	
			
//			LinkedHashMap<String, Double> _1_LPA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			_1_LPA_RelativeAbundances_10.put("[M+H]+", 1.0);
//			_1_LPA_RelativeAbundances_20.put("[M+H]+", 1.0);
//			_1_LPA_RelativeAbundances_40.put("[M+H]+", 1.0);
//			_1_LPA_RelativeAbundances_40.put("[M-H2O+H]+", 1.0);
//			_1_LPA_RelativeAbundances_40.put("Sn1-H", 1.0);
//			_1_LPA_RelativeAbundances_40.put("[M+H]-sn1", 1.0);
//			_1_LPA_RelativeAbundances_40.put("Ion fragment C3H8O4[O+]P", 1.0);
//			_1_lpa.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+H]+", 10), _1_LPA_RelativeAbundances_10));
//			_1_lpa.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+H]+", 20), _1_LPA_RelativeAbundances_20));
//			_1_lpa.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDIC_ACIDS, new FragmentationCondition("[M+H]+", 40), _1_LPA_RelativeAbundances_40));
//
//			
//			
//			LinkedHashMap<String, Double> _1_LPA_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPA_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPA_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			_1_LPA_H_RelativeAbundances_10.put("[M-H]-", 1.0);
//			_1_LPA_H_RelativeAbundances_20.put("[M-H]-", 1.0);
//			_1_LPA_H_RelativeAbundances_40.put("[M-H]-", 1.0);
//			_1_LPA_H_RelativeAbundances_40.put("[M-H]-H2O (-18)", 1.0);
//			_1_LPA_H_RelativeAbundances_40.put("[M-H]-sn1", 1.0);
//			_1_LPA_H_RelativeAbundances_40.put("Sn1 FA-H", 1.0);
//			_1_LPA_H_RelativeAbundances_40.put("Oxophosphinate", 1.0);
//			_1_LPA_H_RelativeAbundances_40.put("[M-H]-sn1-H2O", 1.0);
//			_1_LPA_H_RelativeAbundances_40.put("Phosphate", 1.0);
//			_1_lpa.put("[M-H]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDIC_ACIDS, new FragmentationCondition("[M-H]-", 10), _1_LPA_H_RelativeAbundances_10));
//			_1_lpa.put("[M-H]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDIC_ACIDS, new FragmentationCondition("[M-H]-", 20), _1_LPA_H_RelativeAbundances_20));
//			_1_lpa.put("[M-H]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDIC_ACIDS, new FragmentationCondition("[M-H]-", 40), _1_LPA_H_RelativeAbundances_40));

			
			
			
			classSpecificRelativeAbundances.put(ClassName._1_LYSOPHOSPHATIDIC_ACIDS, _1_lpa);
			
			
			/*
			 * Phosphatidylethanolamines ([M+H]+, [M+Na]+, [M-H]-)
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=718&LM_ID=LMGP02010009&TRACK_ID=78
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=718&LM_ID=LMGP02010009&TRACK_ID=307
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=740&LM_ID=LMGP02010096&TRACK_ID=309
			 * 
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=716&LM_ID=LMGP02010009&TRACK_ID=79
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=716&LM_ID=LMGP02010009&TRACK_ID=306
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=738&LM_ID=LMGP02010096&TRACK_ID=308
			 */		
			
			LinkedHashMap<String,MSPeakRelativeAbundance> pe= new LinkedHashMap<String,MSPeakRelativeAbundance>();
			LinkedHashMap<String, Double> PE_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PE_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PE_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
			
			PE_RelativeAbundances_10.put("[M+H]+", 90.0);
			PE_RelativeAbundances_10.put("[M-H2O+H]+", 15.0);
			PE_RelativeAbundances_10.put("[M+H]-C2H8NO4P (-141)", 15.0);
			
			PE_RelativeAbundances_20.put("[M+H]+", 60.0);
			PE_RelativeAbundances_20.put("[M-H2O+H]+", 15.0);
			PE_RelativeAbundances_20.put("[M+H]-C2H8NO4P (-141)", 90.0);
			PE_RelativeAbundances_20.put("[M+H]-sn1", 15.0);
			PE_RelativeAbundances_20.put("[M+H]-sn2", 15.0);
			PE_RelativeAbundances_20.put("[M+H]-sn1-H2O", 15.0);
			PE_RelativeAbundances_20.put("[M+H]-sn2-H2O", 15.0);	
			PE_RelativeAbundances_20.put("[M+H]-sn1-C2H8NO4P (-141)", 15.0);
			PE_RelativeAbundances_20.put("[M+H]-sn2-C2H8NO4P (-141)", 15.0);
			
			PE_RelativeAbundances_40.put("[M+H]+", 15.0);
			PE_RelativeAbundances_40.put("[M-H2O+H]+", 15.0);
			PE_RelativeAbundances_40.put("[M+H]-C2H8NO4P (-141)", 90.0);
			PE_RelativeAbundances_40.put("[M+H]-sn1", 15.0);
			PE_RelativeAbundances_40.put("[M+H]-sn2", 15.0);
			PE_RelativeAbundances_40.put("[M+H]-sn1-H2O", 15.0);
			PE_RelativeAbundances_40.put("[M+H]-sn2-H2O", 15.0);
			PE_RelativeAbundances_40.put("[M+H]-sn1-C2H8NO4P (-141)", 15.0);
			PE_RelativeAbundances_40.put("[M+H]-sn2-C2H8NO4P (-141)", 15.0);
			
			pe.put("[M+H]+_10",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M+H]+", 10), PE_RelativeAbundances_10));
			pe.put("[M+H]+_20",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M+H]+", 20), PE_RelativeAbundances_20));		
			pe.put("[M+H]+_40",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M+H]+", 40), PE_RelativeAbundances_40));

			
			
//			LinkedHashMap<String, Double> PE_NA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PE_NA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PE_NA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			
//			PE_NA_RelativeAbundances_10.put("[M+Na]+", 1.0);
//			PE_NA_RelativeAbundances_10.put("[M+Na]-C2H8NO4P (-141)", 1.0);			
//			PE_NA_RelativeAbundances_20.put("[M+Na]+", 1.0);
//			PE_NA_RelativeAbundances_20.put("[M+Na]-C2H5N (-43)", 1.0);
//			PE_NA_RelativeAbundances_20.put("[M+Na]-C2H8NO4P (-141)", 1.0);		
//			PE_NA_RelativeAbundances_40.put("[M+Na]+", 1.0);
//			PE_NA_RelativeAbundances_40.put("[M+Na]-C2H5N (-43)", 1.0);
//			PE_NA_RelativeAbundances_40.put("[M+Na]-C2H8NO4P (-141)", 1.0);
//			PE_NA_RelativeAbundances_40.put("[M+Na]-43-sn1-H2O", 1.0);
//			PE_NA_RelativeAbundances_40.put("[M+Na]-43-sn2-H2O", 1.0);
//			PE_NA_RelativeAbundances_40.put("C2H8NO4P (-141) + Na", 1.0);
//			pe.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M+Na]+", 10), PE_NA_RelativeAbundances_10));
//			pe.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M+Na]+", 20), PE_NA_RelativeAbundances_20));
//			pe.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M+Na]+", 40), PE_NA_RelativeAbundances_40));
			
			
			LinkedHashMap<String, Double> PE_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PE_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PE_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
						
			PE_H_RelativeAbundances_10.put("[M-H]-", 90.0);
			PE_H_RelativeAbundances_10.put("Sn1 FA-H", 15.0);
			PE_H_RelativeAbundances_10.put("Sn2 FA-H", 15.0);
			
			PE_H_RelativeAbundances_20.put("[M-H]-", 90.0);
			PE_H_RelativeAbundances_20.put("Sn1 FA-H", 15.0);
			PE_H_RelativeAbundances_20.put("Sn2 FA-H", 15.0);
			
			PE_H_RelativeAbundances_40.put("[M-H]-", 60.0);
			PE_H_RelativeAbundances_40.put("Sn1 FA-H", 90.0);
			PE_H_RelativeAbundances_40.put("Sn2 FA-H", 90.0);
			PE_H_RelativeAbundances_40.put("[M-H]-sn1", 15.0);
			PE_H_RelativeAbundances_40.put("[M-H]-sn2", 15.0);
			PE_H_RelativeAbundances_40.put("[M-H]-sn1-H2O", 15.0);
			PE_H_RelativeAbundances_40.put("[M-H]-sn2-H2O", 15.0);
			PE_H_RelativeAbundances_40.put("Ion fragment C3H6O4[O‐]P", 15.0);
			PE_H_RelativeAbundances_40.put("Ion fragment C5H11NO4[O‐]P", 15.0);
			
			pe.put("[M-H]-_10",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M-H]-_10", 10), PE_H_RelativeAbundances_10));
			pe.put("[M-H]-_20",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M-H]-_20", 20), PE_H_RelativeAbundances_20));
			pe.put("[M-H]-_40",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLETHANOLAMINES, new FragmentationCondition("[M-H]-_40", 40), PE_H_RelativeAbundances_40));
			
			
			
			classSpecificRelativeAbundances.put(ClassName.PHOSPHATIDYLETHANOLAMINES, pe);
			
			/*
			 * Phosphatidylcholines
			 */
			
			LinkedHashMap<String,MSPeakRelativeAbundance> pc= new LinkedHashMap<String,MSPeakRelativeAbundance>();
			
			LinkedHashMap<String, Double> PC_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PC_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PC_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
			PC_H_RelativeAbundances_10.put("[M+H]+", 100.0);
			PC_H_RelativeAbundances_20.put("[M+H]+", 90.0);
			PC_H_RelativeAbundances_20.put("Phosphocholine", 60.0);
			PC_H_RelativeAbundances_40.put("[M+H]+", 15.0);
			PC_H_RelativeAbundances_40.put("[M-H2O+H]", 15.0);
			PC_H_RelativeAbundances_40.put("[M+H]-C3H9N (-59)", 15.0);
			PC_H_RelativeAbundances_40.put("[M+H]-C5H14NO4P (-183)", 15.0);
			PC_H_RelativeAbundances_40.put("[M+H]-sn1", 15.0);
			PC_H_RelativeAbundances_40.put("[M+H]-sn2", 15.0);
			PC_H_RelativeAbundances_40.put("[M+H]-sn1-H2O", 15.0);
			PC_H_RelativeAbundances_40.put("[M+H]-sn2-H2O", 15.0);
			PC_H_RelativeAbundances_40.put("Phosphocholine", 90.0);
			PC_H_RelativeAbundances_40.put("Cyclic 1,2-cyclic phosphate diester", 15.0);
			PC_H_RelativeAbundances_40.put("Choline", 15.0);
			PC_H_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 15.0);

			
			
			pc.put("[M+H]+_10",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_10", 10), PC_H_RelativeAbundances_10));
			pc.put("[M+H]+_20",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_20", 20), PC_H_RelativeAbundances_20));
			pc.put("[M+H]+_40",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_40", 40), PC_H_RelativeAbundances_40));
			
			
//			LinkedHashMap<String, Double> PC_NA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PC_NA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PC_NA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PC_NA_RelativeAbundances_10.put("[M+Na]+", 100.0);
//			PC_NA_RelativeAbundances_20.put("[M+Na]+", 100.0);
//			PC_NA_RelativeAbundances_20.put("[M+Na]-C3H9N (-59)", 33.3);
//			PC_NA_RelativeAbundances_40.put("[M+Na]+", 100.0);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-C3H9N (-59)", 100.0);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-C5H14NO4P (-183)", 66.6);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-sn1", 33.3);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-sn2", 33.3);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-sn1-H2O", 33.3);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-sn2-H2O", 33.3);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-sn1-H2O-C3H9N", 33.3);
//			PC_NA_RelativeAbundances_40.put("[M+Na]-sn2-H2O-C3H9N", 33.3);
//			PC_NA_RelativeAbundances_40.put("Sodiated cyclic 1,2-cyclic phosphate diester", 33.3);
//			
//			pc.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_10", 10), PC_NA_RelativeAbundances_10));
//			pc.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_20", 20), PC_NA_RelativeAbundances_20));
//			pc.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_40", 40), PC_NA_RelativeAbundances_40));
//			
//			
//			
//			LinkedHashMap<String, Double> PC_LI_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PC_LI_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PC_LI_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PC_LI_RelativeAbundances_10.put("[M+Li]+", 100.0);
//			PC_LI_RelativeAbundances_20.put("[M+Li]+", 100.0);
//			PC_LI_RelativeAbundances_20.put("[M+Li]-C3H9N (-59)", 33.3);
//			PC_LI_RelativeAbundances_40.put("[M+Li]+", 66.6);
//			PC_LI_RelativeAbundances_40.put("[M+Li]-C3H9N (-59)", 33.3);
//			PC_LI_RelativeAbundances_40.put("[M+Li]-C5H14NO4P (-183)", 66.6);
//			PC_LI_RelativeAbundances_40.put("[M+H]-C5H14NO4P (-189)", 66.6);
//			
//			PC_LI_RelativeAbundances_40.put("[M+Li]-sn1", 100.0);
//			PC_LI_RelativeAbundances_40.put("[M+Li]-sn2", 100.0);
//			PC_LI_RelativeAbundances_40.put("[M+Li]-sn1-H2O", 33.3);
//			PC_LI_RelativeAbundances_40.put("[M+Li]-sn2-H2O", 33.3);
//			PC_LI_RelativeAbundances_40.put("[M+Li]-sn1-H2O-C3H9N", 33.3);
//			PC_LI_RelativeAbundances_40.put("[M+Li]-sn2-H2O-C3H9N", 33.3);
//			PC_LI_RelativeAbundances_40.put("Phosphocholine", 33.3);			
//			pc.put("[M+Li]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Li]+_10", 10), PC_LI_RelativeAbundances_10));
//			pc.put("[M+Li]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Li]+_20", 20), PC_LI_RelativeAbundances_20));
//			pc.put("[M+Li]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Li]+_40", 40), PC_LI_RelativeAbundances_40));		
			
			LinkedHashMap<String, Double> PC_CL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PC_CL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PC_CL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
			PC_CL_RelativeAbundances_10.put("[M+Cl]-", 90.0);
			PC_CL_RelativeAbundances_10.put("[M-H]-CH3", 15.0);
			
			PC_CL_RelativeAbundances_20.put("[M+Cl]-", 90.0);
			PC_CL_RelativeAbundances_20.put("[M-H]-CH3", 60.0);
			PC_CL_RelativeAbundances_20.put("Sn1 FA - H", 15.0);
			PC_CL_RelativeAbundances_20.put("Sn2 FA - H", 15.0);
			
			PC_CL_RelativeAbundances_20.put("[M-H]-sn1-CH3", 15.0);	
			PC_CL_RelativeAbundances_20.put("[M-H]-sn1-CH3-H2O", 15.0);
			PC_CL_RelativeAbundances_20.put("[M-H]-sn2-CH3", 15.0);
			PC_CL_RelativeAbundances_20.put("[M-H]-sn2-CH3-H2O", 15.0);			
			
			PC_CL_RelativeAbundances_40.put("[M+Cl]-", 60.0);
			PC_CL_RelativeAbundances_40.put("[M-H]-CH3", 90.0);

			PC_CL_RelativeAbundances_40.put("[M-H]-sn1-CH3", 15.0);	
			PC_CL_RelativeAbundances_40.put("[M-H]-sn1-CH3-H2O", 15.0);
			PC_CL_RelativeAbundances_40.put("[M-H]-sn2-CH3", 15.0);
			PC_CL_RelativeAbundances_40.put("[M-H]-sn2-CH3-H2O", 15.0);
			PC_CL_RelativeAbundances_40.put("Sn1 FA - H", 60.0);
			PC_CL_RelativeAbundances_40.put("Sn2 FA - H", 60.0);	
			PC_CL_RelativeAbundances_40.put("[M-H]-sn1-sn2-CH3-H2O", 15.0);
			
			pc.put("[M+Cl]-_10",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_10", 10), PC_CL_RelativeAbundances_10));
			pc.put("[M+Cl]-_20",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_20", 20), PC_CL_RelativeAbundances_20));
			pc.put("[M+Cl]-_40",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_40", 40), PC_CL_RelativeAbundances_40));			
			
			classSpecificRelativeAbundances.put(ClassName.PHOSPHATIDYLCHOLINES, pc);
			
			
			
			/*
			 * Lysophosphatidylcholines
			 */
//			LinkedHashMap<String,MSPeakRelativeAbundance> _1_lpc= new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> _1_LPC_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();		
//			_1_LPC_H_RelativeAbundances_10.put("[M+H]+", 1.0);
//			_1_LPC_H_RelativeAbundances_10.put("[M-H2O+H]+", 1.0);
//			_1_LPC_H_RelativeAbundances_10.put("Phosphocholine", 1.0);
//			_1_LPC_H_RelativeAbundances_10.put("Choline", 1.0);
//			_1_LPC_H_RelativeAbundances_20.put("[M+H]+", 1.0);
//			_1_LPC_H_RelativeAbundances_20.put("[M-H2O+H]+", 1.0);
//			_1_LPC_H_RelativeAbundances_20.put("Phosphocholine", 1.0);
//			_1_LPC_H_RelativeAbundances_20.put("Choline", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("[M+H]+", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("[M-H2O+H]+", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("Phosphocholine", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("Choline", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("[M+H]-sn1 (glycerophosphocholine)", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("[M+H]-sn1-H2O", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 1.0);
//			_1_LPC_H_RelativeAbundances_40.put("1,2-cyclic phosphate diester", 1.0);
//			_1_lpc.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_10", 10), _1_LPC_H_RelativeAbundances_10));
//			_1_lpc.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_20", 20), _1_LPC_H_RelativeAbundances_20));			
//			_1_lpc.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_40", 40), _1_LPC_H_RelativeAbundances_40));	
//			
//			LinkedHashMap<String, Double> _1_LPC_NA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_NA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_NA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();		
//			_1_LPC_NA_RelativeAbundances_10.put("[M+Na]+", 1.0);
//			_1_LPC_NA_RelativeAbundances_10.put("[M+Na]-C3H9N", 1.0);
//			_1_LPC_NA_RelativeAbundances_20.put("[M+Na]+", 1.0);
//			_1_LPC_NA_RelativeAbundances_20.put("[M+Na]-C3H9N", 1.0);
//			_1_LPC_NA_RelativeAbundances_20.put("[M+Na]-Na-C5H13NO4P", 1.0);			
//			_1_LPC_NA_RelativeAbundances_40.put("[M+Na]+", 1.0);
//			_1_LPC_NA_RelativeAbundances_40.put("[M+Na]-C3H9N", 1.0);
//			_1_LPC_NA_RelativeAbundances_40.put("[M+Na]-C5H14[N+]O", 1.0);
//			_1_LPC_NA_RelativeAbundances_40.put("[M+Na]-C5H14NO4P (-183)", 1.0);
//			_1_LPC_NA_RelativeAbundances_40.put("[M+Na]-Na-C5H13NO4P", 1.0);
//			_1_LPC_NA_RelativeAbundances_40.put("Choline", 1.0);
//			_1_LPC_NA_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 1.0);
//			_1_LPC_NA_RelativeAbundances_40.put("Sodiated 1,2-cyclic phosphate diester", 1.0);
//			_1_lpc.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_10", 10), _1_LPC_NA_RelativeAbundances_10));
//			_1_lpc.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_20", 20), _1_LPC_NA_RelativeAbundances_20));
//			_1_lpc.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_40", 40), _1_LPC_NA_RelativeAbundances_40));
//			
//			
//			LinkedHashMap<String, Double> _1_LPC_LI_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_LI_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_LI_RelativeAbundances_40 = new LinkedHashMap<String, Double>();		
//			_1_LPC_LI_RelativeAbundances_10.put("[M+Li]+", 1.0);
//			_1_LPC_LI_RelativeAbundances_10.put("[M+Li]-C3H9N", 1.0);
//			_1_LPC_LI_RelativeAbundances_20.put("[M+Li]+", 1.0);
//			_1_LPC_LI_RelativeAbundances_20.put("[M+Li]-C3H9N", 1.0);
//			_1_LPC_LI_RelativeAbundances_20.put("[M+Li]-Li-C5H13NO4P", 1.0);			
//			_1_LPC_LI_RelativeAbundances_40.put("[M+Li]+", 1.0);
//			_1_LPC_LI_RelativeAbundances_40.put("[M+Li]-C3H9N", 1.0);
//			_1_LPC_LI_RelativeAbundances_40.put("[M+Li]-C5H14[N+]O", 1.0);
//			_1_LPC_LI_RelativeAbundances_40.put("[M+Li]-C5H14NO4P (-183)", 1.0);
//			_1_LPC_LI_RelativeAbundances_40.put("[M+Li]-Li-C5H13NO4P", 1.0);
//			_1_LPC_LI_RelativeAbundances_40.put("Choline", 1.0);
//			_1_LPC_LI_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 1.0);
//			_1_lpc.put("[M+Li]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Li]+_10", 10), _1_LPC_LI_RelativeAbundances_10));
//			_1_lpc.put("[M+Li]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Li]+_20", 20), _1_LPC_LI_RelativeAbundances_20));
//			_1_lpc.put("[M+Li]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Li]+_40", 40), _1_LPC_LI_RelativeAbundances_40));
//
//			
//			LinkedHashMap<String, Double> _1_LPC_CL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_CL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPC_CL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();		
//			_1_LPC_CL_RelativeAbundances_10.put("[M+Cl]-", 1.0);
//			_1_LPC_CL_RelativeAbundances_10.put("[M+Cl]-CH3", 1.0);
//			_1_LPC_CL_RelativeAbundances_10.put("[M+Cl]-Cl-CH3", 1.0);
//			_1_LPC_CL_RelativeAbundances_20.put("[M+Cl]-", 1.0);
//			_1_LPC_CL_RelativeAbundances_20.put("[M+Cl]-CH3", 1.0);
//			_1_LPC_CL_RelativeAbundances_20.put("Sn FA-", 1.0);
//			_1_LPC_CL_RelativeAbundances_40.put("[M+Cl]-", 1.0);
//			_1_LPC_CL_RelativeAbundances_40.put("[M+Cl]-CH3", 1.0);
//			_1_LPC_CL_RelativeAbundances_40.put("Sn FA-", 1.0);
//			_1_LPC_CL_RelativeAbundances_40.put("Ion fragment C7H15NO4[O‐]P", 1.0);
//			_1_LPC_CL_RelativeAbundances_40.put("Ion fragment C6H13NO4[O‐]P", 1.0);
//			
//			
//			_1_lpc.put("[M+Cl]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_10", 10), _1_LPC_CL_RelativeAbundances_10));
//			_1_lpc.put("[M+Cl]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_20", 20), _1_LPC_CL_RelativeAbundances_20));
//			_1_lpc.put("[M+Cl]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_40", 40), _1_LPC_CL_RelativeAbundances_40));
//			classSpecificRelativeAbundances.put(ClassName._1_LYSOPHOSPHATIDYLCHOLINES, _1_lpc);
//			
//			
//			/*
//			 * Phosphatidylserines
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=624&LM_ID=LMGP03010027&TRACK_ID=278
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=760&LM_ID=LMGP03010024&TRACK_ID=281
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=762&LM_ID=LMGP03010024&TRACK_ID=282
//			 * 
//			 */
			LinkedHashMap<String,MSPeakRelativeAbundance> ps = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> PS_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PS_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PS_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PS_RelativeAbundances_10.put("[M+H]+", 1.0);
//			PS_RelativeAbundances_10.put("[M+H]-C3H8NO6P (-185)", 1.0);
//			PS_RelativeAbundances_20.put("[M+H]+", 1.0);
//			PS_RelativeAbundances_20.put("[M+H]-C3H8NO6P (-185)", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]+", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-C3H5NO2 (-87)", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-C3H8NO6P (-185)", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-C3H8NO6P (-185) – sn1", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-C3H8NO6P (-185) – sn2", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-sn1", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-sn2", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-sn1-H2O", 1.0);
//			PS_RelativeAbundances_40.put("[M+H]-sn2-H2O", 1.0);
//			
//			
//			ps.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+H]+_10", 10), PS_RelativeAbundances_10));
//			ps.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+H]+_20", 20), PS_RelativeAbundances_20));
//			ps.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+H]+_40", 40), PS_RelativeAbundances_40));
//			
//			
//
//			LinkedHashMap<String, Double> PS_LI_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PS_LI_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PS_LI_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PS_LI_RelativeAbundances_10.put("[M+Li]+", 1.0);
//			PS_LI_RelativeAbundances_10.put("[M+Li]-C3H5NO2 (-87)", 1.0);
//			PS_LI_RelativeAbundances_10.put("[M+Li]-C3H7NO6PLi", 1.0);
//			
//			PS_LI_RelativeAbundances_20.put("[M+Li]+", 1.0);
//			PS_LI_RelativeAbundances_20.put("[M+Li]-C3H5NO2 (-87)", 1.0);
//			PS_LI_RelativeAbundances_20.put("[M+Li]-C3H7NO6PLi", 1.0);
//			PS_LI_RelativeAbundances_20.put("Ion fragment C3H7NO6P[Li+]", 1.0);
//			
//			PS_LI_RelativeAbundances_40.put("[M+Li]+", 1.0);
//			PS_LI_RelativeAbundances_40.put("[M+Li]-C3H5NO2 (-87)", 1.0);
//			PS_LI_RelativeAbundances_40.put("[M+Li]-C3H7NO6PLi", 1.0);
//			PS_LI_RelativeAbundances_40.put("Ion fragment C3H7NO6P[Li+]", 1.0);
//			
//			
//			ps.put("[M+Li]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+Li]+_10", 10), PS_LI_RelativeAbundances_10));
//			ps.put("[M+Li]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+Li]+_20", 20), PS_LI_RelativeAbundances_20));
//			ps.put("[M+Li]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+Li]+_40", 40), PS_LI_RelativeAbundances_40));
//
//			
//			LinkedHashMap<String, Double> PS_NA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PS_NA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PS_NA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PS_NA_RelativeAbundances_10.put("[M+Na]+", 1.0);
//			PS_NA_RelativeAbundances_10.put("[M+Na]-C3H5NO2 (-87)", 1.0);
//			PS_NA_RelativeAbundances_10.put("[M+Na]-C3H7NO6PNa", 1.0);
//			
//			PS_NA_RelativeAbundances_20.put("[M+Na]+", 1.0);
//			PS_NA_RelativeAbundances_20.put("[M+Na]-C3H5NO2 (-87)", 1.0);
//			PS_NA_RelativeAbundances_20.put("[M+Na]-C3H7NO6PLi", 1.0);
//			
//			PS_NA_RelativeAbundances_40.put("[M+Na]+", 1.0);
//			PS_NA_RelativeAbundances_40.put("[M+Na]-C3H5NO2 (-87)", 1.0);
//			PS_NA_RelativeAbundances_40.put("[M+Na]-C3H7NO6PNa", 1.0);
//			PS_NA_RelativeAbundances_40.put("[M+Na]-H3PO4 (-98)", 1.0);
//			PS_NA_RelativeAbundances_40.put("Ion fragment C3H7NO6P[Na+]", 1.0);
//			PS_NA_RelativeAbundances_40.put("[M+Na]-C3H8NO6P (-185)", 1.0);
//			
//			
//			ps.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+Na]+_10", 10), PS_NA_RelativeAbundances_10));
//			ps.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+Na]+_20", 20), PS_NA_RelativeAbundances_20));
//			ps.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M+Na]+_40", 40), PS_NA_RelativeAbundances_40));
		
			
			
			LinkedHashMap<String, Double> PS_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PS_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PS_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
			PS_H_RelativeAbundances_10.put("[M-H]-", 100.0);
			PS_H_RelativeAbundances_20.put("[M-H]-", 90.0);
			PS_H_RelativeAbundances_20.put("[M-H]-C3H5NO2 (-87)", 15.0);
			PS_H_RelativeAbundances_40.put("[M-H]-", 60.0);
			PS_H_RelativeAbundances_40.put("[M-H]-C3H5NO2 (-87)", 90.0);
			PS_H_RelativeAbundances_40.put("[M-H-87]-sn1", 15.0);
			PS_H_RelativeAbundances_40.put("[M-H-87]-sn2", 15.0);
			PS_H_RelativeAbundances_40.put("[M-H-87]-sn1-H2O", 60.0);
			PS_H_RelativeAbundances_40.put("[M-H-87]-sn2-H2O", 60.0);
			PS_H_RelativeAbundances_40.put("Sn1 FA-H", 60.0);
			PS_H_RelativeAbundances_40.put("Sn2 FA-H", 60.0);
			PS_H_RelativeAbundances_40.put("Ion fragment C3H6O4[O‐]P", 15.0);
			ps.put("[M-H]-_10",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M-H]-_10", 10), PS_H_RelativeAbundances_10));
			ps.put("[M-H]-_20",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M-H]-_20", 20), PS_H_RelativeAbundances_20));
			ps.put("[M-H]-_40",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLSERINES, new FragmentationCondition("[M-H]-_40", 40), PS_H_RelativeAbundances_40));

			classSpecificRelativeAbundances.put(ClassName.PHOSPHATIDYLSERINES, ps);
			
			
			LinkedHashMap<String,MSPeakRelativeAbundance> cer = new LinkedHashMap<String,MSPeakRelativeAbundance>();
			LinkedHashMap<String, Double> CER_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> CER_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> CER_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
			CER_RelativeAbundances_10.put("[M+H]+", 90.0);
			CER_RelativeAbundances_10.put("[M-H2O+H]+", 15.0);
			CER_RelativeAbundances_20.put("[M+H]+", 90.0);
			CER_RelativeAbundances_20.put("[M-H2O+H]+", 15.0);
			CER_RelativeAbundances_20.put("[M-(H2O)2+H]-sn1", 60.0);
			CER_RelativeAbundances_40.put("[M-(H2O)2+H]+", 15.0);			
			CER_RelativeAbundances_40.put("[M+H]+", 15.0);
			CER_RelativeAbundances_40.put("[M-H2O+H]+", 60.0);
			CER_RelativeAbundances_40.put("[M-(H2O)2+H]+", 15.0);		
			CER_RelativeAbundances_40.put("[M-H2O+H]-CH2O", 15.0);
			CER_RelativeAbundances_40.put("[M-(H2O)2+H]-sn1", 90.0);
			CER_RelativeAbundances_40.put("[M-H2O+H]-sn1", 60.0);
			
			
			cer.put("[M+H]+_10",				
					new MSPeakRelativeAbundance(
							ClassName.CERAMIDES, new FragmentationCondition("[M+H]+_10", 10), CER_RelativeAbundances_10));
			cer.put("[M+H]+_20",				
					new MSPeakRelativeAbundance(
							ClassName.CERAMIDES, new FragmentationCondition("[M+H]+_20", 20), CER_RelativeAbundances_20));
			cer.put("[M+H]+_40",				
					new MSPeakRelativeAbundance(
							ClassName.CERAMIDES, new FragmentationCondition("[M+H]+_40", 40), CER_RelativeAbundances_40));
			

//			LinkedHashMap<String, Double> CER_LI_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CER_LI_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CER_LI_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			CER_LI_RelativeAbundances_10.put("[M+Li]+", 1.0);
//			CER_LI_RelativeAbundances_20.put("[M+Li]+", 1.0);
//			CER_LI_RelativeAbundances_40.put("[M+Li]+", 1.0);
//			CER_LI_RelativeAbundances_40.put("[M+Li]-O", 1.0);
//			CER_LI_RelativeAbundances_40.put("[M+Li]-48", 1.0);
//			cer.put("[M+Li]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.CERAMIDES, new FragmentationCondition("[M+Li]+_10", 10), CER_LI_RelativeAbundances_10));
//			cer.put("[M+Li]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.CERAMIDES, new FragmentationCondition("[M+Li]+_20", 20), CER_LI_RelativeAbundances_20));			
//			cer.put("[M+Li]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.CERAMIDES, new FragmentationCondition("[M+Li]+_40", 40), CER_LI_RelativeAbundances_40));		
//			
//			
//			LinkedHashMap<String, Double> CER_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CER_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CER_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			CER_H_RelativeAbundances_10.put("[M-H]-", 1.0);
//			CER_H_RelativeAbundances_20.put("[M-H]-", 1.0);
//			CER_H_RelativeAbundances_20.put("[M-H]-sn1-O", 1.0);
//			CER_H_RelativeAbundances_40.put("[M-H]-", 1.0);
//			CER_H_RelativeAbundances_40.put("[M-H-36]-", 1.0);
//			CER_H_RelativeAbundances_40.put("[M-H]-sn1", 1.0);
//			CER_H_RelativeAbundances_40.put("[M-H]-sn1-O", 1.0);
//			CER_H_RelativeAbundances_40.put("[M-H-282]-", 1.0);
//			CER_H_RelativeAbundances_40.put("[M-H-299.3]-", 1.0);
//			CER_H_RelativeAbundances_40.put("[M-H-327]-", 1.0);		
//			cer.put("[M-H]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.CERAMIDES, new FragmentationCondition("[M-H]-_10", 10), CER_H_RelativeAbundances_10));
//			cer.put("[M-H]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.CERAMIDES, new FragmentationCondition("[M-H]-_20", 20), CER_H_RelativeAbundances_20));
//			cer.put("[M-H]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.CERAMIDES, new FragmentationCondition("[M-H]-_40", 20), CER_H_RelativeAbundances_40));		
			
			classSpecificRelativeAbundances.put(ClassName.CERAMIDES, cer);
			
//			LinkedHashMap<String,MSPeakRelativeAbundance> sm = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> SM_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> SM_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> SM_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			SM_RelativeAbundances_10.put("[M+H]+", 1.0);
//			SM_RelativeAbundances_10.put("Phosphocholine", 1.0);
//			SM_RelativeAbundances_20.put("[M+H]+", 1.0);
//			SM_RelativeAbundances_20.put("Phosphocholine", 1.0);
//			SM_RelativeAbundances_40.put("[M+H]+", 1.0);
//			SM_RelativeAbundances_40.put("[M-H2O+H]+", 1.0);
//			SM_RelativeAbundances_40.put("Phosphocholine", 1.0);
//			sm.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+H]+_10", 10), SM_RelativeAbundances_10));
//			sm.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+H]+_20", 20), SM_RelativeAbundances_20));
//			sm.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+H]+_40", 40), SM_RelativeAbundances_40));
//			
//			
//			
//			LinkedHashMap<String, Double> SM_LI_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> SM_LI_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> SM_LI_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			SM_LI_RelativeAbundances_10.put("[M+Li]+", 1.0);
//			SM_LI_RelativeAbundances_10.put("[M+Li]-C3H9N", 1.0);
//			SM_LI_RelativeAbundances_10.put("[M+Li]-C5H14NO4P (-183)", 1.0);
//
//			SM_LI_RelativeAbundances_20.put("[M+Li]+", 1.0);
//			SM_LI_RelativeAbundances_20.put("[M+Li]-C3H9N", 1.0);
//			SM_LI_RelativeAbundances_20.put("[M+Li]-C5H14NO4P (-183)", 1.0);
//
//			SM_LI_RelativeAbundances_40.put("[M+Li]+", 1.0);
//			SM_LI_RelativeAbundances_40.put("[M+Li]-C3H9N", 1.0);
//			SM_LI_RelativeAbundances_40.put("[M+Li]-C5H14NO4P (-183)", 1.0);
//			SM_LI_RelativeAbundances_40.put("[M+Li]-C5H14NO4PLi", 1.0);
//			SM_LI_RelativeAbundances_40.put("[M+Li]-C5H14NO4PLi - H2O", 1.0);
//			
//			sm.put("[M+Li]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+Li]+_10", 10), SM_LI_RelativeAbundances_10));
//			sm.put("[M+Li]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+Li]+_20", 20), SM_LI_RelativeAbundances_20));			
//			sm.put("[M+Li]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+Li]+_40", 40), SM_LI_RelativeAbundances_40));			
//						
//			LinkedHashMap<String, Double> SM_NA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> SM_NA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> SM_NA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			SM_NA_RelativeAbundances_10.put("[M+Na]+", 90.0);
//			SM_NA_RelativeAbundances_10.put("[M+Na]-C3H9N", 15.0);
//			SM_NA_RelativeAbundances_10.put("[M+Na]-C5H14NO4P (-183)", 15.0);
//
//			SM_NA_RelativeAbundances_20.put("[M+Na]+", 90.0);
//			SM_NA_RelativeAbundances_20.put("[M+Na]-C3H9N", 15.0);
//			SM_NA_RelativeAbundances_20.put("[M+Na]-C5H14NO4P (-183)", 15.0);
//
//			SM_NA_RelativeAbundances_40.put("[M+Na]+", 60.0);
//			SM_NA_RelativeAbundances_40.put("[M+Na]-C3H9N", 90.0);
//			SM_NA_RelativeAbundances_40.put("[M+Na]-C5H14NO4P (-183)", 90.0);
//			SM_NA_RelativeAbundances_40.put("[M+Na]-C5H14NO4PNa", 15.0);
//			SM_NA_RelativeAbundances_40.put("[M+Na]-C5H14NO4PNa - H2O", 15.0);
//			
//			sm.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+Na]+_10", 10), SM_NA_RelativeAbundances_10));
//			sm.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+Na]+_20", 20), SM_NA_RelativeAbundances_20));			
//			sm.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.SPHINGOMYELINS, new FragmentationCondition("[M+Na]+_40", 40), SM_NA_RelativeAbundances_40));	
//			classSpecificRelativeAbundances.put(ClassName.SPHINGOMYELINS, sm);
			
			
			/*
			 * Cardiolipins
			 */
			
			LinkedHashMap<String,MSPeakRelativeAbundance> cl = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> CL_H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CL_H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CL_H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			CL_H_RelativeAbundances_10.put("[M-H]-", 1.0);
//			CL_H_RelativeAbundances_20.put("[M-H]-", 1.0);
//			CL_H_RelativeAbundances_40.put("[M-H]-", 1.0);
//			cl.put("[M-H]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.CARDIOLIPINS, new FragmentationCondition("[M-H]-_10", 10), CL_H_RelativeAbundances_10));
//			cl.put("[M-H]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.CARDIOLIPINS, new FragmentationCondition("[M-H]-_20", 20), CL_H_RelativeAbundances_20));			
//			cl.put("[M-H]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.CARDIOLIPINS, new FragmentationCondition("[M-H]-_40", 40), CL_H_RelativeAbundances_40));	

			LinkedHashMap<String, Double> CL_2H_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> CL_2H_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> CL_2H_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
			CL_2H_RelativeAbundances_10.put("[M-2H](2H)-", 90.0);
			CL_2H_RelativeAbundances_10.put("Sn1 FA", 15.0);
			CL_2H_RelativeAbundances_10.put("Sn2 FA", 15.0);
			CL_2H_RelativeAbundances_10.put("Sn3 FA", 15.0);
			CL_2H_RelativeAbundances_10.put("Sn4 FA", 15.0);		
			CL_2H_RelativeAbundances_20.put("[M-2H](2H)-", 90.0);
			CL_2H_RelativeAbundances_20.put("[M-2H](2-)-sn1 (monolyso)", 15.0);
			CL_2H_RelativeAbundances_20.put("[M-2H](2-)-sn3 (monolyso)", 15.0);
			CL_2H_RelativeAbundances_20.put("[sn1+C3H6PO4-H](-) (137.00)", 15.0);
			CL_2H_RelativeAbundances_20.put("[sn2+C3H6PO4-H](-) (137.00)", 15.0);
			CL_2H_RelativeAbundances_20.put("[sn3+C3H6PO4-H](-) (137.00)", 15.0);
			CL_2H_RelativeAbundances_20.put("[sn4+C3H6PO4-H](-) (137.00)", 15.0);
			CL_2H_RelativeAbundances_20.put("Sn1 FA", 15.0);
			CL_2H_RelativeAbundances_20.put("Sn2 FA", 15.0);
			CL_2H_RelativeAbundances_20.put("Sn3 FA", 15.0);
			CL_2H_RelativeAbundances_20.put("Sn4 FA", 15.0);
			CL_2H_RelativeAbundances_40.put("[M-2H](2H)-", 15.0);
			CL_2H_RelativeAbundances_40.put("[M-2H](2-)-sn1 (monolyso)", 15.0);
			CL_2H_RelativeAbundances_40.put("[M-2H](2-)-sn3 (monolyso)", 15.0);
			CL_2H_RelativeAbundances_40.put("[M-2H](2-) -(sn1+sn2 dilyso)", 15.0);
			CL_2H_RelativeAbundances_40.put("[M-2H](2-) -(sn3+sn4 dilyso)", 15.0);
			CL_2H_RelativeAbundances_40.put("[sn1+C3H6PO4-H](-) (137.00)", 60.0);
			CL_2H_RelativeAbundances_40.put("[sn2+C3H6PO4-H](-) (137.00)", 60.0);
			CL_2H_RelativeAbundances_40.put("[sn3+C3H6PO4-H](-) (137.00)", 60.0);
			CL_2H_RelativeAbundances_40.put("[sn4+C3H6PO4-H](-) (137.00)", 60.0);
			CL_2H_RelativeAbundances_40.put("Sn1 FA", 90.0);
			CL_2H_RelativeAbundances_40.put("Sn2 FA", 90.0);
			CL_2H_RelativeAbundances_40.put("Sn3 FA", 90.0);
			CL_2H_RelativeAbundances_40.put("Sn4 FA", 90.0);
			
			cl.put("[M-2H](2H)-_10",				
					new MSPeakRelativeAbundance(
							ClassName.CARDIOLIPINS, new FragmentationCondition("[M-2H](2H)-_10", 10), CL_2H_RelativeAbundances_10));
			cl.put("[M-2H](2H)-_20",				
					new MSPeakRelativeAbundance(
							ClassName.CARDIOLIPINS, new FragmentationCondition("[M-2H](2H)-_20", 20), CL_2H_RelativeAbundances_20));			
			cl.put("[M-2H](2H)-_40",				
					new MSPeakRelativeAbundance(
							ClassName.CARDIOLIPINS, new FragmentationCondition("[M-2H](2H)-_40", 40), CL_2H_RelativeAbundances_40));					
			
			
			classSpecificRelativeAbundances.put(ClassName.CARDIOLIPINS, cl);
		
		
		
			/*
			 * Acylcarnitines
			 */
//			LinkedHashMap<String,MSPeakRelativeAbundance> car = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> CAR_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CAR_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CAR_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			CAR_RelativeAbundances_10.put("[M+]", 1.0);
//			CAR_RelativeAbundances_10.put("[M+]-C3H9N", 1.0);
//			CAR_RelativeAbundances_20.put("[M+]", 1.0);
//			CAR_RelativeAbundances_20.put("[M+]-C3H9N", 1.0);
//			CAR_RelativeAbundances_40.put("[M+]", 1.0);
//			CAR_RelativeAbundances_40.put("[M+]-C3H9N", 1.0);
//			CAR_RelativeAbundances_40.put("Sn1 FA", 1.0);
//			CAR_RelativeAbundances_40.put("Sn1-H", 1.0);
//			CAR_RelativeAbundances_40.put("[M+]-sn1-H2O", 1.0);
//			CAR_RelativeAbundances_40.put("Ion fragment C4H5O2+", 1.0);
//
//			car.put("[M+]_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.ACYL_CARNITINES, new FragmentationCondition("[M+]_10", 10), CAR_RelativeAbundances_10));
//			car.put("[M+]_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.ACYL_CARNITINES, new FragmentationCondition("[M+]_20", 20), CAR_RelativeAbundances_20));
//			car.put("[M+]_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.ACYL_CARNITINES, new FragmentationCondition("[M+]_40", 40), CAR_RelativeAbundances_40));		
//			
//			
//			LinkedHashMap<String, Double> CAR_LI_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CAR_LI_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> CAR_LI_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			CAR_LI_RelativeAbundances_10.put("[M+Li]+", 1.0);
//			CAR_LI_RelativeAbundances_10.put("[M+Li]-C3H9N", 1.0);						
//			CAR_LI_RelativeAbundances_20.put("[M+Li]+", 1.0);
//			CAR_LI_RelativeAbundances_20.put("[M+Li]-C3H9N", 1.0);
//			CAR_LI_RelativeAbundances_20.put("Sn1 FA", 1.0);		
//			CAR_LI_RelativeAbundances_40.put("[M+Li]+", 1.0);
//			CAR_LI_RelativeAbundances_40.put("[M+Li]-C3H9N", 1.0);
//			CAR_LI_RelativeAbundances_40.put("[M+Li]-H2O-C3H9N (-18-59)", 1.0);
//			CAR_LI_RelativeAbundances_40.put("Sn1 FA", 1.0);
//			CAR_LI_RelativeAbundances_40.put("[M+Li]-Sn1-H2O", 1.0);
//			CAR_LI_RelativeAbundances_40.put("Ion fragment C4H4LiO2+", 1.0);
//			car.put("[M+Li]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.ACYL_CARNITINES, new FragmentationCondition("[M+Li]+_10", 10), CAR_LI_RelativeAbundances_10));	
//			car.put("[M+Li]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.ACYL_CARNITINES, new FragmentationCondition("[M+Li]+_20", 20), CAR_LI_RelativeAbundances_20));
//			car.put("[M+Li]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.ACYL_CARNITINES, new FragmentationCondition("[M+Li]+_40", 40), CAR_LI_RelativeAbundances_40));	
//			
//			classSpecificRelativeAbundances.put(ClassName.ACYL_CARNITINES, car);
			
			
			
			/*
			 * Phosphatidylglycerols
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=747&LM_ID=LMGP04010002&TRACK_ID=296
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=769&LM_ID=LMGP04010036&TRACK_ID=297
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=705&LM_ID=LMGP04010007&TRACK_ID=151
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=783&LM_ID=LMGP04010006&TRACK_ID=134
			 */
			
			LinkedHashMap<String,MSPeakRelativeAbundance> pgl = new LinkedHashMap<String,MSPeakRelativeAbundance>();
			LinkedHashMap<String, Double> PGL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PGL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
			LinkedHashMap<String, Double> PGL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();			
			PGL_RelativeAbundances_10.put("[M-H]-", 100.0);
			PGL_RelativeAbundances_20.put("[M-H]-", 90.0);			
			PGL_RelativeAbundances_20.put("[M-H]-sn1", 15.0);
			PGL_RelativeAbundances_20.put("[M-H]-sn1-H2O", 15.0);
			PGL_RelativeAbundances_20.put("[M-H]-sn2", 15.0);
			PGL_RelativeAbundances_20.put("[M-H]-sn2-H2O", 15.0);
			PGL_RelativeAbundances_20.put("[M-H]-sn1-C3H8O3", 15.0);
			PGL_RelativeAbundances_20.put("[M-H]-sn2-C3H8O3", 15.0);
			PGL_RelativeAbundances_20.put("Sn1 FA-H", 90.0);
			PGL_RelativeAbundances_20.put("Sn2 FA-H", 90.0);
			PGL_RelativeAbundances_20.put("Ion fragment C3H6O4[O‐]P", 15.0);				
			PGL_RelativeAbundances_40.put("[M-H]-", 90.0);
			PGL_RelativeAbundances_40.put("[M-H]-sn1", 15.0);
			PGL_RelativeAbundances_40.put("[M-H]-sn1-H2O", 15.0);
			PGL_RelativeAbundances_40.put("[M-H]-sn2", 15.0);
			PGL_RelativeAbundances_40.put("[M-H]-sn2-H2O", 15.0);
			PGL_RelativeAbundances_40.put("[M-H]-sn1-C3H8O3", 15.0);
			PGL_RelativeAbundances_40.put("[M-H]-sn2-C3H8O3", 15.0);
			PGL_RelativeAbundances_40.put("Sn1 FA-H", 90.0);
			PGL_RelativeAbundances_40.put("Sn2 FA-H", 90.0);
			PGL_RelativeAbundances_40.put("Ion fragment C3H6O4[O‐]P", 15.0);	
			PGL_RelativeAbundances_40.put("Ion fragment C3H8O6P-", 15.0);	
			pgl.put("[M-H]-_10",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLGLYCEROLS, new FragmentationCondition("[M-H]-_10", 10), PGL_RelativeAbundances_10));	
			pgl.put("[M-H]-_20",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLGLYCEROLS, new FragmentationCondition("[M-H]-_20", 20), PGL_RelativeAbundances_20));				
			pgl.put("[M-H]-_40",				
					new MSPeakRelativeAbundance(
							ClassName.PHOSPHATIDYLGLYCEROLS, new FragmentationCondition("[M-H]-_40", 40), PGL_RelativeAbundances_40));
			
			classSpecificRelativeAbundances.put(ClassName.PHOSPHATIDYLGLYCEROLS, pgl);
			
			
			
			/*
			 * Lysophosphatidylglycerols
			 * 
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=509&LM_ID=LMGP04050006&TRACK_ID=484
			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=441&LM_ID=LMGP04050007&TRACK_ID=481
			 */
//			LinkedHashMap<String,MSPeakRelativeAbundance> _1_lpgl = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> _1_LPGL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPGL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_LPGL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();	
//			_1_LPGL_RelativeAbundances_10.put("[M-H]-", 1.0);
//			_1_LPGL_RelativeAbundances_20.put("[M-H]-", 1.0);
//			_1_LPGL_RelativeAbundances_40.put("[M-H]-", 1.0);
//			_1_LPGL_RelativeAbundances_40.put("[M-H]-C3H7O2", 1.0);
//			_1_LPGL_RelativeAbundances_40.put("[M-H]-C3H8O3", 1.0);
//			_1_LPGL_RelativeAbundances_40.put("Sn1 FA-H", 1.0);
//			_1_LPGL_RelativeAbundances_40.put("Ion fragment C6H14O8P-", 1.0);
//			_1_LPGL_RelativeAbundances_40.put("Ion fragment C6H12O7P-", 1.0);
//			_1_LPGL_RelativeAbundances_40.put("Ion fragment C3H6O5P-", 1.0);
//			
//			_1_lpgl.put("[M-H]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLGLYCEROLS, new FragmentationCondition("[M-H]-_10", 10), _1_LPGL_RelativeAbundances_10));	
//			_1_lpgl.put("[M-H]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLGLYCEROLS, new FragmentationCondition("[M-H]-_20", 20), _1_LPGL_RelativeAbundances_20));
//			_1_lpgl.put("[M-H]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_LYSOPHOSPHATIDYLGLYCEROLS, new FragmentationCondition("[M-H]-_40", 40), _1_LPGL_RelativeAbundances_40));	
//
//			classSpecificRelativeAbundances.put(ClassName._1_LYSOPHOSPHATIDYLGLYCEROLS, _1_lpgl);
//			
//			
//			/*
//			 * Plasmanyl-phosphatidylcholines
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=770&LM_ID=LMGP01020053&TRACK_ID=230
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=768&LM_ID=LMGP01020056&TRACK_ID=228
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=804&LM_ID=LMGP01020053&TRACK_ID=229
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=802&LM_ID=LMGP01020056&TRACK_ID=227
//			 * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4331342/
//			 */
//			LinkedHashMap<String,MSPeakRelativeAbundance> plasmapc = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> PLASMAPC_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMAPC_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMAPC_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PLASMAPC_RelativeAbundances_10.put("[M+H]+", 1.0);
//			PLASMAPC_RelativeAbundances_20.put("[M+H]+", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("[M+H]+", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("[M+H]-C3H9N", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("[M+H]-C5H14NO4P (-183)", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("[M+H]-sn1", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("[M+H]-sn2", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("[M+H]-sn1-H2O", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("[M+H]-sn2-H2O", 1.0);			
//			PLASMAPC_RelativeAbundances_40.put("Phosphocholine", 1.0);
//			PLASMAPC_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 1.0);
//			
//			
//			plasmapc.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMANYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_10", 10), PLASMAPC_RelativeAbundances_10));	
//			plasmapc.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMANYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_20", 20), PLASMAPC_RelativeAbundances_20));	
//			plasmapc.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMANYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_40", 40), PLASMAPC_RelativeAbundances_40));
//			
//			
//			LinkedHashMap<String, Double> PLASMAPC_CL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMAPC_CL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMAPC_CL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PLASMAPC_CL_RelativeAbundances_10.put("[M+Cl]-", 1.0);
//			PLASMAPC_CL_RelativeAbundances_20.put("[M+Cl]-", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("[M+Cl]-", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("[M+Cl]-C3H9N", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-CH3", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("[M+H]-C5H13N-Cl (-122)", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-sn2-CH3", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-sn2-H2O-CH3", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("Sn2 FA-H", 1.0);
//			PLASMAPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-sn2-H2O-C5H13N", 1.0);
//			plasmapc.put("[M+Cl]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMANYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_10", 10), PLASMAPC_CL_RelativeAbundances_10));	
//			plasmapc.put("[M+Cl]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMANYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_20", 20), PLASMAPC_CL_RelativeAbundances_20));	
//			plasmapc.put("[M+Cl]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMANYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_40", 40), PLASMAPC_CL_RelativeAbundances_40));				
//			
//			
//			classSpecificRelativeAbundances.put(ClassName.PLASMANYLPHOSPHATIDYLCHOLINES, plasmapc);
//			
//			
//			/*
//			 * Plasemenyl-Phosphatidylcholines
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=794&LM_ID=LMGP01030012&TRACK_ID=232
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=818&LM_ID=LMGP01030014&TRACK_ID=519
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=828&LM_ID=LMGP01030012&TRACK_ID=231
//			 * http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=852&LM_ID=LMGP01030014&TRACK_ID=475
//			 * Fong-Fu Hsu et al. (2014); J. Am. Soc. Mass Spectrom. (2014) 25:1412Y1420; DOI: 10.1007/s13361-014-0908-x
//			 */
//			
//			LinkedHashMap<String,MSPeakRelativeAbundance> plasmneylpc = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> PLASMENYLPC_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMENYLPC_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMENYLPC_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PLASMENYLPC_RelativeAbundances_10.put("[M+H]+", 1.0);
//			PLASMENYLPC_RelativeAbundances_20.put("[M+H]+", 1.0);
//			PLASMENYLPC_RelativeAbundances_20.put("Phosphocholine", 1.0);
//			PLASMENYLPC_RelativeAbundances_20.put("N,N,N-Trimethylethenaminium", 1.0);			
//			PLASMENYLPC_RelativeAbundances_40.put("[M+H]+", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("[M+H]-C3H9N (-59)", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("[M+H]-C5H14NO4P (-183)", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("[M+H]-sn1", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("[M+H]-sn2", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("[M+H]-sn2-H2O", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("[M+H]-sn1-C3H9N (-59)", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("Phosphocholine", 1.0);
//			PLASMENYLPC_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 1.0);
//
//			plasmneylpc.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMENYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_10", 10), PLASMENYLPC_RelativeAbundances_10));	
//			plasmneylpc.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMENYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_20", 20), PLASMENYLPC_RelativeAbundances_20));	
//			plasmneylpc.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMENYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_40", 40), PLASMENYLPC_RelativeAbundances_40));	
//			
//			
//			LinkedHashMap<String, Double> PLASMENYLPC_CL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMENYLPC_CL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> PLASMENYLPC_CL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			PLASMENYLPC_CL_RelativeAbundances_10.put("[M+Cl]-", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_20.put("[M+Cl]-", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("[M+Cl]-", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-CH3", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("[M+Cl]-C5H13N-Cl (-122)", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-sn2-CH3", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-sn2-H2O-CH3", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("Sn2 FA-H", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-sn1-CH3", 1.0);
//			PLASMENYLPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-sn1-H2O-CH3", 1.0);
//			
//			
//			plasmneylpc.put("[M+Cl]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMENYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_10", 10), PLASMENYLPC_CL_RelativeAbundances_10));	
//			plasmneylpc.put("[M+Cl]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMENYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_20", 20), PLASMENYLPC_CL_RelativeAbundances_20));				
//			plasmneylpc.put("[M+Cl]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName.PLASMENYLPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_40", 40), PLASMENYLPC_CL_RelativeAbundances_40));				
//			
//			classSpecificRelativeAbundances.put(ClassName.PLASMENYLPHOSPHATIDYLCHOLINES, plasmneylpc);
//			
//			
//			/*
//			 * 1-alkyl-lysophosphatidylcholines
//			 * 10 eV
//			 * 		[M+H]+   http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=508&LM_ID=LMGP01060039&TRACK_ID=509
//			 * 20 eV
//			 * 		[M+H]+ http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=510&LM_ID=LMGP01060014&TRACK_ID=507
//			 * 40 eV
//			 * 		[M+H]+  http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=482&LM_ID=LMGP01060010&TRACK_ID=237
//			 * 		[M+H]+  http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=482&LM_ID=LMGP01060010&TRACK_ID=505
//			 * 		[M+Cl]-  http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=544&LM_ID=LMGP01060014&TRACK_ID=506
//			 * 		[M+Cl]-  http://lipidmaps.org/data/standards/fetch_gif_mult.php?MASS=542&LM_ID=LMGP01060039&TRACK_ID=508
//			 */
//			
//			LinkedHashMap<String,MSPeakRelativeAbundance> _1alkyl_lysopc = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			_1_ALKYL_LYSOPC_RelativeAbundances_10.put("[M+H]+", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_10.put("[M-H+H2O]+", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_10.put("Ion fragment C8H19[N+]O5P", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_10.put("Phosphocholine", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_10.put("N,N,N-Trimethylethenaminium", 1.0);					
//			_1_ALKYL_LYSOPC_RelativeAbundances_20.put("[M+H]+", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_20.put("[M-H+H2O]+", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_20.put("Ion fragment C8H19[N+]O5P", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_20.put("Phosphocholine", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_20.put("N,N,N-Trimethylethenaminium", 1.0);	
//			_1_ALKYL_LYSOPC_RelativeAbundances_40.put("[M+H]+", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_40.put("[M-H+H2O]+", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_40.put("Ion fragment C8H19[N+]O5P", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_40.put("1,2-cyclic phosphate diester", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_40.put("Phosphocholine", 1.0);
//			_1_ALKYL_LYSOPC_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 1.0);			
//			_1alkyl_lysopc.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_10", 10), _1_ALKYL_LYSOPC_RelativeAbundances_10));	
//			_1alkyl_lysopc.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_20", 20), _1_ALKYL_LYSOPC_RelativeAbundances_20));	
//			_1alkyl_lysopc.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_40", 40), _1_ALKYL_LYSOPC_RelativeAbundances_40));	
//			
//			
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_CL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_CL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_CL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();			
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_10.put("[M+Cl]-", 1.0);			
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_20.put("[M+Cl]-", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_20.put("[M+Cl]-Cl-CH3", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_20.put("[M+Cl]-sn1-H2O-CH3", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_20.put("Ion fragment C4H11NO4P-", 1.0);			
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-CH3", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-C5H14N", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-C5H14N-H2O", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-sn1-H2O-CH3", 1.0);
//			_1_ALKYL_LYSOPC_CL_RelativeAbundances_40.put("Ion fragment C4H11NO4P-", 1.0);
//			_1alkyl_lysopc.put("[M+Cl]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_10", 10), _1_ALKYL_LYSOPC_CL_RelativeAbundances_10));	
//			_1alkyl_lysopc.put("[M+Cl]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_20", 20), _1_ALKYL_LYSOPC_CL_RelativeAbundances_20));				
//			_1alkyl_lysopc.put("[M+Cl]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_40", 40), _1_ALKYL_LYSOPC_CL_RelativeAbundances_40));				
//			
//			
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_NA_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_NA_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKYL_LYSOPC_NA_RelativeAbundances_40 = new LinkedHashMap<String, Double>();			
//			_1_ALKYL_LYSOPC_NA_RelativeAbundances_10.put("[M+Na]+", 1.0);
//			_1_ALKYL_LYSOPC_NA_RelativeAbundances_20.put("[M+Na]+", 1.0);
//			_1_ALKYL_LYSOPC_NA_RelativeAbundances_40.put("[M+Na]+", 1.0);
//			_1_ALKYL_LYSOPC_NA_RelativeAbundances_40.put("[M+Na]-H2O", 1.0);
//			_1_ALKYL_LYSOPC_NA_RelativeAbundances_40.put("[M+Na]-C3H9N", 1.0);
//			_1alkyl_lysopc.put("[M+Na]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_10", 10), _1_ALKYL_LYSOPC_NA_RelativeAbundances_10));	
//			_1alkyl_lysopc.put("[M+Na]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_20", 20), _1_ALKYL_LYSOPC_NA_RelativeAbundances_20));	
//			_1alkyl_lysopc.put("[M+Na]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Na]+_40", 40), _1_ALKYL_LYSOPC_NA_RelativeAbundances_40));	
//				
//			
//			classSpecificRelativeAbundances.put(ClassName._1_ALKYL_LYSOPHOSPHATIDYLCHOLINES, _1alkyl_lysopc);
//			
//			
//			LinkedHashMap<String,MSPeakRelativeAbundance> _1alkenyl_lysopc = new LinkedHashMap<String,MSPeakRelativeAbundance>();
//			LinkedHashMap<String, Double> _1_ALKENYL_LYSOPC_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKENYL_LYSOPC_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKENYL_LYSOPC_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("[M+H]+", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("[M-H2O+H]+", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("[M+H]-sn1-H2O", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("[M+H]-sn1-(H2O)2", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("Phosphocholine", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("[M+H]-sn1-H2O-C3H9N", 1.0);		
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("Choline", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_10.put("N,N,N-Trimethylethenaminium", 1.0);			
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("[M+H]+", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("[M-H2O+H]+", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("[M+H]-sn1-H2O", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("[M+H]-sn1-(H2O)2", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("Phosphocholine", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("[M+H]-sn1-H2O-C3H9N", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("Choline", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_20.put("N,N,N-Trimethylethenaminium", 1.0);	
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("[M+H]+", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("[M-H2O+H]+", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("[M+H]-sn1-H2O", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("[M+H]-sn1-(H2O)2", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("Phosphocholine", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("[M+H]-sn1-H2O-C3H9N", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("Choline", 1.0);
//			_1_ALKENYL_LYSOPC_RelativeAbundances_40.put("N,N,N-Trimethylethenaminium", 1.0);	
//			_1alkenyl_lysopc.put("[M+H]+_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKENYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_10", 10), _1_ALKENYL_LYSOPC_RelativeAbundances_10));	
//			_1alkenyl_lysopc.put("[M+H]+_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKENYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_20", 20), _1_ALKENYL_LYSOPC_RelativeAbundances_20));	
//			_1alkenyl_lysopc.put("[M+H]+_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKENYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+H]+_40", 40), _1_ALKENYL_LYSOPC_RelativeAbundances_40));	
//			
//			
//			LinkedHashMap<String, Double> _1_ALKENYL_LYSOPC_CL_RelativeAbundances_10 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKENYL_LYSOPC_CL_RelativeAbundances_20 = new LinkedHashMap<String, Double>();
//			LinkedHashMap<String, Double> _1_ALKENYL_LYSOPC_CL_RelativeAbundances_40 = new LinkedHashMap<String, Double>();
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_10.put("[M+Cl]-", 1.0);	
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_20.put("[M+Cl]-", 1.0);
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_20.put("[M+Cl]-Cl-CH3", 1.0);
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_20.put("[M+Cl]-Cl-C5H14N-H2O (-139)", 1.0);	
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_20.put("[M+Cl]-sn1-H2O-CH3", 1.0);
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_20.put("Ion fragment C8H19[N+]O5P", 1.0);
//
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-", 1.0);
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-CH3", 1.0);
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-Cl-C5H14N-H2O (-139)", 1.0);
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_40.put("[M+Cl]-sn1-H2O-CH3", 1.0);
//			_1_ALKENYL_LYSOPC_CL_RelativeAbundances_40.put("Ion fragment C8H19[N+]O5P", 1.0);
//			
//			
//			_1alkenyl_lysopc.put("[M+Cl]-_10",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKENYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_10", 10), _1_ALKENYL_LYSOPC_CL_RelativeAbundances_10));	
//			_1alkenyl_lysopc.put("[M+Cl]-_20",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKENYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_20", 20), _1_ALKENYL_LYSOPC_CL_RelativeAbundances_20));				
//			_1alkenyl_lysopc.put("[M+Cl]-_40",				
//					new MSPeakRelativeAbundance(
//							ClassName._1_ALKENYL_LYSOPHOSPHATIDYLCHOLINES, new FragmentationCondition("[M+Cl]-_40", 40), _1_ALKENYL_LYSOPC_CL_RelativeAbundances_40));	
//			
//			
//			
//			
//			classSpecificRelativeAbundances.put(ClassName._1_ALKENYL_LYSOPHOSPHATIDYLCHOLINES, _1alkenyl_lysopc);
		}
}
