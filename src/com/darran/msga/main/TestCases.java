package com.darran.msga.main;

import static java.lang.System.out;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;

import com.darran.msga.population.Alignment;
import com.darran.msga.population.Individual;
import com.darran.msga.scorer.DNAScorer;

/**
 * Sample problems to compare GAMSA to ClustalW
 * @author L.danran
 * Created on Apr 22, 2014
 */
public class TestCases
{
	//Scorers vary by problem, so we need to reset the scorer used through the InfoCenter.
	private static InfoCenter i_center = InfoCenter.getCenter();
	
	/**
	 * Select which test to run
	 */
	public static void main(String[] args) throws IOException
	{
		MultiSeqAligner aligner = new MultiSeqAligner();
		String[] inputSequences = null;
		Alignment solution = null;
		
//		out.println("Please select a test to run");
//		out.println("  0) Amino Acid Test 1");
//		out.println("  1) Amino Acid Test 2");		
//		out.println("  2) DNA Mickey Mouse Sequences");
//		out.println("  3) DNA MYH16 Sequences");
//		out.println("  4) DNA Beta Globin Sequences");
//		out.println("  5) DNA HIV Sequences");
//		out.println("  6) DNA BRCA1 Sequences\n");
//		int c = System.in.read();
//		switch (c)
//		{
//			case '0':
//			{
//				i_center.setScorer(new Blosum62Scorer());				
//				inputSequences = SequenceUtils.getAAHemoGlobinSequences();
//				solution = SequenceUtils.getAAHemoGlobinClustalWSolution();
//				break;
//			}
//			case '1':
//			{
//				i_center.setScorer(new Blosum62Scorer());				
//				inputSequences = SequenceUtils.getAAGrowthHormoneSequences();
//				solution = SequenceUtils.getAAGrowthHormoneClustalWSolution();
//				break;
//			}			
//			case '2':
//			{
//				i_center.setScorer(new DNAScorer());				
//				inputSequences = SequenceUtils.getDNAMickeyMouseSequences();
//				solution = SequenceUtils.getDNAMickeyMouseClustalWSolution();
//				break;
//			}			
//			case '3':
//			{
//				i_center.setScorer(new DNAScorer());				
//				inputSequences = SequenceUtils.getDNAMYH16Sequences();
//				solution = SequenceUtils.getDNAMYH16ClustalWSolution();
//				break;
//			}
//			case '4':
//			{
				i_center.setScorer(new DNAScorer());				
				inputSequences = SequenceUtils.getDNABetaGlobinSequences();
				solution = SequenceUtils.getDNABetaGlobinClustalWSolution();
//				break;
//			}
//			case '5':
//			{
//				i_center.setScorer(new DNAScorer());				
//				inputSequences = SequenceUtils.getDNAHIVSequences();
//				solution = SequenceUtils.getDNAHIVClustalWSolution();
//				break;
//			}
//			case '6':
//			{
//				i_center.setScorer(new DNAScorer());				
//				inputSequences = SequenceUtils.getDNABRCA1Sequences();
//				solution = SequenceUtils.getDNABRCA1ClustalWSolution();
//				break;
//			}			
//		}		
		
		long startTime = System.currentTimeMillis();		
		Individual best = aligner.findSolution(inputSequences);
		long endTime = System.currentTimeMillis();		
					
		Date date = new Date(endTime - startTime);
		SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss:SS");
		sdf.setTimeZone(new java.util.SimpleTimeZone(0, "UTC"));
		
		out.println("******");
		out.println("MSGA Best " + best);				 
//		out.println();
//		out.println("ClustalW Solution " + solution); 
		out.println();
		out.println("Statistics:");
		out.println("Time Elapsed (HH:mm:ss:SS): " + sdf.format(date));
//		Properties properties = MultiSeqAligner.getDefaultProperties();
//		out.println("Population Size: " + properties.getProperty("populationSize"));
//		out.println("Number of Unchanged Rounds: " + properties.getProperty("unchangedRoundsNeeded"));
//		out.println("Number of Total Rounds: " + aligner.getNumberOfRoundsExecuted());
//		out.println("Goodbye, and thank you for flying GAMSA.");
		int seqNum = inputSequences.length;
		int seqAvgLength = 0;
		for(String seq : inputSequences){
			seqAvgLength += seq.length();
		}
		seqAvgLength /= seqNum;
		out.println("Number of sequences: " + seqNum);
		out.println("Sequence average length: " + seqAvgLength);
		out.println("******");
	}

}