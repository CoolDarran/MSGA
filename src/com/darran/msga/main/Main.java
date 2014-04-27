package com.darran.msga.main;

import static java.lang.System.out;

import java.text.SimpleDateFormat;
import java.util.Date;

import mpi.MPI;

import com.darran.msga.constant.Constant;
import com.darran.msga.population.Alignment;
import com.darran.msga.population.Individual;
import com.darran.msga.population.Population;
import com.darran.msga.scorer.DNAScorer;

public class Main {
	private static InfoCenter i_center = InfoCenter.getCenter();
	private static MultiSeqAligner aligner = new MultiSeqAligner();
	
	public static void main(String[] args) throws Exception{		
		String[] inputSequences = null;
		Alignment solution = null;		
		
		i_center.setScorer(new DNAScorer());				
		inputSequences = SequenceUtils.getDNABetaGlobinSequences();
		solution = SequenceUtils.getDNABetaGlobinClustalWSolution();
		
		doMPIAligner(args,i_center,inputSequences,aligner,solution);
	}
	
	private static void doMPIAligner(String[] args, InfoCenter i_center2,
			String[] inputSequences, MultiSeqAligner aligner, Alignment solution) {
		
		MPI.Init(args);
		int me = MPI.COMM_WORLD.Rank();
		int size = MPI.COMM_WORLD.Size();
		
		// MPI的对象传递需通过array，传递的对象需要进行Serializable
		Alignment r[] = new Alignment[20];
		switch(me){
		// HRM
		case 0:{
			long startTime = System.currentTimeMillis();
			
			int[] bGAs = {3,4};
			int[] lRMs = {1,2};
			HRM hrm = new HRM(bGAs, lRMs);
			Individual best = hrm.mainDo(inputSequences);
			
			long endTime = System.currentTimeMillis();
			Date date = new Date(endTime - startTime);
			SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss:SS");
			sdf.setTimeZone(new java.util.SimpleTimeZone(0, "UTC"));
			
			out.println("******");
			out.println("MSGA Best " + best);							 
			out.println();
//			out.println("ClustalW Solution " + solution); 			
//			out.println();
			out.println("Statistics:");
			out.println("Time Elapsed (HH:mm:ss:SS): " + sdf.format(date));
//			Properties properties = MultiSeqAligner.getDefaultProperties();
//			out.println("Population Size: " + properties.getProperty("populationSize"));
//			out.println("Number of Unchanged Rounds: " + properties.getProperty("unchangedRoundsNeeded"));
//			out.println("Number of Total Rounds: " + aligner.getNumberOfRoundsExecuted());
			int seqNum = inputSequences.length;
			int seqAvgLength = 0;
			for(String seq : inputSequences){
				seqAvgLength += seq.length();
			}
			seqAvgLength /= seqNum;
			out.println("Number of sequences: " + seqNum);
			out.println("Sequence average length: " + seqAvgLength);
			out.println("******");
			break;			
		}
		// LRM1
		case 1:{					
			int[] bGAs = {5,6};
			int hRM = 0, rlrm = 2,llrm = 2;
			LRM lrm = new LRM(bGAs, hRM, rlrm, llrm);
			lrm.mainDo(inputSequences);
			break;
		}
		// LRM2
		case 2:{
			int[] bGAs = {7,8};
			int hRM = 0, rlrm = 1,llrm = 1;
			LRM lrm = new LRM(bGAs, hRM, rlrm, llrm);
			lrm.mainDo(inputSequences);
			break;
		}
		// BGA1 H
		case 3:{
			int hRM = 0;
			BGA4H(hRM);
			break;
		}
		// BGA2 H
		case 4: {
			int hRM = 0;
			BGA4H(hRM);
			break;
		}
		// BGA1 L1
		case 5: {
			int lRM = 1;
			BGA4L(lRM);
			break;
		}
		// BGA2 L1
		case 6: {
			int lRM = 1;
			BGA4L(lRM);
			break;
		}
		// BGA1 L2
		case 7: {
			int lRM = 2;
			BGA4L(lRM);
			break;
		}
		// BGA2 L2
		case 8: {
			int lRM = 2;
			BGA4L(lRM);
			break;
		}
		}
		MPI.Finalize();
	}

	private static void BGA4L(int lRM) {
		int RECVSize = 2;
		int generations = 10;
		Population<Alignment> pop = new Population<Alignment>();
		// 接收初始种群
		Population<Alignment> popArray[] = new Population[1];
		while(true){
			pop.setIndividualsNum(20);
			long currentTime = System.currentTimeMillis();
			long timeInterval = currentTime - Constant.startTime;
			// generations depend on internal or external interval
			generations = timeInterval%120>20 ? 50:10;
			MPI.COMM_WORLD.Recv(popArray, 0, 1, MPI.OBJECT, lRM, 98);
			pop = popArray[0];
			Alignment[] ins = aligner.doGA(pop, generations, RECVSize);
			MPI.COMM_WORLD.Send(ins, 0, RECVSize, MPI.OBJECT, lRM, 99);
		}
	}

	private static void BGA4H(int hRM) {
		int RECVSize = 2;
		int generations = 10;
		Population<Alignment> pop = new Population<Alignment>();
		// 接收初始种群
		Population<Alignment> popArray[] = new Population[1];
		MPI.COMM_WORLD.Recv(popArray, 0, 1, MPI.OBJECT, hRM, 98);
		pop = popArray[0];
		// GA，返回给HRM
		Alignment[] ins = aligner.doGA(pop, generations, RECVSize);
		MPI.COMM_WORLD.Send(ins, 0, RECVSize, MPI.OBJECT, hRM, 99);
		// 接收个体
		while(true){
			pop.setIndividualsNum(20);
			long currentTime = System.currentTimeMillis();
			long timeInterval = currentTime - Constant.startTime;
			// generations depend on internal or external interval
			generations = timeInterval%120>20 ? 50:10;
			MPI.COMM_WORLD.Recv(ins, 0, RECVSize, MPI.OBJECT, hRM, 99);
			pop.addIndividuals(ins);
			ins = aligner.doGA(pop, generations, RECVSize);
			MPI.COMM_WORLD.Send(ins, 0, RECVSize, MPI.OBJECT, hRM, 99);
		}
	}
}
