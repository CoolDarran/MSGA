package com.darran.msga.main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import mpi.MPI;

import org.apache.commons.lang3.ArrayUtils;

import com.darran.msga.constant.Constant;
import com.darran.msga.population.Alignment;
import com.darran.msga.population.Population;

/**
 * LRM for MBSA
 * 处理部分个体，将最好的结果提交到HRM里进行处理
 * 
 * @author danran
 * Created on 2014年4月22日
 */

/**
 * 伪代码
1. begin
2. receive_input_sequences(seqs);
3. pop=generate_initial_population(seqs);
4. while termination_condition_not_attained do
5. 	begin
6. 		if external_migration_interval
7. 			begin
8. 				migrate(individuals, LRMs);
9. 				migrate(individuals, HRM);
10. 			receive(individuals, LRMs);
11. 			receive(individuals, HRM);
12. 			select_best_individuals();
13. 			pop=generate_new_population();
14. 		end
15. 	cut_individuals_into_parts();
16. 	migrate_parts(BGAs);
17. 	receive_best_parts(BGAs);
18. 	merge_parts();
19. 	apply_gap_operators();
20. 	evaluate_individuals();
21. 	generate_new_population();
22. end
23.end
 */
public class LRM {
	// 传递到下一个种群的个体个数
	private static final int NEWSize = 5;
	// 接收的个体个数
	private static final int RECVSize = 2;
	// 发送到HRM的个体个数
	private static final int MIGRSize = 2;
	// 输入序列
	private String[] inputSequences;
	// 用来存储接收到的个体
	private Population<Alignment> pop;
	// 在BGA和HRM之间接收和发送的个体
	private Alignment individualList[];
	//
	private MultiSeqAligner aligner = new MultiSeqAligner();
	// 存储BGA、HRM及LRMs的rank号
	private int[] BGAs;
	private int HRM;
	private int rlrm; // right lrm
	private int llrm; // left lrm		
	
	public LRM(int[] bGAs, int hRM, int rlrm, int llrm) {
		super();
		BGAs = bGAs;
		HRM = hRM;
		this.rlrm = rlrm;
		this.llrm = llrm;
	}
	
	public void mainDo(String[] inputSeqs){
		recvInputSeq(inputSeqs);
		genInitPopulation();
		Constant.setStartTime(System.currentTimeMillis());
		do{
			long currentTime = System.currentTimeMillis();
			long timeInterval = currentTime - Constant.startTime;
			// external migration interval
			if(timeInterval%120>20){
				// migrate to right lrm
				migrateIndividual(MIGRSize, rlrm);
				// migrate to hrm
				migrateIndividual(MIGRSize, HRM);
				// receive from left lrm
				receive(RECVSize, llrm);
				// receive from hrm
				receive(RECVSize, HRM);
				// using the bests to generate the new population
				genNewPoplation(bests(NEWSize));
			}
			// migrate individual parted to BGAs
			migrateParts2BGAs(cutIndividualsIntoParts());
			// receive individual part from BGAs
			receiveBGAParts(RECVSize, BGAs);
			// merge the individual received from BGAs
			mergeBGAParts(RECVSize, BGAs);
			// gap operation
			gapOperators();
			// using the bests to generate the new population
			genNewPoplation(bests(NEWSize));
		}while (true);
	}

	public void recvInputSeq(String[] inputSeqs){
		this.inputSequences = inputSeqs;
	}
	
	public void genInitPopulation(){
		pop = aligner.generateInitialPopulation(this.inputSequences,180);
	}
	
	/**
	 * cut individuals into parts
	 * @return
	 */
	public List<String[]> cutIndividualsIntoParts(){
		double part = BGAs.length;
		List<String[]> partSeqs = new ArrayList<String[]>();
		// cut the sequences into n part
		for(int i = 0; i < part; i++){
			String[] temp = new String[inputSequences.length];
			int j = 0;
			for(String seq : inputSequences){
				int partLen = (int) Math.floor(seq.length()/part);
				String reg = "(?<=\\G.{"+ partLen +"})";
				temp[j] = seq.split(reg)[i];
				j++;
			}
			partSeqs.add(temp);
		}
		return partSeqs;
	}
	
	/**
	 * migrate individuals parted to BGAs
	 * @param partSeqs
	 */
	public void migrateParts2BGAs(List<String[]> partSeqs){
		@SuppressWarnings("unchecked")
		Population<Alignment> popArray[] = new Population[partSeqs.size()];
		int i = 0;
		for(String[] seqs : partSeqs){
			Population<Alignment> pop = aligner.generateInitialPopulation(seqs);
			popArray[i] = pop;
			MPI.COMM_WORLD.Send(popArray, i, 1, MPI.OBJECT, BGAs[i], 98);
			i++;
		}
	}
	
	/**
	 * receive BGA parts
	 * @param recvSize 
	 * @param bGAs
	 */
	public void receiveBGAParts(int recvSize, int[] bGAs){
		individualList = new Alignment[recvSize*bGAs.length];
		for(int i = 0; i<bGAs.length; i++){
			MPI.COMM_WORLD.Recv(individualList, i*recvSize, recvSize, MPI.OBJECT, bGAs[i], 99);
		}
	}
	
	/**
	 * merge parts received from BGAs
	 * @param recvSize 
	 * @param bGAs 
	 */
	public void mergeBGAParts(int recvSize, int[] bGAs){
		String[] merged = null;
		List<Alignment> as = new ArrayList<Alignment>();
		for(int i=0;i<recvSize;i++){
			for(int j=0;j<bGAs.length;j++)
				merged = ArrayUtils.addAll(merged, individualList[i+j*recvSize].getSequences());
			Alignment in = new Alignment(merged);
			as.add(in);
		}
		individualList = new Alignment[as.size()];
		int n = 0;
		for(Alignment a : as){
			individualList[n] = a;
			n++;
		}
			
	}
	
	/**
	 * add gap into the new individual merged
	 */
	public void gapOperators(){
		
	}
	
	/**
	 * select the best n individuals including origin and received
	 * @param newSize 
	 * @return 
	 */
	public Alignment[] bests(int newSize){
		List<Alignment> origin = pop.getIndividuals();
		Alignment[] temp = new Alignment[origin.size()];
		for(int i = 0;i<origin.size();i++)
			temp[i] = origin.get(i);
		individualList = ArrayUtils.addAll(individualList, temp);
		Arrays.sort(individualList,Collections.reverseOrder());
		Alignment[] temp1 = new Alignment[newSize];
		if(individualList.length >= newSize){
			for(int i = 0;i<newSize;i++)
				temp1[i] = individualList[i];
			return temp1;
		}
		else
			return individualList;
		
	}
	
	/**
	 * generate a new population
	 * by adding the best into current population
	 * @param alignments 
	 */
	public void genNewPoplation(Alignment[] alignments){
		pop = new Population<Alignment>();
		pop.addIndividuals(alignments);
		pop.setIndividualsNum(180);
		// 重新设置individual list
		individualList = alignments;
	}

	/**
	 * receive individuals from left lrm or hrm
	 * 
	 * @param recvsize
	 * @param id
	 */
	private void receive(int recvsize, int id) {
		Alignment[] temp = new Alignment[recvsize];
		MPI.COMM_WORLD.Recv(temp, 0, recvsize, MPI.OBJECT, id, 99);
		individualList = ArrayUtils.addAll(individualList, temp);
	}

	/**
	 * migrate individual to right lrm or hrm
	 * 
	 * @param migrsize
	 * @param id
	 */
	private void migrateIndividual(int migrsize, int id) {
		Arrays.sort(individualList,Collections.reverseOrder());
		MPI.COMM_WORLD.Isend(individualList, 0, migrsize, MPI.OBJECT, id, 99);
		//MPI.COMM_WORLD.Send(individualList, 0, migrsize, MPI.OBJECT, id, 99);
	}
}
