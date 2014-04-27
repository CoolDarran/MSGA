package com.darran.msga.main;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import mpi.MPI;

import com.darran.msga.population.Alignment;
import com.darran.msga.population.Individual;
import com.darran.msga.population.Population;

/**
 * HRM for MBSA
 * �������и��壬���Ľ���ɴ˲���
 * @author danran
 * Created on 2014��4��22��
 */
/**
 * α����
1. begin
2. receive_input_sequences(seqs);
3. pop=generate_initial_population(seqs);
4. distribute_initial_population(pop, BGAs);
5. while termination_condition_not_reached do
6. begin
7. 		receive(individuals, BGAs, i);
8. 		migrate_individuals (i, LRMs);
9. 		receive(individuals, LRMs, i);
10. 	distribute_individuals (i, BGAs);
11.end
12.receive(individuals, BGAs, i);
13.output = obtain_best_individual();
14.end
 */
public class HRM {
	// ���յĸ������
	private static final int RECVSize = 2;
	// ���͵�LRM�ĸ������
	private static final int MIGRSize = 2;
	// �ַ���BGA�ĸ������
	private static final int DESRSize = 2;
	// ��������
	private String[] inputSequences;
	// �����洢���յ��ĸ���
	private Population<Alignment> pop;
	// ��BGA��LRM֮����պͷ��͵ĸ���
	private Alignment individualList[];
	// 
	private MultiSeqAligner aligner = new MultiSeqAligner();
	// �洢BGA��LRM��rank��
	private int[] BGAs;	
	private int[] LRMs;

	public HRM(int[] bGAs, int[] lRMs) {
		BGAs = bGAs;
		LRMs = lRMs;
	}

	public void recvInputSeq(String[] inputSeqs){
		this.inputSequences = inputSeqs;
	}
	
	public void genInitPopulation(){
		pop = aligner.generateInitialPopulation(this.inputSequences,180);
	}
	
	public void distributePoplation2BGA(){
		int size = this.BGAs.length;
		double individualSize = this.pop.getPopulationSize();
		int bgaSize = (int) Math.floor(individualSize/size);
		Population<Alignment> popArray[] = new Population[size];
		List<Alignment> individualList = pop.getIndividuals();
		for(int i=0;i<size;i++){
			Population<Alignment> temp = new Population<Alignment>();
			if(	individualSize%bgaSize != 0 && i == size -1)
				temp.addIndividuals(individualList.subList(i*bgaSize, individualList.size()-1));
			else
				temp.addIndividuals(individualList.subList(i*bgaSize, (i+1)*bgaSize-1));
			popArray[i] = temp;
		}
		for(int i=0;i<size;i++)
			MPI.COMM_WORLD.Send(popArray, i, 1, MPI.OBJECT, BGAs[i], 98);
	}
	
	public Individual mainDo(String[] inputSeqs){
		Individual best = null;
		recvInputSeq(inputSeqs);
		genInitPopulation();
		distributePoplation2BGA();
		do{
			// receive best i individuals from BGAs
			receive(BGAs,RECVSize);
			// migrates the best i individuals to LRMs
			migrateIndividuals (MIGRSize, LRMs);
			// receive best i individuals from LRMs
			receive(LRMs,RECVSize);
			// distributes received individuals from LRMs to BGAs
			distributeIndividuals2BGA(DESRSize,BGAs);
		}while (!aligner.hasSolution(pop));
		// receive best i individuals from BGAs for the last time
		receive(BGAs,RECVSize);
		// select the best individual
		best = pop.best();
		return best;
	}

	/**
	 * distribute individuals to BGA
	 * @param dstrSize
	 * @param bGAs
	 */
	private void distributeIndividuals2BGA(int dstrSize, int[] bGAs) {
		Arrays.sort(individualList,Collections.reverseOrder());
		for(int i = 0; i<bGAs.length;i++){
			MPI.COMM_WORLD.Send(individualList, i*dstrSize, dstrSize, MPI.OBJECT, bGAs[i], 99);
		}
	}

	/**
	 * migrate individuals to LRMs
	 * @param migrSize
	 * @param lRMs
	 */
	private void migrateIndividuals(int migrSize, int[] lRMs) {
		Arrays.sort(individualList,Collections.reverseOrder());
		for(int i = 0; i<lRMs.length;i++){
			MPI.COMM_WORLD.Send(individualList, i*migrSize, migrSize, MPI.OBJECT, lRMs[i], 99);
		}
	}

	/**
	 * receive individuals from BGAs or LRMs
	 * 
	 * @param individualList
	 * @param ids
	 * @param recvSize
	 */
	private void receive(int[] ids, int recvSize) {
		individualList = new Alignment[recvSize*ids.length];
		for(int i = 0; i<ids.length; i++){
			MPI.COMM_WORLD.Recv(individualList, i*recvSize, recvSize, MPI.OBJECT, ids[i], 99);
		}
		pop.addIndividuals(individualList);
		pop.setIndividualsNum(180);
	}
}
