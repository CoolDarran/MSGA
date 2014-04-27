package com.darran.msga.main;

import java.util.Collections;
import java.util.List;
import java.util.Properties;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.darran.msga.operator.BlockShuffleMutation;
import com.darran.msga.operator.Crossover;
import com.darran.msga.operator.GapColumnDeletionMutation;
import com.darran.msga.operator.GapDeletionMutation;
import com.darran.msga.operator.GapInsertionMutation;
import com.darran.msga.operator.Mutation;
import com.darran.msga.operator.NoopCrossover;
import com.darran.msga.operator.NoopMutation;
import com.darran.msga.operator.OnePointCrossoverGapsBeginning;
import com.darran.msga.operator.OnePointCrossoverGapsEnd;
import com.darran.msga.operator.OnePointCrossoverGapsMiddle;
import com.darran.msga.operator.Operator;
import com.darran.msga.population.Alignment;
import com.darran.msga.population.Individual;
import com.darran.msga.population.Population;
import com.darran.msga.scorer.Scorer;

/**
 * Basic engine to run the genetic algorithm.
 * @author L.danran
 * Created on Apr 22, 2014
 */
public class MultiSeqAligner
{
	private static final Level LOG_LEVEL = Level.INFO;	
	private static final boolean DEBUG = false;	
	
	private static Logger logger = Logger.getLogger(MultiSeqAligner.class.getName());
	
	private Population<Crossover> i_crossoverOps;
	private Population<Mutation> i_mutationOps;
	
	//Used if crossovers and mutations are lumped together.
	private Population<Operator> i_operators;
	
	//Settings for how sequences are processed.
	private boolean i_mergeOperators;
	private int i_populationSize;
	private double i_percentageIncrease;
	private int i_unchangedRoundsNeeded;
	
	//Variables to check for termination condition
	private long i_bestLastRound;
	private int i_numTimesUnchanged;
	private int i_maxRounds;
	private int i_numRounds;
	
	/**
	   Constructor.  Default arguments will be used.
	 */
	public MultiSeqAligner()
	{
		this(getDefaultProperties());
	}
	
	/**
	   Specifies default properties, in case there is no config.
	 */
	protected static Properties getDefaultProperties()
	{
		Properties p = new Properties();
		
		p.put("mergeOps", "false");
		p.put("seqScorer", "com.darran.msga.scorer.DNAScorer");
		//p.put("seqScorer", "gamsa.scorer.Blosum62Scorer");
		p.put("percentageIncrease", "15.0");
		p.put("populationSize", "50");
		p.put("unchangedRoundsNeeded", "50");
		p.put("maxRounds", "200");
		
		return p;
	}
	
	/**
	 * Returns the number of rounds that were executed.
	 * @return
	 */
	protected int getNumberOfRoundsExecuted()
	{
		return i_numRounds;
	}
	
	/**
	   Constructor.  Properties will override default values.
	 */
	public MultiSeqAligner(Properties p)
	{
		//Set the log level
		for (Handler h : Logger.getLogger("").getHandlers())
		{
			h.setLevel(LOG_LEVEL);
		}
		logger.setLevel(LOG_LEVEL);
		
		//Specify sequence scorer to be used.
		String seqScorerClassName = p.getProperty("seqScorer");
		try
		{
			Class c = Class.forName(seqScorerClassName);
			Scorer scorer = (Scorer) c.newInstance();
			InfoCenter.getCenter().setScorer(scorer);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			logger.severe("Error loading class '" + seqScorerClassName + "'.  Exiting.");
			System.exit(1);
		}
		
		//Specify the size of the population and the number of rounds with no change in score
		i_populationSize = Integer.parseInt(p.getProperty("populationSize"));
		i_unchangedRoundsNeeded = Integer.parseInt(p.getProperty("unchangedRoundsNeeded"));
		
		//A negative value means that maxRounds will be disabled.
		i_maxRounds = -999;
		if (p.containsKey("maxRounds"))
			i_maxRounds = Integer.parseInt(p.getProperty("maxRounds"));
		
		//Specify the amount of spacing that should be in the alignments.
		i_percentageIncrease = Double.parseDouble(p.getProperty("percentageIncrease"));
		
		//Determine if crossovers and mutations should be handled in separate steps.
		i_mergeOperators = Boolean.parseBoolean(p.getProperty("mergeOps"));
		
		//Load operators...Nice to specify these from a config too.
		i_crossoverOps = new Population<Crossover>();
		i_crossoverOps.addIndividual(new NoopCrossover());
		i_crossoverOps.addIndividual(new OnePointCrossoverGapsBeginning());
		i_crossoverOps.addIndividual(new OnePointCrossoverGapsMiddle());
		i_crossoverOps.addIndividual(new OnePointCrossoverGapsEnd());
		
		i_mutationOps = new Population<Mutation>();
		i_mutationOps.addIndividual(new NoopMutation());
		i_mutationOps.addIndividual(new BlockShuffleMutation());
		i_mutationOps.addIndividual(new GapInsertionMutation());
		i_mutationOps.addIndividual(new GapDeletionMutation());
		i_mutationOps.addIndividual(new GapColumnDeletionMutation());		
		
		if (i_mergeOperators)
		{
			i_operators = new Population<Operator>();
			i_operators.addAllIndividuals(i_crossoverOps, i_mutationOps);
		}
	}
	
	/**
	   Returns the Optimal solution.
	 */
	public Individual findSolution(String[] inputSequences) 
	{
		logger.info("Starting search for solution.");
		
		i_bestLastRound = Long.MIN_VALUE;
		i_numTimesUnchanged=0;
		i_numRounds = 0;
		double bestScore = Double.NEGATIVE_INFINITY;
		double tempBest;
		
		Population<Alignment> pop = this.generateInitialPopulation(inputSequences);
		logger.fine("Initial population:");
		for (Alignment a : pop)
			logger.fine(a.toString());
		
		int counter = 0;
		Alignment best;
		do
		{
			best = pop.best();
			
			pop = this.performReproduction(pop);
			
			if (i_mergeOperators)
			{
				pop = this.performOperators(pop);
			}
			else
			{
				pop = this.performCrossover(pop);
				pop = this.performMutation(pop);
			}
			
			//Preserve the last best individual for another generation.
			pop.addIndividual(best);
			
			// keep track of the overall best score
			tempBest = best.getFitness();
			if (Double.compare(tempBest, bestScore) > 0)
				bestScore = tempBest;
			
			if (counter%10 == 0)
			{
				logger.info("Number of Rounds: " + i_numRounds);
				logger.info("Best so far in population: " + pop.best().getFitness());
				logger.info("Best ever: " + bestScore);
				logger.info("" + best);
				logger.fine("Crossover Operators: " + i_crossoverOps.getDistributionString());
				logger.fine("Mutation Operators: " + i_mutationOps.getDistributionString());
			}
			counter++;
		} while (!hasSolution(pop));
		
		logger.info("Found solution");
		logger.info("Number of Rounds: " + i_numRounds);
		logger.info("Best ever: " + bestScore);
		logger.info("Crossover Operators: " + i_crossoverOps.getDistributionString());
		logger.info("Mutation Operators: " + i_mutationOps.getDistributionString());
		
		// Return the best scoring individual in the population.
		return pop.best();
	}

	/**
	   Generates a population of Alignments.
	 */
	public Population<Alignment> generateInitialPopulation(String[] inputSequences)
	{
		Population<Alignment> pop = new Population<Alignment>();
		if (inputSequences != null)
		{
			int maxSequenceSize = 0;
			
			//find the max sequence size
			for (int i=0; i<inputSequences.length; i++)
			{
				if (inputSequences[i].length() > maxSequenceSize)
					maxSequenceSize = inputSequences[i].length(); 
			}
			
			// create a new sequence length based on the given percentage
			int alignmentLength = 0;
			if (Math.round(i_percentageIncrease) == 0)
				alignmentLength = maxSequenceSize;
			else
				alignmentLength = new Double(maxSequenceSize * (1 + i_percentageIncrease * .01)).intValue();
		
			// create the individuals
			for (int i=0; i<i_populationSize; i++)
			{
				pop.addIndividual(new Alignment(inputSequences, alignmentLength));
			}
		}
		return pop;
	}
	
	/**
	   Generates a population of Alignments.
	 */
	public Population<Alignment> generateInitialPopulation(String[] inputSequences, int popSize)
	{
		i_populationSize = popSize;
		Population<Alignment> pop = new Population<Alignment>();
		if (inputSequences != null)
		{
			int maxSequenceSize = 0;
			
			//find the max sequence size
			for (int i=0; i<inputSequences.length; i++)
			{
				if (inputSequences[i].length() > maxSequenceSize)
					maxSequenceSize = inputSequences[i].length(); 
			}
			
			// create a new sequence length based on the given percentage
			int alignmentLength = 0;
			if (Math.round(i_percentageIncrease) == 0)
				alignmentLength = maxSequenceSize;
			else
				alignmentLength = new Double(maxSequenceSize * (1 + i_percentageIncrease * .01)).intValue();
		
			// create the individuals
			for (int i=0; i<i_populationSize; i++)
			{
				pop.addIndividual(new Alignment(inputSequences, alignmentLength));
			}
		}
		return pop;
	}
	
	/**
	   Creates and returns the new generation of alignments, based on chance
	   and the fitness values of the old alignments.  One spot is left for
	   the best individual after all transformations.
	 */
	private Population<Alignment> performReproduction(Population<Alignment> oldPop)
	{
		Population<Alignment> newPop = new Population<Alignment>();
		Alignment alignment;
		Alignment clonedAlignment;
		
		//The new population will be the same size as the old.
		for (int i=0; i<i_populationSize-1; i++)
		{
			alignment = oldPop.getRandomIndividual();
			clonedAlignment = alignment.clone();

			newPop.addIndividual(clonedAlignment);
		}
		oldPop = null;
		return newPop;
	}
	
	/**
	   Perform Crossover operation, and score the results.
	   Return a new, post-crossover population.
	 */
	private Population<Alignment> performCrossover(Population<Alignment> oldPop)
	{
		Population<Alignment> newPop = new Population<Alignment>();
		Alignment mother;
		Alignment father = null;
		int tries = 0;
		Crossover xover;
		Alignment[] children;
		
		// The next population size should be the same size as the current generation.
		while (newPop.getPopulationSize() < oldPop.getPopulationSize())
		{
			//Choose a mother and father.
			mother = oldPop.getRandomIndividual();
			
			//The father and mother should not be the same.			
			tries = 0;
			do
			{
				father = oldPop.getRandomIndividual();
			} while (tries++<10 && father == mother);
			
			xover = i_crossoverOps.getRandomIndividual();
			//logger.info("Before XOver: " + xover.getClass().getName());
			children = xover.perform(mother, father);
			//logger.info("After XOver: " + xover.getClass().getName());
			
			if (DEBUG)
			{
				for (int i=0; i<children.length; i++)
				{
					if (!children[i].verifyAlignment(mother))
					{
						logger.warning("Child does not match parent after " + xover);
						logger.warning(mother.toString());
						logger.warning(children[i].toString());
						//replace bad child with mother
						children[i] = mother;
					}
				}
			}

			newPop.addIndividual(children[0]);
			newPop.addIndividual(children[1]);
			
			//Update the success of this operator, and past operators.
			rateOperators(xover, children, mother, father);
			i_crossoverOps.recalculateProbabilities();
			logger.finer(i_crossoverOps.getDistributionString());
		}
		oldPop = null;
		return newPop;
	}
	
	/**
	   Perform mutation and score the results.
	   Return a new, post-mutation population.
	 */
	private Population<Alignment> performMutation(Population<Alignment> oldPop)
	{
		Population<Alignment> newPop = new Population<Alignment>();
		Alignment parent;
		Mutation mute;
		Alignment child;
		Alignment[] children = new Alignment[1];
		
		// The next population size should be the same size as the current generation.
		while (newPop.getPopulationSize() < oldPop.getPopulationSize())
		{
			//Choose a parent.
			parent = oldPop.getRandomIndividual();			
			
			mute = i_mutationOps.getRandomIndividual();
			//logger.info("Before Mutation: " + mute.getClass().getName());
			child = mute.perform(parent);
			//logger.info("After Mutation: " + mute.getClass().getName());
			
			if (DEBUG)
			{
				if (!child.verifyAlignment(parent))
				{
					logger.warning("Child does not match parent after " + mute);
					logger.warning(parent.toString());
					logger.warning(child.toString());
					//replace bad child with parent
					child = parent;
				}
			}
			
			newPop.addIndividual(child);
			
			//Update the success of this operator, and past operators.
			children[0] = child;
			rateOperators(mute, children, parent);
			i_mutationOps.recalculateProbabilities();
			logger.finer(i_mutationOps.getDistributionString());
		}
		oldPop = null;
		return newPop;
	}
	
	/**
	   This is used if mutations and crossovers should be combined together.
	 */
	private Population<Alignment> performOperators(Population<Alignment> oldPop)
	{
		Population<Alignment> newPop = new Population<Alignment>();
		Operator op;
		Alignment mother;
		Alignment father;
		Alignment[] children = null;
		Alignment[] parents = null;
		int tries = 0;
		
		// The next population size should be the same size as the current generation.
		while (newPop.getPopulationSize() < oldPop.getPopulationSize())
		{
			op = i_operators.getRandomIndividual();
			
			//Both mutations and crossovers need at least one parent
			mother = oldPop.getRandomIndividual();			
			
			if (op instanceof Crossover)
			{
				//The father and mother should not be the same.				
				tries=0;
				do
				{
					father = oldPop.getRandomIndividual();
				} while (tries++<10 && father == mother);
				
				children = op.performOp(mother, father);
				parents = new Alignment[2];
				parents[0] = mother;
				parents[1] = father;
			}
			else
			{
				children = op.performOp(mother);
				parents = new Alignment[1];
				parents[0] = mother;
			}

			for (Alignment child : children)
			{
				if (DEBUG && !child.verifyAlignment(parents[0]))
				{
					logger.warning("Child does not match parent after " + op);
					logger.warning(parents[0].toString());
					logger.warning(child.toString());
					//Replace bad child with parent
					child = parents[0];
				}
				newPop.addIndividual(child);
			}
			
			//Update the success of this operator, and past operators.
			rateOperators(op, children, parents);
			i_operators.recalculateProbabilities();
			logger.finer(i_operators.getDistributionString());
		}
		
		return newPop;
	}

	/**
	   Rate the success of the operation, plus operations leading up
	   to this operation.
	 */
	private void rateOperators(Operator op, Alignment[] children, Alignment... parents)
	{
		logger.finest(op.toString());
		
		//Rate the children
		double childrenFitness = 0;
		for (Alignment child : children) 
			childrenFitness += child.getFitness();
		
		//Rate the parents
		double parentFitness = 0;
		for (Alignment parent : parents) 
			parentFitness += parent.getFitness();
		
		//The success of the operation is based on the difference between the children's
		// fitness and the parents' fitness.
		double opSuccess = childrenFitness - parentFitness;
		
		op.addScore(opSuccess);
	}
	
	/**
	   Returns true if we have arrived at a solution.
	   This is based on the best scoring individual in
	   the population.  If this has not changed for
	   a certain number of generations, we assume that
	   the population has reached a stable point.
	 */
	public boolean hasSolution(Population<Alignment> p)
    {
//		if (i_maxRounds > 0)
//		{
//			logger.fine("NumRounds: " + i_numRounds + "/" + i_maxRounds);
//			if (i_numRounds >= i_maxRounds)
//				return true;
//			else
//				i_numRounds++;
//		}
//		else
//		{
//			i_numRounds++;
//		}
		
        Alignment newBest = p.best();
        
//		if (Math.round(newBest.getFitness()) == i_bestLastRound)
//            i_numTimesUnchanged++;
//        else
//        {
//            i_numTimesUnchanged = 0;
//            i_bestLastRound = Math.round(newBest.getFitness());
//        }
//        
//		if (i_unchangedRoundsNeeded > i_numTimesUnchanged)
//            return false;
//        else
//            return true;
		if (newBest.getFitness() < 349)
            return false;
        else
            return true;
    }
	
	/**
	   Returns the Optimal solution.
	 */
	public Individual findLastSolution(List<Alignment> aliList) 
	{
		logger.info("Starting search for solution.");
		
		i_bestLastRound = Long.MIN_VALUE;
		i_numTimesUnchanged=0;
		i_numRounds = 0;
		double bestScore = Double.NEGATIVE_INFINITY;
		double tempBest;
		
		Population<Alignment> pop = new Population<Alignment>();
		for(Alignment al : aliList)
			pop.addIndividual(al);
		logger.fine("Initial population:");
		for (Alignment a : pop)
			logger.fine(a.toString());
		
		int counter = 0;
		Alignment best;
		do
		{
			best = pop.best();
			
			pop = this.performReproduction(pop);
			
			if (i_mergeOperators)
			{
				pop = this.performOperators(pop);
			}
			else
			{
				pop = this.performCrossover(pop);
				pop = this.performMutation(pop);
			}
			
			//Preserve the last best individual for another generation.
			pop.addIndividual(best);
			
			// keep track of the overall best score
			tempBest = best.getFitness();
			if (Double.compare(tempBest, bestScore) > 0)
				bestScore = tempBest;
			
			if (counter%10 == 0)
			{
				logger.info("Number of Rounds: " + i_numRounds);
				logger.info("Best so far in population: " + pop.best().getFitness());
				logger.info("Best ever: " + bestScore);
				logger.info("" + best);
				logger.fine("Crossover Operators: " + i_crossoverOps.getDistributionString());
				logger.fine("Mutation Operators: " + i_mutationOps.getDistributionString());
			}
			counter++;
		} while (!hasSolution(pop));
		
		logger.info("Found solution");
		logger.info("Number of Rounds: " + i_numRounds);
		logger.info("Best ever: " + bestScore);
		logger.info("Crossover Operators: " + i_crossoverOps.getDistributionString());
		logger.info("Mutation Operators: " + i_mutationOps.getDistributionString());
		
		// Return the best scoring individual in the population.
		return pop.best();
	}
	
	/**
	 * do ga
	 * @param pop
	 * @param generations
	 * @param size 
	 * @return
	 */
	public Alignment[] doGA(Population<Alignment> pop, int generations, int size) 
	{
		
		i_bestLastRound = Long.MIN_VALUE;
		i_numTimesUnchanged=0;
		i_numRounds = 0;
		double bestScore = Double.NEGATIVE_INFINITY;
		double tempBest;
		// 记录操作次数
		int counter = 0;
		Alignment best;
		do
		{
			best = pop.best();
			// 从旧种群里根据概率和分数来选取生成新的种群，分数越好的个体越有可能被选取
			pop = this.performReproduction(pop);
			if (i_mergeOperators)
				// 合并交叉变异操作
				pop = this.performOperators(pop);
			else
			{
				// 交叉
				pop = this.performCrossover(pop);
				// 变异
				pop = this.performMutation(pop);
			}
			// 将最好的个体保存下来，进入到下一次种群生成
			pop.addIndividual(best);
			// 记录整个种群里最好的分数
			tempBest = best.getFitness();
			if (Double.compare(tempBest, bestScore) > 0)
				bestScore = tempBest;
			if (counter%10 == 0)
			{
				logger.info("Number of Rounds: " + i_numRounds);
				logger.info("Best so far in population: " + pop.best().getFitness());
				logger.info("Best ever: " + bestScore);
				logger.info("" + best);
				logger.fine("Crossover Operators: " + i_crossoverOps.getDistributionString());
				logger.fine("Mutation Operators: " + i_mutationOps.getDistributionString());
			}
			counter++;
		} while (counter < generations); //当操作次数大于预定义的次数时停止遗传操作
		// 返回最好的n个个体
		List<Alignment> ins = pop.getIndividuals();
		Collections.sort(ins,Collections.reverseOrder());
		Alignment[] alis = new Alignment[size];
		for(int i = 0;i < size; i++)
			alis[i] = ins.get(i);
		return alis;
	}
	
}