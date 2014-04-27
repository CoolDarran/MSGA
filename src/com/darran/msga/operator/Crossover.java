package com.darran.msga.operator;
import com.darran.msga.population.Alignment;

/**
 * Represents a crossover operation between two parent
   alignments, resulting in two new alignments.
 * @author danran
 * Created on 2014Äê4ÔÂ26ÈÕ
 */
public abstract class Crossover extends Operator
{
	/**
	   Takes 2 parent alignments and returns two child alignments.
	 */
    public abstract Alignment[] perform(Alignment mother, Alignment father);
    
    /**
       Wrapper around perform method.
       @see gamsa.operator.Operator#performOp(gamsa.population.Alignment[])
     */
    public Alignment[] performOp(Alignment... parents)
    {
    		if (parents.length != 2)
    			throw new IllegalArgumentException("Crossovers can only have two parents");
    		
    		return perform(parents[0], parents[1]);
    }
}
