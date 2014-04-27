package com.darran.msga.operator;

import com.darran.msga.population.Alignment;

/**
   Mutation that does nothing.
   
   @author Tom Austin and Amie Radenbaugh
 */
public class NoopMutation extends Mutation
{
	/**
	   Performs the mutation
	   @see gamsa.operator.Mutation#perform(gamsa.population.Alignment)
	 */
	@Override
	public Alignment perform(Alignment parent)
	{
		return parent;
	}
}
