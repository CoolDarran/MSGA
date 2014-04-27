package com.darran.msga.constant;

/**
 * 记录全局开始时间
 * @author danran
 * Created on 2014年4月26日
 */
public class Constant {
	public static long startTime;

	public static long getStartTime() {
		return startTime;
	}

	public static void setStartTime(long startTime) {
		Constant.startTime = startTime;
	}
}
