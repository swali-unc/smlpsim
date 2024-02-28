#pragma once

struct request;

struct task {
	double cost;
	double period;
	double deadline;
	request* r;
};

struct job {
	double start;
	double remCost;
	double deadline;
	double totalcost;
	int cpu;
	request* r;
};

struct request {
	double* Li;
	int SMcount;
	double reqsize;
	job* j;
	double issuetime;
	double piblockingtime;
	int k;
};

enum acttype {
	RELEASE, UNLOCK, COMPLETE, SLICE
};

struct act {
	acttype type;
	double time;
	union {
		job* j;
		request* r;
	};
};