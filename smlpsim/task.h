#pragma once

struct task {
	double cost;
	double period;
	double deadline;
};

struct job {
	double start;
	double remCost;
	double deadline;
	double totalcost;
	int cpu;
};

struct request {
	double* Li;
	int SMcount;
	double reqsize;
	job* j;
	double issuetime;
};

enum acttype {
	RELEASE, UNLOCK, COMPLETE
};

struct act {
	acttype type;
	double time;
	union {
		job* j;
		request* r;
	};
};