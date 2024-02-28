#include <stdio.h>
#include <stdlib.h>

#include <random>
#include <vector>
#include <numeric>
#include <utility>
#include <queue>

#include "task.h"

using std::default_random_engine;
using std::uniform_real_distribution;
using std::vector;
using std::priority_queue;
using std::queue;

double plock = 0.1;
double tmin = 10;
double tmax = 100;
double tg = 1;
int n = 100;
const int M = 16;
double U = 0.8;
int h = 2;
int H = 16;
double Omega = 10;
double smmult = 0.1;
double paraloss = 0.006;
double minReqSizeMultiplier = 0.4;

double sample_Ti();
double* generate_utilizations();
request* generate_request(job* j);
int z(request* r, int It);

default_random_engine gen;

void main() {
	gen.seed(time(nullptr));

	double* u = generate_utilizations();
	int hyperperiod;
	vector<task*> taskSet;
	double totalutil = 0;
	for (int i = 0; i < n; ++i) {
		auto Ti = sample_Ti();
		taskSet.push_back(new task{ std::min( M * u[i] * Ti, Ti ), Ti, Ti });
		totalutil += taskSet[i]->cost / taskSet[i]->period;
		//printf("(%.2f,%.2f)\n", taskSet[i]->cost, taskSet[i]->period);
	}
	delete[] u;

	// calculate hyperperiod
	hyperperiod = std::accumulate(taskSet.begin(), taskSet.end(), 1, [](int a, task* b) { return std::lcm(a,(int)b->period); });
	printf("hyperperiod: %d\n", hyperperiod);
	printf("total util: %.2f (want %.2f)\n", totalutil, M * U);

	// simulate
	vector<job*> jobs;
	auto actcmp = [](act* a, act* b) { return a->time > b->time; };
	priority_queue<act*,std::vector<act*>,decltype(actcmp)> events;

	for (auto i : taskSet) {
		for (double t = 0; t < hyperperiod; t += i->period) {
			jobs.push_back(new job{ t, i->cost, t + i->period, i->cost, -1 });
			events.push(new act { RELEASE, t, jobs.back() });
		}
	}

	uniform_real_distribution<double> dist(0, 1);
	int deadline_miss_count = 0;
	double closest_miss = -1;

	job* cpu[M];
	for (int i = 0; i < M; ++i)
		cpu[i] = nullptr;
	auto jobcmp = [](job* a, job* b) { return a->deadline > b->deadline; };
	priority_queue<job*, vector<job*>, decltype(jobcmp)> cpuPQ(jobcmp);

	auto lowestcpu = [&cpu] {
		int lowest = 0;
		for (int i = 0; i < M; ++i) {
			if (cpu[i] == nullptr)
				return i;
			if (cpu[i]->deadline > cpu[lowest]->deadline)
				lowest = i;
		}
		return lowest;
	};

	auto highestjob = [&cpuPQ] {
		if( !cpuPQ.size() )
			return (job*)nullptr;
		auto j = cpuPQ.top();
		cpuPQ.pop();
		return j;
	};

	int It = H;
	auto reqcmp = [It](request* a, request* b) { return a->j->deadline > b->j->deadline; };
	priority_queue<request*, vector<request*>, decltype(reqcmp)> reqPQ(reqcmp);
	queue<request*> reqFQ;
	int numInSQ = 0;

	auto satisfy = [&It, &reqFQ, &reqPQ, &events, &numInSQ](double t) {
		if (reqFQ.size() + numInSQ < M && reqPQ.size()) {
			reqFQ.push(reqPQ.top());
			reqPQ.pop();
		}

		if (It < h) return;

		request* r = nullptr;
		if (reqFQ.size()) {
			r = reqFQ.front();
			reqFQ.pop();
		}
		else if (reqPQ.size()) {
			r = reqPQ.top();
			reqPQ.pop();
		}

		if (!r) return;

		r->SMcount = z(r, It);
		It -= r->SMcount;
		events.push(new act{ UNLOCK, t + r->Li[r->SMcount], (job*)r });
		++numInSQ;
	};

	for (double t = 0; t < hyperperiod; ) {
		if (!events.size())
			break;
		act* a = events.top();
		auto tprev = t;
		bool nopop = false;

		for (auto& i : cpu) {
			if (i && i->remCost < a->time - t) {
				a = new act{ COMPLETE, t + i->remCost, i };
				nopop = true;
			}
		}

		if( !nopop )
			events.pop();
		t = a->time;

		for (auto& i : cpu) {
			if(i) i->remCost -= t - tprev;
		}

		if (a->type == RELEASE) {
			//printf("release job at %.2f\n", t);
			if (dist(gen) < plock) {
				// issue request
				auto r = generate_request(a->j);
				r->issuetime = t;
				reqPQ.push(r);
				satisfy(t);
			}
			else {
				int lowest = lowestcpu();
				if (!cpu[lowest] || cpu[lowest]->deadline > a->j->deadline ) {
					cpu[lowest] = a->j;
					a->j->cpu = lowest;
				}
				else
					cpuPQ.push(a->j);
			}
		}
		else if (a->type == UNLOCK) {
			printf("unlock job at %.2f\n", t);
			--numInSQ;
			It += a->r->SMcount;
			satisfy(t);

			int lowest = lowestcpu();
			if (!cpu[lowest] || cpu[lowest]->deadline > a->r->j->deadline) {
				cpu[lowest] = a->r->j;
				a->r->j->cpu = lowest;
			}
			else
				cpuPQ.push(a->r->j);
			delete[] a->r->Li;
			delete a->r;
		}
		else if (a->type == COMPLETE) {
			if (t > a->j->deadline) {
				printf("deadline missed at %.2f\n", t);
				++deadline_miss_count;
			}
			else {
				//printf("job completed at %.2f\n", t);
				if( closest_miss < 0 || a->j->deadline - t < closest_miss )
					closest_miss = a->j->deadline - t;
			}
			cpu[a->j->cpu] = highestjob();
			if( cpu[a->j->cpu] ) cpu[a->j->cpu]->cpu = a->j->cpu;
			a->j->cpu = -1;
		}
		else {
			printf("unknown event at %.2f\n", t);
		}

		delete a;
	}

	printf("number of deadlines missed: %d\n", deadline_miss_count);
	printf("closest miss: %.2f\n", closest_miss);
	printf("simulation complete\n");

	// clean up
	for (auto i : taskSet)
		delete i;
	for (auto i : jobs)
		delete i;
}

request* generate_request(job* j) {
	request* r = new request;
	//r->Li = new int[H / h];
	r->Li = new double[H+1];
	r->SMcount = 0;
	r->j = j;
	double maxRequestSize = h * j->totalcost / smmult;
	uniform_real_distribution<double> dist(minReqSizeMultiplier * maxRequestSize, maxRequestSize);
	r->reqsize = dist(gen);
	memset(r->Li, 0, sizeof(double) * (H+1));
	for (int i = h; i <= H; ++i)
		r->Li[i] = smmult * ceil(r->reqsize / i) * (1 + (i - h) * paraloss);
	//for( int i = 0; i <= H; ++i )
	//	printf("%.2f ", r->Li[i]);
	//printf("\n");
	//for (int i = 0; i < H / h; ++i)
	//	r->Li[i] = smmult * ceil(r->reqsize / (h*(i + 1))) * (1 + (h*(i + 1) - h) * paraloss);
	return r;
}

int z(request* r, int It) {
	auto cyc = [r](int x) {
		int roundval = (int)ceil(r->reqsize / x);
		return ((roundval + h - 1) / h) * h;
	};

	int smallest = It;
	for (int i = It; i > 0; i -= h) {
		if (cyc(i) <= cyc(smallest))
			smallest = i;
	}

	//printf("zi(t)= %d when I(t)=%d\n", smallest, It);
	return smallest;
}

double* generate_utilizations(int n, double u) {
	double* r = new double[n];

	{
		uniform_real_distribution<double> dist(0, 1);
		for (int i = 1; i < n; ++i)
			r[i] = dist(gen);
	}

	double* s = new double[n + 1];
	s[0] = 0;
	s[n] = u;

	for (int i = n; i > 1; --i)
		s[i - 1] = s[i] * powl(r[i - 1], (double)1 / (i - 1));
	delete[] r;

	double* ui = new double[n];
	for (int i = 0; i < n; ++i)
		ui[i] = s[i + 1] - s[i];

	delete[] s;
	return ui;
}

double sample_Ti(double ri, double tg) {
	return floor(exp(ri) / tg) * tg;
}

double sample_ri(double tmin, double tmax, double tg) {
	uniform_real_distribution<double> dist(log10(tmin), log10(tmax + tg));
	return dist(gen);
}

double sample_ri() {
	return sample_ri(tmin, tmax, tg);
}

double sample_Ti(double ri) {
	return sample_Ti(ri, tg);
}

double sample_Ti() {
	return sample_Ti(sample_ri());
}

double* generate_utilizations() {
	return generate_utilizations(n, U);
}