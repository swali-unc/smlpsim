#include <stdio.h>
#include <stdlib.h>
#include <conio.h>

#include <random>
#include <vector>
#include <numeric>
#include <utility>
#include <queue>
#include <thread>
#include <future>

#include "task.h"

using std::default_random_engine;
using std::uniform_real_distribution;
using std::vector;
using std::priority_queue;
using std::deque;

double plock = 0.25;
double tmin = 10;
double tmax = 100;
double tg = 1;
int n = 100;
int M = 16;
double U = 0.5;
int h = 2;
int H = 16;
double Omega = 10; // set to -1 to disable slicing
double slicemargin = 0.01; // ensure CS ends well before slice
double smmult = 0.05;
double paraloss = 0.006;
//double paraloss = 0;
double minReqSizeMultiplier = 0.05;
default_random_engine gen;

double sample_Ti();
double* generate_utilizations();
request* generate_request(task* j);
int z_omlp(request* r, int It);
int z_smlp(request* r, int It);
int z_k(request* r, int It);
int zmin_omlp(request* r);
int zmin_smlp(request* r);
int zmin_k(request* r);

typedef int (*zfunc)(request*, int);
typedef int (*zminfunc)(request*);
double simulate(int hyperperiod, vector<task*>& taskSet, zfunc z, zminfunc zmin);
vector<task*>* make_ts(double& totalutil, double& Lmax, double& LHmax);

void main() {
	FILE* f = nullptr;
	int tsindex = 0;
	std::string tsize;

	gen.seed(time(nullptr));

	auto runsim = [&](int hyperperiod, vector<task*>& taskSet) {
		auto omlp_ret = std::async([hyperperiod, &taskSet] { return simulate(hyperperiod, taskSet, z_omlp, zmin_omlp); });
		auto smlp_ret = std::async([hyperperiod, &taskSet] { return simulate(hyperperiod, taskSet, z_smlp, zmin_smlp); });
		auto k_ret = std::async([hyperperiod, &taskSet] { return simulate(hyperperiod, taskSet, z_k, zmin_k); });
		double omlp_piblock = omlp_ret.get();
		double smlp_piblock = smlp_ret.get();
		double k_piblock = k_ret.get();

		//printf("SMLP pi-blocking: %.2f\n", smlp_piblock);
		//printf("OMLP pi-blocking: %.2f\n", omlp_piblock);
		//printf("K exclusion pi-blocking: %.2f\n", k_piblock);
		fprintf(f,"%d,%.2f,%d,%d,%s,%d,%d,%.2f,%.2f,%.2f,%.2f\n", tsindex, U, M, n, tsize.c_str(), hyperperiod, H, Omega, omlp_piblock, smlp_piblock, k_piblock);
	};

	int hyperperiod;
	vector<task*>* taskSet;
	double totalutil = 0;
	double Lmax = 0;
	double LHmax = 0;
	auto Tvalues = { std::pair<double,double>{3,33}, std::pair<double,double>{10,100}, std::pair<double,double>{50,200} };
	auto Tnames = {"small", "medium", "large"};
	printf("ts,U,M,n,size,HP,H,Omega\n");
	for (U = 0.2; U <= 0.9; U += 0.1) {
		for (M = 4; M <= 16; M *= 2) {
			{
				std::uniform_int_distribution<int> dist(2 * M, 150);
				n = dist(gen);
			}

			{
				if (f) fclose(f);

				char tmp[260];
				sprintf_s(tmp, "%.2f-%d.csv\0", U, M);
				if (fopen_s(&f, tmp, "w")) {
					printf("Failed to open file %s\n", tmp);
					return;
				}
				fprintf(f, "TS,U,M,n,Tsize,HP,H,Omega,OMLP,SMLP,K\n");
			}

			for( int j = 0; j < 3; ++j ) {
				tmin = (Tvalues.begin() + j)->first;
				tmax = (Tvalues.begin() + j)->second;
				tsize = *(Tnames.begin() + j);

				for (tsindex = 0; tsindex < 100; ++tsindex) {
					taskSet = make_ts(totalutil, Lmax, LHmax);
					//printf("Got task set: %d\n", taskSet->size());

					// calculate hyperperiod
					hyperperiod = std::min( std::accumulate(taskSet->begin(), taskSet->end(), 1, [](int a, task* b) { return std::lcm(a, (int)b->period); }), 5000 );
					//printf("hyperperiod: %d\n", hyperperiod);
					//printf("total util: %.2f (want %.2f)\n", totalutil, M * U);
					//printf("Lmax, LHmax: %.2f, %.2f\n", Lmax, LHmax);

					for (H = 4; H <= 64; H *= 2) {
						Omega = -1;
						//printf("Unsliced\n");
						//printf("TS %d, M %d, H %d, Omega %.2f                     \r", tsindex, M, H, Omega);
						printf("%d,%.2f,%d,%d,%s,%d,%d,%.2f                         \r", tsindex, U, M, n, tsize.c_str(), hyperperiod, H, Omega );
						runsim(hyperperiod, *taskSet);
						for (Omega = 0.5; Omega <= 2; Omega += 0.5) {
							if (Omega < LHmax) // cannot slice less than the highest request length
								continue;
							printf("%d,%.2f,%d,%d,%s,%d,%d,%.2f                         \r", tsindex, U, M, n, tsize.c_str(), hyperperiod, H, Omega);
							//printf("Sliced: %.2f\n", Omega);
							runsim(hyperperiod, *taskSet);
						}
					}

					for (auto& i : *taskSet) {
						if (i->r) {
							delete[] i->r->Li;
							delete i->r;
						}
						delete i;
					}
					delete taskSet;
				}
			}
		}
	}

	printf("Done.\n");
}

vector<task*>* make_ts(double& totalutil, double& Lmax, double& LHmax) {
	uniform_real_distribution<double> dist(0, 1);
	double* u = generate_utilizations();
	vector<task*>* taskSet = new vector<task*>();
	for (int i = 0; i < n; ++i) {
		auto Ti = sample_Ti();
		taskSet->push_back(new task{ std::min(M * u[i] * Ti, Ti), Ti, Ti, nullptr });
		totalutil += taskSet->at(i)->cost / taskSet->at(i)->period;
		if (dist(gen) < plock) {
			taskSet->at(i)->r = generate_request(taskSet->back());
			Lmax = std::accumulate(taskSet->at(i)->r->Li, taskSet->at(i)->r->Li + H + 1, Lmax, [](double a, double b) { return std::max(a, b); });
			LHmax = std::max(LHmax, taskSet->at(i)->r->Li[H]);
		}
	}
	delete[] u;

	return taskSet;
}

template <class T, class S, class C>
S& PQContainer(priority_queue<T, S, C>& q) {
	struct HackedQueue : private priority_queue<T, S, C> {
		static S& Container(priority_queue<T, S, C>& q) {
			return q.* & HackedQueue::c;
		}
	};
	return HackedQueue::Container(q);
}

int Z(request* r, double t, int It, zfunc& z, zminfunc& zmin);
double roundUpMult(double val, double mult);

double simulate(int hyperperiod, vector<task*>& taskSet, zfunc z, zminfunc zmin) {
	// simulate
	vector<job*> jobs;
	// releases can happen at the same time as slice, so make sure slice gets priority
	auto actcmp = [](act* a, act* b) { return abs(a->time - b->time) < 0.001 ? a->type != SLICE : a->time > b->time; };
	//auto actcmp = [](act* a, act* b) { return a->time > b->time; };
	priority_queue<act*, std::vector<act*>, decltype(actcmp)> events;

	for (auto i : taskSet) {
		for (double t = 0; t < hyperperiod; t += i->period) {
			jobs.push_back(new job{ t, i->cost, t + i->period, i->cost, -1, i->r });
			events.push(new act{ RELEASE, t, jobs.back() });
		}
	}

	if (Omega > 0) {
		for (double t = 0; t < hyperperiod; t += Omega)
			events.push(new act{ SLICE, t, nullptr });
	}

	int deadline_miss_count = 0;
	double closest_miss = -1;
	double worst_piblock = 0;

	job** cpu = new job * [M];
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
	deque<request*> reqFQ;
	int numInSQ = 0;

	auto& reqpq_vec = PQContainer(reqPQ);

	auto satisfy = [&It, &reqFQ, &reqPQ, &events, &numInSQ, &z, &zmin](double t) {
		if (reqFQ.size() + numInSQ < M && reqPQ.size()) {
			reqFQ.push_back(reqPQ.top());
			reqPQ.pop();
		}

		if (It < h)
			return false;

		request* r = nullptr;
		int fqpos = 0;
		int fq_erase = -1;
		vector<request*> requests_removedPQ;
		for(;;) {
			r = nullptr;

			// skip ahead mechanism (possibly)
			if (reqFQ.size() > fqpos) {
				fq_erase = fqpos;
				r = reqFQ[fqpos++];
			}
			else fq_erase = -1;

			// find the next request in PQ
			if (!r) {
				if (!reqPQ.size())
					break;
				r = reqPQ.top();
				requests_removedPQ.push_back(r);
				reqPQ.pop();
			}

			r->SMcount = Z(r, t, It, z, zmin);
			if (r->SMcount > 0)
				break;

			// only if we're specifically fz-blocked do we do the skipahead
			if (Omega < 0 || It < zmin(r)) {
				r = nullptr;
				break; // can't be fz-blocked
			}
			int zpossible = z(r, It);
			if (r->Li[zpossible] + t + slicemargin < roundUpMult(t, Omega)) {
				r = nullptr;
				break; // not fz-blocked
			}
		}

		for (auto& i : requests_removedPQ) {
			if( i != r )
				reqPQ.push(i);
		}

		if (!r) return false;

		if (fq_erase >= 0)
			reqFQ.erase(reqFQ.begin() + fq_erase);

		It -= r->SMcount;
		r->j->remCost -= r->Li[r->SMcount];
		events.push(new act{ UNLOCK, t + r->Li[r->SMcount], (job*)r });
		++numInSQ;
		/*if (Omega > 0) {
			if (t + r->Li[r->SMcount] > roundUpMult(t, Omega)) {
				printf("bad allocation at %.2f %.2f\n", t, roundUpMult(t, Omega));
			}
			else printf("Req at %.2f will complete %.2f (next slice %.2f)\n",
				t, t + r->Li[r->SMcount], Omega < 0 ? -1 : roundUpMult(t, Omega));
		}*/

		return true;
	};

	for (double t = 0; t < hyperperiod * 2; ) {
		//printf("Progress: %.2f%%\r", t / hyperperiod * 100);
		if (!events.size())
			break;
		act* a = events.top();
		auto tprev = t;
		bool nopop = false;

		int smallest_rem = -1;
		for(int i = 0; i < M; ++i) {
			if (cpu[i] && (smallest_rem < 0 || cpu[i]->remCost < cpu[smallest_rem]->remCost))
				smallest_rem = i;
		}

		if (smallest_rem > 0 && cpu[smallest_rem]->remCost < a->time - t) {
			a = new act{ COMPLETE, t + cpu[smallest_rem]->remCost, cpu[smallest_rem] };
			nopop = true;
		}

		if( !nopop )
			events.pop();
		t = a->time;

		// determine pi-blocking time in this duration
		// if a request is in the FQ/PQ and in the top M priorities, then it is pi-blocked
		for (auto& i : reqFQ) {
			auto lowest = lowestcpu();
			if (!cpu[lowest] || cpu[lowest]->deadline > i->j->deadline)
				i->piblockingtime += t - tprev;
		}
		for (auto& i : reqpq_vec) {
			auto lowest = lowestcpu();
			if (!cpu[lowest] || cpu[lowest]->deadline > i->j->deadline)
				i->piblockingtime += t - tprev;
		}

		//for (auto& i : cpu) {
			//if(i) i->remCost -= t - tprev;
		//}
		for(int i = 0; i < M; ++i)
			if( cpu[i] ) cpu[i]->remCost -= t - tprev;

		if (a->type == SLICE) {
			if (numInSQ > 0) {
				printf("SQ %d non-empty at slice: %.2f\n", numInSQ, a->time);

			}
			while (satisfy(t));
		}
		else if (a->type == RELEASE) {
			if (a->j->r) {
				auto r = new request;
				memcpy(r, a->j->r, sizeof(request));
				r->issuetime = t;
				r->j = a->j;
				r->piblockingtime = 0;
				reqPQ.push(r);
				satisfy(t);
			}
			else {
				int lowest = lowestcpu();
				if (!cpu[lowest] || cpu[lowest]->deadline > a->j->deadline) {
					cpu[lowest] = a->j;
					a->j->cpu = lowest;
				}
				else
					cpuPQ.push(a->j);
			}
		}
		else if (a->type == UNLOCK) {
			if( a->r->piblockingtime > worst_piblock )
				worst_piblock = a->r->piblockingtime;
			--numInSQ;
			It += a->r->SMcount;
			//if( Omega > 0 )
				//printf("%.2f: Req %.2f length responded in %.2f\n", t, a->r->Li[a->r->SMcount], t - a->r->issuetime);
			while (satisfy(t));

			int lowest = lowestcpu();
			if (!cpu[lowest] || cpu[lowest]->deadline > a->r->j->deadline) {
				cpu[lowest] = a->r->j;
				a->r->j->cpu = lowest;
			}
			else
				cpuPQ.push(a->r->j);
			delete a->r;
		}
		else if (a->type == COMPLETE) {
			if (t > a->j->deadline) {
				//printf("deadline missed at %.2f (dl=%.2f)\n", t, a->j->deadline);
				++deadline_miss_count;
			}
			else {
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
	//printf("\n");
	//printf("worst observed pi-blocking: %.2f\n", worst_piblock);
	//printf("number of deadlines missed: %d\n", deadline_miss_count);
	//printf("closest non-miss: %.2f\n", closest_miss);
	//printf("simulation complete\n");
	//printf("non-empty requests: %d %d\n", reqFQ.size(), reqPQ.size());

	// clean up
	for (auto i : jobs)
		delete i;
	while (events.size()) {
		auto a = events.top();
		events.pop();
		delete a;
	}

	return worst_piblock;
}

int Z(request* r, double t, int It, zfunc& z, zminfunc& zmin) {
	double nextSlice = Omega < 0 ? DBL_MAX : roundUpMult(t, Omega);

	for (int i = z(r, It); i >= zmin(r) && i <= It; i -= h) {
		if (Omega < 0 || r->Li[i] + t + slicemargin < nextSlice)
			return i;
	}

	return 0;
}

request* generate_request(task* t) {
	request* r = new request;
	r->Li = new double[H+1];
	r->SMcount = 0;
	r->j = nullptr;
	double maxRequestSize = h * t->cost / smmult;
	uniform_real_distribution<double> dist(minReqSizeMultiplier * maxRequestSize, maxRequestSize);
	r->reqsize = dist(gen);
	r->piblockingtime = 0;
	memset(r->Li, 0, sizeof(double) * (H+1));
	for (int i = h; i <= H; ++i) {
		r->Li[i] = smmult * ceil(r->reqsize / i) * (1 + (i - h) * paraloss);
		//printf("%.2f ", r->Li[i]);
	}
	//printf("\n");
	//uniform_real_distribution<double> dist2(0.5, 1);
	//r->k = z_smlp(r, (int)ceil(dist2(gen) * (double)H));
	r->k = z_smlp(r, H);
	//printf("value for k: %d\n", r->k);
	return r;
}

int z_smlp(request* r, int It) {
	auto cyc = [r](int x) {
		int roundval = (int)ceil(r->reqsize / x);
		return ((roundval + h - 1) / h) * h;
	};

	int smallest = It;
	for (int i = It; i > 0; i -= h) {
		if (cyc(i) <= cyc(smallest))
			smallest = i;
	}

	return smallest;
}

int z_omlp(request* r, int It) {
	return H;
}

int z_k(request* r, int It) {
	return r->k;
}

int zmin_omlp(request* r) {
	return H;
}

int zmin_smlp(request* r) {
	return h;
}

int zmin_k(request* r) {
	return r->k;
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

double roundUpMult(double val, double mult) {
	return std::floor( (val+mult) / mult ) * mult;
}