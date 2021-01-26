#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <atomic>
#include <mutex>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "get_time.h"
#include "utilities.h"


using namespace std;

const  double PI = 3.1415926;
const  double EPSINON = 0.0000001;

inline uint32_t hash1(uint32_t a) {
	a = (a + 0x7ed55d16) + (a << 12);
	a = (a ^ 0xc761c23c) ^ (a >> 19);
	a = (a + 0x165667b1) + (a << 5);
	a = (a + 0xd3a2646c) ^ (a << 9);
	a = (a + 0xfd7046c5) + (a << 3);
	a = (a ^ 0xb55a4f09) ^ (a >> 16);
	if (a < 0) a = -a;
	return a;
}

inline uint32_t hash2(uint32_t a) {
	a = (a + 0x1a976ce5) + (a << 15);
	a = (a ^ 0x98bf7c21) ^ (a >> 13);
	a = (a + 0xc092738a) + (a << 7);
	a = (a + 0xbbc2976e) ^ (a << 9);
	a = (a + 0xf6efd123) + (a << 11);
	a = (a ^ 0x987abc01) ^ (a >> 2);
	if (a < 0) a = -a;
	return a;
}

int hash3(uint32_t a, uint32_t b, int size) {
	a = (hash2(a) ^ hash2(b) << 1)% size;
	if (a < 0) a = -a;
	return (a);
}

struct point
{
	double x, y;
	int pivot, n;
	bool operator < (const point &o) const
	{
		return pivot < o.pivot;
	}
	bool operator == (const point &o) const
	{
		return pivot == o.pivot;
	}/**/
};

typedef struct facet
{
	point v1, v2;
	bool operator < (const facet &o) const
	{
		if (v1.pivot == o.v1.pivot)
			return v2.pivot < o.v2.pivot;
		return v1.pivot < o.v1.pivot;
	}
	bool operator == (const facet &o) const
	{
		return v1.pivot == o.v1.pivot && v2.pivot == o.v2.pivot;
	}
}Facet;

struct mapCkey
{
	int p1, p2;
};

struct mapCvalue
{
	point *it;
	int size;
};

facet* H;
point* p;
point* space;
bool* small;
bool* mid;
bool* large;
pair <int, int> *R;
pair <facet, mapCvalue> *hashC;
int* hashlock;
//point* pmerge;

int hashfacet(facet o) {
	int a = (hash1(o.v1.pivot) ^ hash1(o.v2.pivot) << 1) % (o.v1.n * 10);
	if (a < 0) a = -a;
	return (a);
}


bool mapCinsert(facet t, mapCvalue val) {
	int i = hashfacet(t);
	while (!pbbs::atomic_compare_and_swap(&hashlock[i], 0, 1)) {
		i = (i + 1) % (t.v1.n * 10);
	}
	hashC[i] = pair<facet, mapCvalue>(t, val);
	return true;
}

mapCvalue mapCfind(facet t) {
	int i = hashfacet(t);
	while (!(hashC[i].first == t)) {
		i = (i + 1) % (t.v1.n * 10);
	}
	mapCvalue out = hashC[i].second;
	return out;
}

bool InsertandSetR(facet t) {
	int i = hash1(t.v1.pivot) % (t.v1.n * 10);
	while (!pbbs::atomic_compare_and_swap(&R[i], pair <int, int>(0, 0), pair<int, int>(t.v1.pivot, t.v2.pivot))) {
		if (R[i].first == t.v1.pivot) {
			return false;
		}
		i = (i + 1) % (t.v1.n * 10);
	}
	return true;
}

bool InsertandSetL(facet t) {
	int i = hash1(t.v2.pivot) % (t.v2.n * 10);
	while (!pbbs::atomic_compare_and_swap(&R[i], pair <int, int>(0, 0), pair<int, int>(t.v2.pivot, t.v1.pivot))) {
		if (R[i].first == t.v2.pivot) {
			return false;
		}
		i = (i + 1) % (t.v2.n * 10);
	}
	return true;
}

facet getvalueL(facet t) {
	int i = hash1(t.v2.pivot) % (t.v2.n * 10);
	while (!(R[i].first == t.v2.pivot)) {
		i = (i + 1) % (t.v2.n * 10);
	}
	facet out = { p[R[i].first], p[R[i].second] };
	return out;
}

facet getvalueR(facet t) {
	int i = hash1(t.v1.pivot) % (t.v1.n * 10);
	while (!(R[i].first == t.v1.pivot)) {
		i = (i + 1) % (t.v1.n * 10);
	}
	facet out = { p[R[i].second], p[R[i].first] };
	return out;
}

void RandmOnCircle(point *p, int n) {
	for (int i = 0; i < n; i++) {
		p[i].x = (double)n + n * cos((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].y = (double)n + n * sin((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].pivot = i;
		p[i].n = n;
	}
	return;
}

bool visible(point v, facet t) {
	if (EPSINON < (t.v1.x - v.x) * (t.v2.y - v.y) - (t.v1.y - v.y) * (t.v2.x - v.x)) return true;
	else return false;
}

bool start(point v1, point v2, point v3) {
	if (EPSINON < (v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x)) return true;
	else return false;
}

bool collinear(point v1, point v2, point v3) {
	if ((v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) < EPSINON && (v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) > -EPSINON) return true;
	else return false;
}

int computevisblepoint1(point* A, facet t,int n) {
	if (n <= 256) {
		int ret = 0;
		for (int i = 0; i < n; i++) {
			if (visible(A[i], t)) ret++;
		}
		return ret;
	}
	int L, R;
	L = cilk_spawn computevisblepoint1(A, t, n / 2);
	R = computevisblepoint1(A + n / 2, t, n - n / 2);
	cilk_sync;
	return L + R;
}

void pointmerge(mapCvalue v1, mapCvalue v2, facet t) {
	int pointer0 = 0;
	if (v1.size == 0) {
		pointer0 = computevisblepoint1(v2.it,t,v2.size);
	}
	else if (v2.size == 0) {
		pointer0 = computevisblepoint1(v1.it, t, v1.size);
	}
	else {
		int pointer1 = 0;
		int pointer2 = 0;
		while (pointer1 < v1.size && pointer2 < v2.size)
		{
			if (v1.it[pointer1] < v2.it[pointer2]) {
				if (visible(v1.it[pointer1], t)) {
					pointer0++;
				}
				pointer1++;
			}
			else if (v2.it[pointer2] < v1.it[pointer1]) {
				if (visible(v2.it[pointer2], t)) {
					pointer0++;
				}
				pointer2++;
			}
			else {
				if (visible(v1.it[pointer1], t)) {
					pointer0++;
				}
				pointer1++;
				pointer2++;
			}
		}
		if (pointer1 != v1.size) {
			for (int i = pointer1; i < v1.size; i++) {
				if (visible(v1.it[i], t)) {
					pointer0++;
				}
			}
		}
		else if (pointer2 != v2.size) {
			for (int i = pointer2; i < v2.size; i++) {
				if (visible(v2.it[i], t)) {
					pointer0++;
				}
			}
		}
	}
	point *pmerge = new point[pointer0];//CHalloc(pointer0,hash3(t.v1.pivot,t.v2.pivot, t.v1.n * 32),t.v1.n*32);
	mapCvalue tempval;
	tempval.it = pmerge;
	pointer0 = 0;
	if (v1.size == 0) {
		for (int i = 0; i < v2.size; i++) {
			if (visible(v2.it[i], t)) {
				pmerge[pointer0] = v2.it[i];
				pointer0++;
			}
		}
	}
	else if (v2.size == 0) {
		for (int i = 0; i < v1.size; i++) {
			if (visible(v1.it[i], t)) {
				pmerge[pointer0] = v1.it[i];
				pointer0++;
			}
		}
	}
	else {
		int pointer1 = 0;
		int pointer2 = 0;
		while (pointer1 < v1.size && pointer2 < v2.size)
		{
			if (v1.it[pointer1] < v2.it[pointer2]) {
				if (visible(v1.it[pointer1], t)) {
					pmerge[pointer0] = v1.it[pointer1];
					pointer0++;
				}
				pointer1++;
			}
			else if (v2.it[pointer2] < v1.it[pointer1]) {
				if (visible(v2.it[pointer2], t)) {
					pmerge[pointer0] = v2.it[pointer2];
					pointer0++;
				}
				pointer2++;
			}
			else {
				if (visible(v1.it[pointer1], t)) {
					pmerge[pointer0] = v1.it[pointer1];
					pointer0++;
				}
				pointer1++;
				pointer2++;
			}
		}
		if (pointer1 != v1.size) {
			for (int i = pointer1; i < v1.size; i++) {
				if (visible(v1.it[i], t)) {
					pmerge[pointer0] = v1.it[i];
					pointer0++;
				}
			}
		}
		else if (pointer2 != v2.size) {
			for (int i = pointer2; i < v2.size; i++) {
				if (visible(v2.it[i], t)) {
					pmerge[pointer0] = v2.it[i];
					pointer0++;
				}
			}
		}
	}
	
	tempval.size = pointer0;
	mapCinsert(t, tempval);
	return;
}

void prefixsum(int* In, int* Out, int* B, int* C, int n) {
	if (n == 0) return;
	if (n == 1) {
		Out[0] = In[0];
		return;
	}
	cilk_for(int i = 0; i < n / 2; i++)
		B[i] = In[2 * i] + In[2 * i + 1];

	prefixsum(B, C, B + n / 2, C + n / 2, n / 2);
	Out[0] = In[0];

	cilk_for(int i = 1; i < n; i++) {
		if (i % 2) Out[i] = C[i / 2];
		else Out[i] = C[i / 2 - 1] + In[i];
	}
}

void pcopy(point *s, int size, point *news){
	for (int i = 0; i < size; i++) {
		news[i] = s[i];
	}
}

void parallelmerge(mapCvalue v1, mapCvalue v2, facet t, facet t2) {
	int n1 = v1.size / 1000 + 1;
	int n2 = v2.size / 1000 + 1;
	int *count1 = new int[n1];
	int *count2 = new int[n2];
	int *count1o = new int[n1];
	int *count2o = new int[n2];
	int *count1b = new int[n1];
	int *count2b = new int[n2];
	int *count1c = new int[n1];
	int *count2c = new int[n2];
	
	count1[0] = 0;
	for (int i = 0; i < v1.size; i += 1000) {
		int number = 0;
		for (int j = i; j < i + 1000 && j < v1.size; j++) {
			if (visible(v1.it[j], t) && visible(v1.it[j], t2)) {
				point temp = v1.it[j];
				v1.it[j] = v1.it[i + number];
				v1.it[i + number] = temp;
				number++;
			}
		}
		count1[i / 1000] = number;
	}
	count2[0] = 0;
	for (int i = 0; i < v2.size; i += 1000) {
		int number = 0;
		for (int j = i; j < i + 1000 && j < v2.size; j++) {
			if (visible(v2.it[j], t)) {
				point temp = v2.it[j];
				v2.it[j] = v2.it[i + number];
				v2.it[i + number] = temp;
				number++;
			}
		}
		count2[i / 1000] = number;
	}
	prefixsum(count2, count2o, count2b, count2c, n2);
	prefixsum(count1, count1o, count1b, count1c, n1);
	int psize = count1[sizeof(count1) / 4 - 1] + count2[sizeof(count2) / 4 - 1];
	point *pmerge = new point[psize];
	for (int i = 0; i < sizeof(count1) / 4; i++) {
		if (i == 0) {
			pcopy(v1.it, count1o[0], pmerge);
		}
		else {
			pcopy(v1.it + i * 1000, i * 1000 + count1o[i] - count1o[i - 1], pmerge + count1o[i - 1]);
		}
	}
	for (int i = 0; i < sizeof(count2) / 4; i++) {
		if (i == 0) {
			pcopy(v2.it, count2o[0], pmerge + count1o[sizeof(count1) / 4]);
		}
		else {
			pcopy(v2.it + i * 1000, i * 1000 + count2o[i] - count2o[i - 1], pmerge + count1o[sizeof(count1) / 4] + count2o[i - 1]);
		}
	}
}

void ProcessRidge(facet t1, point r, facet t2) {
	mapCvalue t1val = mapCfind(t1);
	mapCvalue t2val = mapCfind(t2);
	if ((t1val.size == 0) && (t2val.size == 0)) {
		H[t1.v1.pivot] = t1;
	}
	else if (!(t2val.size == 0) && !(t1val.size == 0) && t1val.it[0] == t2val.it[0]) {
	}
	else if ((t1val.size == 0) || (t2val.size != 0 && (t2val.it[0] < t1val.it[0]))) {//左
		facet t;
		t.v1 = r;
		t.v2 = t2val.it[0];
		//pointmerge(t1val, t2val, t);
		parallelmerge(t1val, t2val, t, t2);
		cilk_for(int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t1, r, t);
			}
			else if (!InsertandSetL(t)) {
				facet oldt = getvalueL(t);
				ProcessRidge(t, t.v2, oldt);
			}
		}
	}
	else {//右
		facet t;
		t.v1 = t1val.it[0];
		t.v2 = r;
		//pointmerge(t1val, t2val, t);
		parallelmerge(t1val, t2val, t, t2);
		cilk_for(int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t, r, t2);
			}
			else if (!InsertandSetR(t)) {
				facet oldt = getvalueR(t);
				ProcessRidge(oldt, t.v1, t);
			}
		}
	}
	return;
}

void initpoints(int in, int n) {
	if (in == 0) {
		for (int i = 0; i < n; i++) {
			p[i].x = (hash1(i)) % (n * 2);
			p[i].y = (hash2(i)) % (n * 2);
			p[i].pivot = i;
			p[i].n = n;
		}
	}
	else if (in == 1) {
		RandmOnCircle(p, n);
	}
	else 
		cout << "Usage: ./qsort [num_elements] [0 or 1]" << endl;

}

void init(int n) {
	mapCvalue tempval;
	facet t;
	if (start(p[0], p[1], p[2])) {
		int initpointer = 0;
		t = { p[0],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep1 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep1[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep1;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep2 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep2[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep2;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep3 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep3[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep3;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

	}
	else {
		int initpointer = 0;
		t = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep1 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep1[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep1;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep2 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep2[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep2;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep3 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep3[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep3;
		tempval.size = initpointer;
		mapCinsert(t, tempval);
	}
}

int main(int argc, char** argv) {

	if (argc != 3) {
		cout << "Usage: ./qsort [num_elements] [type of input]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	int type_of_input = atoi(argv[2]);
	p = new point[n];
	space = new point[96*n];
	small = new bool[32 * n];
	mid = new bool[32 * n / 100];
	large = new bool[32 * n / 10000];
	R = new pair<int, int>[10 * n];
	H = new facet[n];
	hashC = new pair<facet, mapCvalue>[10 * n];
	hashlock = new int[10 * n];
	timer timer1;
	facet facet11, facet12, facet21, facet22, facet31, facet32;
	initpoints(type_of_input, n);
	timer1.start();
	
	init(n);
	if (start(p[0], p[1], p[2])) {
		facet11 = { p[1],p[0] }; facet12 = { p[0],p[2] };
		facet21 = { p[2],p[1] }; facet22 = { p[1],p[0] };
		facet31 = { p[0],p[2] }; facet32 = { p[2],p[1] };
		cilk_for(int i = 0; i < 3; i++) {
			if (i == 0) {
				ProcessRidge(facet11, p[0], facet12);
			}
			if (i == 1) {
				ProcessRidge(facet21, p[1], facet22);
			}
			if (i == 2) {
				ProcessRidge(facet31, p[2], facet32);
			}
		}
	}
	else {
		facet11 = { p[2],p[0] }; facet12 = { p[0],p[1] };
		facet21 = { p[0],p[1] }; facet22 = { p[1],p[2] };
		facet31 = { p[1],p[2] }; facet32 = { p[2],p[0] };
		cilk_for(int i = 0; i < 3; i++) {
			if (i == 0) {
				ProcessRidge(facet11, p[0], facet12);
			}
			if (i == 1) {
				ProcessRidge(facet21, p[1], facet22);
			}
			if (i == 2) {
				ProcessRidge(facet31, p[2], facet32);
			}
		}
	}
	timer1.stop();
	int Hsize = 0;
	for (int i = 0; i < n; i++) {
		if (H[i].v1.n != 0) {
			cout << H[i].v1.pivot << "," << H[i].v2.pivot << endl;
			Hsize++;
		}

	}

	cout << "total time: " << timer1.get_total() << endl;
	cout << "output size: " << Hsize << endl;
	return 0;
}
